#!/usr/bin/env python2.7
"""Methods for simple projection in the z-direction"""
import sys
from math import sqrt
import scipy
import scipy.ndimage as nd
import numpy as np
import numpy.linalg as la


def flush_message(message):
    print(message),
    sys.stdout.flush()


def max_indices_z(A, clip_bottom=True):
    """
    From a 3D matrix representing an image stack,
    find the depth of the voxel with maximum intensity in the z-direction.

    Args:
        A (numpy.ndarray): 3D volumetric image data

        clip_bottom (Optional[bool]): For columns which have uniform signal,
                                      whether to return the index of the top or
                                      bottom (default) of the image stack

    Returns:
        numpy.ndarray[int]: 2D array containing the depth of the projection surface
    """

    h = np.argmax(A, 2)
    if clip_bottom:
        h[np.logical_and(h == 0, A[:, :, 0] == A[:, :, -1])] = A.shape[2] - 1
    return h


def threshold_indices_z(A, b, c, clip_bottom=True):
    """
    From a 3D matrix representing an image stack,
    find the depth of the first voxel which exceeds a threshold intensity
    (depending on the mean intensity of the pixels in this "pillar")

    Args:
        A (numpy.ndarray): 3D volumetric image data
        b, c (float): Threshold is taken to be b*mean of the pillar + c

        clip_bottom (Optional[bool]): For columns which have uniform signal,
                                      whether to return the index of the top or
                                      bottom (default) of the image stack

    Returns:
        numpy.ndarray[int]: 2D array containing the depth of the projection surface
    """

    t = b * np.mean(A, axis=2) + c
    # now wish to find first occurrence
    return max_indices_z(A > t[:, :, np.newaxis], clip_bottom)


def sheet_indices_z(A, niter=1000, dt=0.05, a=3.0, b=3.0):
    """
    From a 3D matrix representing an image stack, find the projection
    surface using an active contour. The projection surface starts at the
    smallest z-depth (h=0), and moves in the direction of increasing z;
    this downwards speed is reduced in the presence of signal. A smoothing
    term is also applied to the surface, which encourages nearby points
    to be at similar depths.

    Args:
        A (numpy.ndarray): 3D volumetric image data
        niter (int): Number of timesteps for the evolution of the
                     projection surface
        dt (float): Timestep for the surface evolution.
                    May need to be limited to maintain stability
                    of the numerical method.
        a (float): Smoothing parameter - larger values penalize differences
                   in the depth of adjacent points.
        b (float): Parameter controlling how much the image intensity slows the
                   downwards movement of the projection surface; larger values
                   result in the image intensity slowing the projection
                   surface more.

    Returns:
        numpy.ndarray[float]: 2D array containing the depth of the projection surface
    """

    h0 = np.zeros(A.shape[0:2])
    i, j = np.meshgrid(np.arange(A.shape[0]), np.arange(A.shape[1]), indexing='ij')
    d2 = np.zeros(h0.shape)
    for k in range(niter):
        d2[1:-1, 1:-1] = h0[:-2, 1:-1] + h0[2:, 1:-1] + h0[1:-1, :-2] \
                         + h0[1:-1, 2:] - 4 * h0[1:-1, 1:-1]
        A_h = A[i, j, h0.astype(int)]
        h0 += dt / (1 + b * A_h) + d2 * a * dt
        h0 = np.clip(h0, 0, A.shape[2] - 1)
    return h0


def projection_from_surface_z(ma, surface, dm=0, dp=10, op=np.max):
    """
    Project the 3D image stack onto a surface (orthogonally, in the
    z-direction. This involves applying an operator (maximum, or mean)
    to the pixels in a band (in the z-direction) about the projection surface.

    Args:
        ma: 3D image stack.
        surface: 2D indices of the projection surface.
        dm: Number of pixels considered "above" (negative z) the surface
        dp: Number of pixels (-1) considered "below" the surface
        op: Operator used to calculate the intensity of the projected
            image im[i,j] from the band of pixels
            ma[i, j, surface-dm:surface+dp]

    Returns:
        numpy.ndarray: Projected image
    """
    xmax, ymax, zmax = ma.shape
    res = np.zeros([xmax, ymax], dtype=ma.dtype)
    z0 = np.clip(surface - dm, 0, zmax - 1).astype(int)
    z1 = np.clip(surface + dp, 1, zmax).astype(int)
    for x in range(0, xmax):
        for y in range(0, ymax):
            res[x, y] = op(ma[x, y, z0[x, y]:z1[x, y]])
    return res


def projection_from_surface_normal_z(ma, surface, dm=0, dp=10,
                                     op=np.max, gradient_radius=0.5):
    """
    Project the 3D image stack, examining a ray of pixels *normal* to
    the projection surface. For each point in 2D, we consider a single ray
    passing through this point, normal to the surface. Points are sampled
    along this ray (from a slightly blurred version of the 3D stack), and
    an operator (default maximum) applied to these samples to give the
    intensity of the 2D image.

    Args:
        ma: 3D image stack.
        surface: depth of the projection surface (2D array)
        dm: Number of pixels considered "above" (negative z) the surface
        dp: Number of pixels (-1) considered "below" the surface
        op: Operator used to calculate the intensity of the projected
            image im[i,j] from the band of pixels
            ma[i, j, surface-dm:surface+dp]
        gradient_radius: Radius of the gaussian operator used to find the slope
                         of the projection surface.

    Returns:
        numpy.ndarray: Projected image
    """

    imax, jmax, kmax = ma.shape
    res = np.zeros([imax, jmax], dtype=ma.dtype)
    surface_di = scipy.ndimage.gaussian_filter1d(surface, 0.5, axis=0, order=1)
    surface_dj = scipy.ndimage.gaussian_filter1d(surface, 0.5, axis=1, order=1)

    def clip(i_, j_, k_):
        return np.array((max(0, min(i_, imax - 1)),
                         max(0, min(j_, jmax - 1)),
                         max(0, min(k_, kmax - 1))))

    for i in range(0, imax):
        for j in range(0, jmax):
            si = -surface_di[i, j]
            sj = -surface_dj[i, j]
            k = surface[i, j]
            sk = 1.0
            h = sqrt(si * si + sj * sj + sk * sk)
            si /= h
            sj /= h
            sk /= h
            p_start = clip(i - dm * si, j - dm * sj, k - dm * sk)
            p_end = clip(i + dp * si, j + dp * sj, k + dp * sk)
            ns = int(la.norm(p_end - p_start) + 1)
            tt = np.linspace(0, 1, ns)
            res[i, j] = op(scipy.ndimage.map_coordinates(ma,
                                                         ([(p_start[0] * (1 - t) + p_end[0] * t) for t in tt],
                                                          [(p_start[1] * (1 - t) + p_end[1] * t) for t in tt],
                                                          [(p_start[2] * (1 - t) + p_end[2] * t) for t in tt]),
                                                         order=1))
    return res
