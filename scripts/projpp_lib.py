#!/usr/bin/env python2.7
"""projection post processing tools"""
import os
import sys
import math
import scipy.ndimage as nd
# import Image
import numpy as np


# def numpy_draw_pil(r):
#    i = Image.fromarray(r)
#    i.show()


def print_transfer_fn(function):
    """prints the transfer function fro 0 to 255 in steps of ten"""
    for i in range(0, 255, 10):
        print(("%d: %d" % (i, function[i])))


def gen_transfer_fn():
    """returns a list of values for the transfer step function from 0 to 256"""
    values = {}
    threshold = 75
    for i in range(0, 256):
        if i < threshold:
            values[i] = 0
        else:
            values[i] = 255

    return values


def sigmoidal_trans(midp, scale_factor):
    """retuens a list of values for the sigmoidal transfer function from 0 to 256"""
    values = {}
    for i in range(0, 256):
        values[i] = 255 - 255 / (1 + math.exp((i - midp) / scale_factor))
    return values


def np_iter(numpy_array):
    """numpy iterator for images?"""
    xmax, ymax = numpy_array.shape

    return ((x, y) for x in range(0, xmax) for y in range(0, ymax))


def np_e(numpy_array):
    """return empty array the same size as numpy_array"""
    xmax, ymax = numpy_array.shape
    return np.empty([xmax, ymax], dtype=np.uint8)


def apply_transfer_fn(numpy_array, function):
    """applies the transfer function to the numpy array """
    numpy_iterator = np_iter(numpy_array)

    numpy_out = np_e(numpy_array)

    # TODO - nice fp way of doing this with real functions
    for x_coord, y_coord in numpy_iterator:
        numpy_out[x_coord, y_coord] = function[numpy_array[x_coord, y_coord]]

    return numpy_out


def proj_filter(raw, smooth, mid, steep):
    """applies median filter and sigmoidal threshold to raw image"""
    medf = nd.median_filter(raw, smooth)
    fns = sigmoidal_trans(mid, steep)
    thresh = apply_transfer_fn(medf, fns)
    return thresh


def main():
    """test stuff"""
    # filename = 'output/surface-g3d-30-30-7.png'
    try:
        filename = sys.argv[1]
    except IndexError:
        print(("Usage: %s filename" % (os.path.basename(sys.argv[0]))))
        sys.exit(1)

        # raw = scipy.misc.imread(filename)
        # numpy_draw_pil(raw)

        # print_transfer_fn(fns)
        # sys.exit(0)

        # fn = gen_transfer_fn()
        # print_transfer_fn(fn)

        # medf = nd.median_filter(raw, 3)
        # numpy_draw_pil(medf)

        # fns = sigmoidal_trans(60, 15)
        # thresh = apply_transfer_fn(medf, fns)
        # thresh = proj_filter(raw, 3, 60, 15)
        # numpy_draw_pil(thresh)

        # alpha = 1
        # fl = nd.gaussian_filter(raw, 1.0)
        # sharpened = thresh + alpha * (thresh - fl)
        # numpy_draw_pil(sharpened)

        # alpha = 1
        # fl = nd.gaussian_filter(thresh, 3.0)
        # sharpened = thresh + alpha * (thresh - fl)
        # numpy_draw_pil(sharpened)


if __name__ == "__main__":
    main()
