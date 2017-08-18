import os
import sys
import re
import argparse

import numpy as np
import scipy
from scipy.ndimage import imread
import skimage.exposure as skime

import proj_lib as proj
import projpp_lib as projpp


def flush_message(message):
    print(message),
    sys.stdout.flush()


def load_image_stack(path):
    tiff_files = [f for f in os.listdir(path) if f.endswith("tiff")]
    im = imread(os.path.join(path, "stack0.tiff"))

    shape = (im.shape[0], im.shape[1], len(tiff_files), 3)

    stack = np.zeros(shape, dtype=np.uint8)

    if not os.path.isfile(os.path.join(path, "order.txt")):
        for tiff_file in tiff_files:
            z_id = int(re.findall('\d+', str(tiff_file))[0])
            stack[:, :, z_id, :] = imread(os.path.join(path, tiff_file))
    else:
        for tiff_file in tiff_files:
            z_id = int(re.findall('\d+', str(tiff_file))[0])
            stack[:, :, len(tiff_files) - z_id - 1, :] = imread(os.path.join(path, tiff_file))

    return stack


def gauss_surface(image_stack, sdx, sdy, sdz, sds):
    """blurs the image stack and returns the projection surface"""
    flush_message("Finding Surface (Gaussian Blur method)... ")
    blurred_image = scipy.ndimage.gaussian_filter(image_stack, [sdx, sdy, sdz])
    projection_surface = proj.max_indices_z(blurred_image, clip_bottom=True)
    smoothed_projection_surface = scipy.ndimage.gaussian_filter(projection_surface, sds)
    print("Done")
    return smoothed_projection_surface


def main():

    sdx, sdy, sdz = 6, 6, 5
    sds = 6

    exp_dir = os.path.abspath(args.exp_dir)
    stack_dir = os.path.join(exp_dir, "stack")
    output_dir = exp_dir
    outname = "_".join(os.path.basename(exp_dir).split('_')[1:-1])

    image3d = load_image_stack(stack_dir)
    greyscale_image_stack = np.amax(image3d, 3)

    max_proj = np.amax(greyscale_image_stack, axis=2)

    mpfilename = os.path.join(output_dir, "%s_max-proj.png" % outname)
    scipy.misc.imsave(mpfilename, max_proj)

    sps = gauss_surface(greyscale_image_stack, sdx, sdy, sdz, sds)

    flush_message("Projecting... ")
    # vis_factor = 255 / z_size
    sfilename = os.path.join(output_dir, "%s_surface-g3d.png" % outname)
    scipy.misc.imsave(sfilename, sps)
    res = proj.projection_from_surface_z(greyscale_image_stack, sps, dm=3, dp=0)
    # res_scale = 255 / np.amax(res)
    filename = os.path.join(output_dir, "%s_proj-g3d.png" % outname)
    scipy.misc.imsave(filename, res)
    print("Done")

    flush_message("Projecting... ")
    # vis_factor = 255 / z_size
    res_rev = proj.projection_from_surface_z(greyscale_image_stack, sps, dm=0, dp=3)
    # res_scale = 255 / np.amax(res)
    filename = os.path.join(output_dir, "%s_proj-g3d_rev.png" % outname)
    scipy.misc.imsave(filename, res_rev)
    print("Done")

    # flush_message("Calculating angle...")
    # norm_angle = proj.normal_angle_from_surface(greyscale_image_stack, sps)
    # norm_angle = float(255) * (norm_angle/float(90))
    # filename = os.path.join(output_dir, "%s_angle.png" % outname)
    # scipy.misc.imsave(filename, norm_angle)
    # print("Done")

    # flush_message("Post processing...")
    # post_processed_image = projpp.proj_filter(res, 3, 60, 15)
    # filename = os.path.join(output_dir, '%s_proj-pp-g3d.png' % outname)
    # # pp_scale = 255 / np.amax(post_processed_image)
    # scipy.misc.imsave(filename, post_processed_image)
    # print("Done")

    flush_message("Equalizing... ")
    filename = os.path.join(output_dir, "%s_proj-pp-clahe-g3d.png" % outname)
    clahe = skime.equalize_adapthist(post_processed_image, clip_limit=0.01)
    # clahe_scale = 255 / np.amax(clahe)
    scipy.misc.imsave(filename, clahe)
    print("Done")
    print("------------")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('exp_dir', help="path to directory containing image stack")

    args = parser.parse_args()

    # load_image_stack(os.path.join(args.exp_dir, "stack"))

    main()
