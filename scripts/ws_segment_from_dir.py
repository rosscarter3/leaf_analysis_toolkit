import os
import argparse

import numpy as np
import skimage.morphology
import scipy.ndimage as nd
import matplotlib.pyplot as plt

import common_functions as cf


def get_im_paths(exp_dir):
    im_path = []
    for dirpath, _, filenames in os.walk(exp_dir):
        for f in filenames:
            if f.endswith('_proj-g3d.png') or f.endswith('_proj-g3d_rev.png'):
                im_path.append(os.path.abspath(os.path.join(dirpath, f)))
    return im_path


def get_seeds_path(exp_dir):
    for dirpath, _, filenames in os.walk(exp_dir):
        for f in filenames:
            if "seeds" in f:
                seeds_path = os.path.abspath(os.path.join(dirpath, f))
    return seeds_path


def watershed(im_path, seeds_path):
    im = nd.imread(im_path)
    seeds = nd.imread(seeds_path)

    seed_array = np.array(seeds)

    seed_array_bool = np.ones_like(im)
    for i in xrange(seed_array.shape[0]):
        for j in xrange(seed_array.shape[1]):
            if seed_array[i, j, 0] == seed_array[i, j, 1] == seed_array[i, j, 2]:
                seed_array_bool[i, j] = 0

    seed_array_bool = skimage.morphology.label(seed_array_bool)
    seg = skimage.morphology.watershed(im, seed_array_bool)
    seg[np.where(seg == seg[0, 0])] = 0

    seg = cf.id_array2rgb(seg)

    seg_path = os.path.join(im_path + "ws_seg.png")

    plt.imsave(seg_path, seg)

    return 0


def main():
    exp_dir = args.dir_path

    im_paths = get_im_paths(exp_dir)
    seeds_path = get_seeds_path(exp_dir)

    for im_path in im_paths:
        print "Segmenting: " + os.path.basename(im_path)
        watershed(im_path, seeds_path)

    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("dir_path", help="")
    args = parser.parse_args()
    main()
