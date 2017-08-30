import os
import argparse

import numpy as np
import skimage.morphology
import scipy.ndimage as nd
import matplotlib.pyplot as plt
from PIL import Image, ImageFile

ImageFile.LOAD_TRUNCATED_IMAGES = True
import common_functions as cf


def get_im_paths(exp_dir):
    im_path = []
    for dirpath, _, filenames in os.walk(exp_dir):
        for f in filenames:
            # if f.endswith('_proj-g3d.png') or f.endswith('_proj-g3d_rev.png'):
            #     im_path.append(os.path.abspath(os.path.join(dirpath, f)))
            if f.endswith('_proj-g3d_rev.png') or f.endswith('_proj_g3d_rev.png'):
                im_path.append(os.path.abspath(os.path.join(dirpath, f)))
    print "got image paths!"
    return im_path


def get_seeds_path(exp_dir):
    seeds_path = None
    for dirpath, _, filenames in os.walk(exp_dir):
        for f in filenames:
            if "manual" in f:
                seeds_path = os.path.abspath(os.path.join(dirpath, f))
    print "Got seeds path!"
    return seeds_path


def watershed(im_path, seeds_path):
    im = Image.open(im_path)
    im = np.array(im)
    if len(im.shape) == 3:
        im = np.mean(im, axis=2)

    print im

    seeds = nd.imread(seeds_path)

    seed_array = seeds

    seed_array[:, :, 1] = 0
    seed_array[:, :, 2] = 0

    # TODO fix shape of seeds

    seed_array_bool = seed_array[:, :, 0] > 254

    seed_array_bool = skimage.morphology.label(seed_array_bool)
    seed_array_bool = skimage.morphology.remove_small_objects(seed_array_bool, 9)

    seg = skimage.morphology.watershed(im, seed_array_bool)

    # TODO remove small cells!!

    seg[np.where(seg == seg[0, 0])] = 0

    seg_col = seg
    color_path = os.path.join(im_path + "colorful.png")
    # rand_col = cf.rand_cmap(len(np.unique(seg_col)), verbose=False)
    plt.imshow(im, cmap='gray_r')
    plt.imshow(seg_col, alpha=0.6, cmap='prism')
    plt.subplots_adjust(left=0.04, bottom=0.01, right=0.9, top=0.96, wspace=0.2, hspace=0.2)
    plt.savefig(color_path, dpi=400)
    # plt.show()

    seg = cf.id_array2rgb(seg)
    seg_path = os.path.join(im_path + "ws_seg.png")
    plt.imsave(seg_path, seg)

    return 0


def main(args):
    print "Script running!"
    exp_dir = args.dir_path

    print "experiment directory: ", exp_dir

    im_paths = get_im_paths(exp_dir)
    seeds_path = get_seeds_path(exp_dir)

    print "im path is: ", im_paths
    print "seeds path is: ", seeds_path

    if seeds_path is None:
        print "No manual seeds image found\n"
        return

    for im_path in im_paths:
        print "Segmenting: " + os.path.basename(im_path)
        watershed(im_path, seeds_path)

    print "finished"
    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("dir_path", help="")
    args = parser.parse_args()
    main(args)
