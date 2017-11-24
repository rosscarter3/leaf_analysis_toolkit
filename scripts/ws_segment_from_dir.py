#! /usr/bin/python
import os
import argparse

import numpy as np
import skimage.morphology
import matplotlib.pyplot as plt
from PIL import Image, ImageFile

import common_functions as cf

ImageFile.LOAD_TRUNCATED_IMAGES = True


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

    seeds = Image.open(seeds_path)
    seed_array = np.array(seeds, dtype=np.int64)

    if seed_array.shape[2] == 3:
        red_condition = [255, 0, 0]
    elif seed_array.shape[2] == 4:
        red_condition = [255, 0, 0, 255]

    seed_array_bool = np.ones([seed_array.shape[0], seed_array.shape[1]])

    # plt.imshow(seed_array)
    # plt.show()

    for column in range(seed_array.shape[0]):
        for row in range(seed_array.shape[1]):
            if list(seed_array[column, row]) != red_condition:
                seed_array_bool[column, row] = 0

    plt.imshow(seed_array_bool)
    plt.show()

    seed_array_bool = skimage.morphology.label(seed_array_bool)

    # seed_array_bool = skimage.morphology.remove_small_objects(seed_array_bool, 2)

    seg = skimage.morphology.watershed(im, seed_array_bool, mask=np.ones_like(im))
    seg = skimage.morphology.remove_small_objects(seg, 20, in_place=True)

    seg[np.where(seg == seg[0, 0])] = 0

    seg_col = seg
    color_path = os.path.join(im_path + "colorful.png")
    rand_col = cf.rand_cmap(len(np.unique(seg_col)), verbose=False)
    plt.imshow(im, cmap='gray_r')
    plt.imshow(seg_col, alpha=0.6,cmap=rand_col)
    plt.imshow(seeds, alpha = 0.6)
    plt.subplots_adjust(left=0.04, bottom=0.01, right=0.9, top=0.96, wspace=0.2, hspace=0.2)
    # plt.show()
    plt.savefig(color_path, dpi=400)

    seg = cf.id_array2rgb(seg)
    seg_path = os.path.join(im_path + "ws_seg.png")
    plt.imsave(seg_path, seg)

    return 0


def auto_watershed(im_path):
    from scipy import ndimage as ndi
    from skimage.feature import peak_local_max
    from skimage.util import invert
    from skimage.filters import threshold_otsu, threshold_local
    from skimage.feature import peak_local_max

    im = Image.open(im_path)
    im = np.array(im)
    if len(im.shape) == 3:
        im = np.mean(im, axis=2)
    
    inv = invert(im)

    plt.imshow(inv)
    plt.show()

    otsu = inv > threshold_otsu(inv)

    im = skimage.morphology.binary_closing(otsu, selem=skimage.morphology.diamond(1))
    im = skimage.morphology.binary_erosion(im, selem=None, out=None)

    plt.imshow(im)
    plt.show()

    label_otsu =

    # block_size = 35
    # adaptive_thresh = threshold_local(inv, block_size, offset=10)
    # binary_adaptive = inv > adaptive_thresh
    
    # distance = ndi.distance_transform_edt(binary_adaptive)
    # local_maxi = peak_local_max(invert(distance), indices=False, footprint=np.ones((3, 3)),
    #                             labels=binary_adaptive)
    # plt.hist(distance)
    # plt.show()

    # dist_threshold = threshold_otsu(distance)
    # dist_threshold = threshold_local(distance, block_size, offset=10)
    # dist_otsu = distance > dist_threshold
    # markers = ndi.label(dist_otsu)[0]
    # labels = skimage.morphology.watershed(im, markers)
    # seg = labels
    #
    # plt.subplot(131)
    # plt.imshow(distance)
    # plt.subplot(132)
    # plt.imshow(markers)
    # plt.subplot(133)
    # plt.imshow(labels)
    # plt.show()
    #
    # seg[np.where(seg == seg[0, 0])] = 0
    #
    # seg_col = seg
    # color_path = os.path.join(im_path + "colorful.png")
    # rand_col = cf.rand_cmap(len(np.unique(seg_col)), verbose=False)
    # plt.imshow(im, cmap='gray_r')
    # plt.imshow(seg_col, alpha=0.6,cmap=rand_col)
    # plt.imshow(seeds, alpha = 0.6)
    # plt.subplots_adjust(left=0.04, bottom=0.01, right=0.9, top=0.96, wspace=0.2, hspace=0.2)
    # plt.show()
    # # plt.savefig(color_path, dpi=400)
    #
    # seg = cf.id_array2rgb(seg)
    # seg_path = os.path.join(im_path + "ws_seg.png")

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
        print "No manual seeds image found, attempting automating seeding\n"
        auto_watershed(im_paths[0])
        return

    # auto_watershed(im_paths[0])

    for im_path in im_paths:
        print "Segmenting: " + os.path.basename(im_path)
        watershed(im_path, seeds_path)

    print "finished"
    return 0


def watershed_test():
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import ndimage as ndi

    from skimage.morphology import watershed
    from skimage.feature import peak_local_max

    # Generate an initial image with two overlapping circles
    x, y = np.indices((80, 80))
    x1, y1, x2, y2 = 28, 28, 44, 52
    r1, r2 = 16, 20
    mask_circle1 = (x - x1) ** 2 + (y - y1) ** 2 < r1 ** 2
    mask_circle2 = (x - x2) ** 2 + (y - y2) ** 2 < r2 ** 2
    image = np.logical_or(mask_circle1, mask_circle2)

    # Now we want to separate the two objects in image
    # Generate the markers as local maxima of the distance to the background
    distance = ndi.distance_transform_edt(image)
    local_maxi = peak_local_max(distance, indices=False, footprint=np.ones((3, 3)),
                                labels=image)
    markers = ndi.label(local_maxi)[0]
    labels = watershed(-distance, markers, mask=image)

    fig, axes = plt.subplots(ncols=3, figsize=(9, 3), sharex=True, sharey=True,
                             subplot_kw={'adjustable': 'box-forced'})
    ax = axes.ravel()

    ax[0].imshow(image, cmap=plt.cm.gray, interpolation='nearest')
    ax[0].set_title('Overlapping objects')
    ax[1].imshow(-distance, cmap=plt.cm.gray, interpolation='nearest')
    ax[1].set_title('Distances')
    ax[2].imshow(labels, cmap=plt.cm.spectral, interpolation='nearest')
    ax[2].set_title('Separated objects')

    for a in ax:
        a.set_axis_off()

    fig.tight_layout()
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("dir_path", help="")
    args = parser.parse_args()
    main(args)
    # auto_watershed(args)

