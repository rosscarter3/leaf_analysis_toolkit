import os
import sys
import argparse

import numpy as np
import skimage.morphology
import matplotlib.pyplot as plt
from PIL import Image, ImageFile

# from openalea.image.algo.basic import logicalnot
# from openalea.image.all import imread, imsave
# from openalea.image.spatial_image import SpatialImage
# from vplants.asclepios.vt_exec.connexe import hysteresis, connected_components
# from vplants.asclepios.vt_exec.regionalmax import regionalmax
# from vplants.asclepios.vt_exec.watershed import watershed
# from vplants.mars_alt.mars import segmentation
from scipy.ndimage.filters import gaussian_filter
from tissueviewer.tvtiff import tiffread, tiffsave

import common_functions as cf

from scipy.ndimage import imread



from scipy import ndimage as nd
from scipy.ndimage.morphology import grey_erosion, grey_dilation


def filtering(img, filter_type="gaussian", filter_value=0.5):
    """
    :Parameters:
    - `image` (openalea.image.SpatialImage) - image

    - `filter_type`(str) - denoising method used for filtering ("gaussian" or "asf" for alternate sequential filter).
                           default is "gaussian".

    - `filter_value` (float for "gaussian" filter or int for alternate sequential filter) - value used for the filtering :
                                            * for a Gaussian filtering, the "filter_value" corresponds to the standard deviation.
                                            * for a Alternate Sequential Filter, the "filter_value" corresponds to the number of
                                              succession of morphological opening and closing operations.
    """
    # if not isinstance(img, SpatialImage):
    #     img = SpatialImage(img)
    #
    # if filter_type == 'gaussian':
    #     if not isinstance(filter_value, float):
    #         raise RuntimeError, 'value used for Gaussian filtering must be a float type'
    #     else:
    #         img = recfilters(img, filter_type="sigma", filter_value=filter_value, Trueyz=(0, 0, 0))

    if filter_type == 'asf':
        if not isinstance(filter_value, int):
            raise RuntimeError('value used for the Alternate Sequential Filter must be a integer type')
        else:
            for rad in range(1, ((filter_value + 1) / 2) + 1):
                print("closing operations with structuring elements of size %s" % rad)
                struct = euclidean_sphere(rad)
                # ~ s=(rad,rad,rad)
                img = grey_dilation(img, footprint=struct)
                img = grey_erosion(img, footprint=struct)

                if filter_value >= rad * 2:
                    print("opening operations with structuring elements of size %s" % rad)
                    img = grey_erosion(img, footprint=struct)
                    img = grey_dilation(img, footprint=struct)
    else:
        raise RuntimeError('filter type not supported')
    return img


def euclidean_sphere(size):
    """
    Generate a euclidean sphere for binary morphological operations

    :Parameters:
        - `size` (int) - the shape of the euclidean sphere = 2*size + 1.

    :Returns:
        - Euclidean sphere which may be used for binary morphological operations, with shape equal to 2*size + 1.
    """
    n = 2*size + 1
    sphere = np.zeros((n,n,n),np.bool)
    for x in range(n):
        for y in range(n):
            for z in range(n):
                if abs(x-size)+abs(y-size)+abs(z-size)<=2*size:
                    sphere[x,y,z]=True
    return sphere


def logicalnot (img):
    """
    """
    d = img.dtype
    vmax = np.iinfo(d).max
    im_target = vmax * np.ones(img.shape, img.dtype)

    image = np.bitwise_xor(img,im_target)
    return image


def max_proj(array3d, direction=0):
    return np.amax(array3d, axis=direction)


def stdev_proj(array3D, direction=0):
    return np.std(array3D, axis=direction)


def segment_tv(imageFName, hmin, asf, sigma, resolution=(1, 1, 1), checkBackground=True, background=1):
    from skimage import exposure
    import libtiff
    import skimage.morphology as morphology
    from skimage.filters import rank

    image = libtiff.TiffFile(imageFName)
    image = image.get_tiff_array()
    image = np.array(image, dtype='uint8')

    plt.imshow(max_proj(image))
    plt.show()

    im_filtered = filtering(image, filter_type='asf', filter_value=int(2))
    im_tmp = logicalnot(im_filtered)

    #im_tmp = regionalmax(im_tmp, hmin)
    # im_tmp = hysteresis(im_tmp, 1, hmin, connectivity=6)
    # seeds = connected_components(im_tmp, 1)
    # seg = watershed(seeds, im_filtered)
    plt.subplot(121)
    plt.imshow(stdev_proj(im_tmp))
    plt.subplot(122)
    plt.imshow(max_proj(im_filtered))
    plt.show()

    return

    maxi_proj = max_proj(image)
    slice = 100

    im_filtered = gaussian_filter(image, sigma)
    max_proj_gauss = max_proj(im_filtered)

    img_rescale = exposure.equalize_hist(im_filtered)

    # selem = morphology.cube(25)
    # img_eq = rank.equalize(im_filtered, selem)

    # plt.subplot(121)
    # plt.imshow(max_proj(img_rescale))
    # plt.subplot(122)
    # plt.imshow(max_proj(img_eq))
    # plt.show()

    plt.subplot(121)
    plt.imshow(maxi_proj)
    plt.subplot(122)
    plt.imshow(max_proj_gauss)
    plt.show()

    # plt.subplot(121)
    # plt.imshow(image[slice])
    # plt.imshow(122)
    # plt.imshow(im_filtered[slice])
    # plt.show()

    # im_filtered = segmentation.filtering(im_filtered, filter_type="asf", filter_value=asf)
    # im_tmp = logicalnot(im_filtered)
    # im_tmp = regionalmax(im_tmp, hmin)
    # im_tmp = hysteresis(im_tmp, 1, hmin, connectivity=6)
    # seeds = connected_components(im_tmp, 1)
    # seg = watershed(seeds, im_filtered)
    # seg.resolution = resolution

    # if checkBackground:
    #     maxLabel = np.max(seg)
    #     print "Checking the background label ..."
    #     largestCellLabel = np.bincount(seg.flatten()).argmax()
    #     if largestCellLabel != background:
    #         np.putmask(seg, seg == background, maxLabel + 1)
    #         np.putmask(seg, seg == largestCellLabel, background)
    # segImageFName = imageFName[:-4] + "_hmin_%d_asf_%1.2f_s_%1.2f.tif" % (hmin, asf, sigma)
    # tiffsave(makeImageLabelsConsecutive(seg)[0], segImageFName, **tags)
    # return segImageFName


def segment(tiff_file):
    print("Segmenting, ", tiff_file)

    return 0


def open_alea_ex(args):
    # from openalea.image import imread, display
    import openalea.image
    im3 = imread(args.tiff_file)
    w3 = display(im3)


def main(args):
    print("Segmentation running!")

    print("im path is: ", args.tiff_file)

    segment_tv(args.tiff_file, 0.3, 0.3, 2)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("tiff_file", help="path to tiff file for segmentation")
    args = parser.parse_args()
    main(args)
    #open_alea_ex(args)
