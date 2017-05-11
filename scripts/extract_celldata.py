#!/usr/bin/env python
"""extracts quantitative values from segmented images and writes to .csv file"""
import sys
import os
import csv
import datetime
import argparse

import numpy as np
from PIL import Image
import skimage.measure as skim

import tensortools.functions as ttf


def flush_message(message):
    """flushes message to command line"""
    print message,
    sys.stdout.flush()


def read_dims(idir):
    """reads the dimensions file located in idir and returns
       image and voxel dimensions"""

    with open(os.path.join(idir, "dims.txt"), 'r') as file_handle:
        dimstr = file_handle.read()
        dimstr = [float(s) for s in dimstr.split(",")]
        size_x, size_y, size_z = dimstr[0], dimstr[1], dimstr[2]
        voxel_x, voxel_y, voxel_z = dimstr[3], dimstr[4], dimstr[5]
    return size_x, size_y, size_z, voxel_x, voxel_y, voxel_z


def read_tip(idir):
    """reads the tip file and returns tip/petiole coordinates"""
    with open(os.path.join(idir, "tip.txt"), 'r') as file_handle:
        dimstr = file_handle.read()
    dimstr = [float(s) for s in dimstr.split(",")]
    tip_x, tip_y, base_x, base_y = dimstr[0], dimstr[1], dimstr[2], dimstr[3]
    return tip_x, tip_y, base_x, base_y


def main():
    """main function extracts and writes data"""
    start_time = datetime.datetime.now()
    idir = args.input_directory
    # impath = os.path.join(idir, "00000_prefiltered.png")
    impath = os.path.join(idir, "00000.png")

    seg_image = Image.open(impath)
    seg_image_array_rgb = np.asarray(seg_image, dtype="uint64")

    cid_array = ttf.path2id_array(impath)

    print np.unique(cid_array)

    _, _, _, voxel_x, voxel_y, _ = read_dims(idir)
    tip_x, tip_y, base_x, base_y = read_tip(idir)

    theta = np.arctan((base_x-tip_x)/(base_y-tip_y))

    vxrot = voxel_x * np.cos(theta) - voxel_y * np.sin(theta)
    vyrot = voxel_x * np.sin(theta) + voxel_y * np.cos(theta)

    flush_message("Writing cell data... ")
    with open(os.path.join(idir, "output_test.csv"), "w") as file_handle:
        writer = csv.writer(file_handle)

        writer.writerow(["cid",
                         "x_centroid_pix",
                         "y_centroid_pix",
                         "area_pix",
                         "x_centroid_re",
                         "y_centroid_re",
                         "area_re",
                         "dist_mv_re",
                         "dist_tip_re",
                         "perimeter_pix",
                         "eccentricity",
                         "solidity"])

        for cell in skim.regionprops(cid_array):
            cid = cell.label

            xcentroid_pix = cell.centroid[1]
            ycentroid_pix = cell.centroid[0]

            xcentroid_re = xcentroid_pix * voxel_x
            ycentroid_re = ycentroid_pix * voxel_y
            area_re = cell.area * voxel_x * voxel_y
            xtr_pix = xcentroid_pix - tip_x
            ytr_pix = ycentroid_pix - tip_y
            xrot_pix = xtr_pix * np.cos(theta) - ytr_pix * np.sin(theta)
            yrot_pix = xtr_pix * np.sin(theta) + ytr_pix * np.cos(theta)
            xrot_re = xrot_pix * vxrot
            yrot_re = yrot_pix * vyrot

            writer.writerow([cid,
                             xcentroid_pix,
                             ycentroid_pix,
                             cell.area,
                             xcentroid_re,
                             ycentroid_re,
                             area_re,
                             xrot_re,
                             yrot_re,
                             cell.perimeter,
                             cell.eccentricity,
                             cell.solidity])

    print "Done"
    print "Time: ", (datetime.datetime.now() - start_time)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_directory", help="Input Directory")
    args = parser.parse_args()

    main()
