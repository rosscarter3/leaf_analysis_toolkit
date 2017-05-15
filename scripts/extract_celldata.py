#!/usr/bin/env python
"""extracts quantitative values from segmented images and writes to .csv file"""

import os
import csv
import json
import argparse

import numpy as np
import skimage.measure as skim

import common_functions as cf


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
    im_path = args.im_path
    im_dir = os.path.dirname(im_path)
    id_array = cf.path2id_array(im_path)

    _, _, _, voxel_x, voxel_y, _ = read_dims(im_dir)
    # tip_x, tip_y, base_x, base_y = read_tip(idir)

    # theta = np.arctan((base_x-tip_x)/(base_y-tip_y))

    # vxrot = voxel_x * np.cos(theta) - voxel_y * np.sin(theta)
    # vyrot = voxel_x * np.sin(theta) + voxel_y * np.cos(theta)

    cell_ids = np.unique(id_array)
    if 0 in cell_ids:
        del cell_ids[cell_ids.index(0)]

    cell_data_dict = {}

    cell_props = skim.regionprops(id_array)

    for cell in cell_props:
        cell_id = cell['label']
        cell_info = {'Area_Real': cell['area'] * voxel_x * voxel_y,
                     'Centroid_x_Pixels': cell['centroid'][1],
                     'Centroid_y_Pixels': cell['centroid'][0],
                     'Perimeter_Real': cell['perimeter'] * voxel_x}

        #ADD MORE STUFF HERE IF NEEDED

        cell_data_dict[cell_id] = cell_info

    csv_headings = [s for s in cell_info.iterkeys()]
    csv_headings.insert(0, 'Cell_ID')

    csv_path = os.path.join(im_dir, "data.csv")
    with open(csv_path, 'wb') as csv_file:
        csvwriter = csv.writer(csv_file, delimiter=',')
        csvwriter.writerow(csv_headings)
        for cell_id, data_dict in cell_data_dict.iteritems():
            data_list = [cell_id]
            for data_type in csv_headings[1:]:
                data_list.append(data_dict[data_type])
            csvwriter.writerow(data_list)

    json_path = os.path.join(im_dir, "data.json")
    with open(json_path, 'w') as json_file:
        json.dump(cell_data_dict, json_file, indent=2)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("im_path", help="Segmented Image")
    args = parser.parse_args()

    main()
