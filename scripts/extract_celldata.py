#!/usr/bin/env python
"""extracts quantitative values from segmented images and writes to .csv file"""

import os
import csv
import json
import argparse

import numpy as np
from scipy import stats
import skimage.measure as skim

import common_functions as cf


def read_tip(tip_path):
    """reads the tip file and returns tip/petiole coordinates"""
    with open(tip_path, 'r') as file_handle:
        dimstr = file_handle.read()
    dimstr = [float(s) for s in dimstr.split(",")]
    tip_x, tip_y, base_x, base_y = dimstr[0], dimstr[1], dimstr[2], dimstr[3]
    return (tip_x, tip_y), (base_x, base_y)


def do_kde(size, voxel, data_dict):
    size_x = size[0]
    size_y = size[1]
    voxel_x = voxel[0]
    voxel_y = voxel[1]

    x_centroid = []
    y_centroid = []
    for cid, data, in data_dict.iteritems():
        x_centroid.append(data['Centroid-x_um'])
        y_centroid.append(data['Centroid-y_um'])
        values = np.vstack([x_centroid, y_centroid])

    grid_x_points, grid_y_points = np.mgrid[0:size_x * voxel_x:np.complex(0, size_x),
                                            0:size_y * voxel_y:np.complex(0, size_y)]
    # grid_x_points, grid_y_points = np.mgrid[0:size_x * voxel_x:500j,
    #                                         0:size_y * voxel_y:500j]

    positions = np.vstack([grid_x_points.ravel(), grid_y_points.ravel()])
    print "Calculating KDE\n"
    kernel = stats.gaussian_kde(values)
    density = np.reshape(kernel(positions).T, grid_x_points.shape)
    flipped_density = np.rot90(density, 3)
    relative_density = flipped_density / np.max(flipped_density)
    return np.fliplr(relative_density)


def main():
    """main function extracts and writes data"""
    exp_dir = args.exp_dir

    seg_path = cf.get_seg_path(exp_dir)
    if not os.path.exists(seg_path):
        print "No segmented image found\n"
        return
    id_array = cf.path2id_array(seg_path)

    size, voxel = cf.read_dims(exp_dir)

    tip_path = os.path.join(exp_dir, "tip.txt")
    if not os.path.exists(tip_path):
        print "\"tips.txt\" not found\n"
        return

    print "extracting data from: ", os.path.basename(seg_path)

    tip, base = read_tip(tip_path)

    theta = np.arctan((base[1]-tip[1])/(base[0]-tip[0]))

    vxrot = voxel[1] * np.cos(theta) - voxel[0] * np.sin(theta)
    vyrot = voxel[1] * np.sin(theta) + voxel[0] * np.cos(theta)

    # calcualte cell level data

    cell_data_dict = {}

    cell_props = skim.regionprops(id_array)

    for cell in cell_props:
        cell_id = cell['label']
        area_real = cell['area'] * voxel[0] * voxel[1]
        centroid_x_pixels = cell['centroid'][1]
        centroid_y_pixels = cell['centroid'][0]

        centroid_x_real = centroid_x_pixels * voxel[0]
        centroid_y_real = centroid_y_pixels * voxel[1]

        perimeter = cell['perimeter'] * voxel[0]

        xtr_pix = centroid_x_pixels - tip[1]
        ytr_pix = centroid_y_pixels - tip[0]
        xrot_pix = xtr_pix * np.cos(theta) - ytr_pix * np.sin(theta)
        yrot_pix = xtr_pix * np.sin(theta) + ytr_pix * np.cos(theta)
        xrot_re = xrot_pix * vxrot
        yrot_re = yrot_pix * vyrot

        circularity = (4 * np.pi * area_real) / perimeter ** 2

        cell_info = {'Area_um2': area_real,
                     'Centroid-x_pixels': int(centroid_x_pixels),
                     'Centroid-y_pixels': int(centroid_y_pixels),
                     'Centroid-x_um': float(centroid_x_real),
                     'Centroid-y_um': float(centroid_y_real),
                     'Perimeter_um': float(perimeter),
                     'Distance-from-mv_pixels': np.abs(xrot_pix),
                     'Distance-from-tip_pixels': np.abs(yrot_pix),
                     'Distance-from-mv_um': np.abs(float(xrot_re)),
                     'Distance-from-tip_um': np.abs(float(yrot_re)),
                     'Circularity_none': circularity}

        # TODO ADD MORE STUFF HERE IF NEEDED
        # Name of the dictionary key of the form: "Name_Units"

        cell_data_dict[cell_id] = cell_info

    density_array = do_kde(size, voxel, cell_data_dict)
    for cell_id in cell_data_dict.iterkeys():
        av_density = np.mean(density_array[id_array == int(cell_id)])
        cell_data_dict[cell_id]['Relative-Cell-Density_none'] = av_density

    # automatic csv writing for each key in cell data dictionary

    csv_headings = [s for s in cell_info.iterkeys()]
    csv_headings.insert(0, 'Cell_ID')

    csv_path = os.path.join(exp_dir, "data.csv")
    with open(csv_path, 'wb') as csv_file:
        csvwriter = csv.writer(csv_file, delimiter=',')
        csvwriter.writerow(csv_headings)
        for cell_id, data_dict in cell_data_dict.iteritems():
            data_list = [cell_id]
            for data_type in csv_headings[1:]:
                data_list.append(data_dict[data_type])
            csvwriter.writerow(data_list)

    # export the dictionary as .json file

    json_path = os.path.join(exp_dir, "data.json")
    with open(json_path, 'w') as json_file:
        json.dump(cell_data_dict, json_file, indent=2)

    # calculate leaf level data
    # TODO needs work for handling of directory name

    leaf_data = {}

    dir_list = args.exp_dir.split('/')[-2].split('_')

    leaf_data["timepoint"] = dir_list[1]
    leaf_data["genotype"] = dir_list[3].split('.')[0]
    leaf_data["no-of-cells"] = len(cell_data_dict)
    area = 0
    for data_dict in cell_data_dict.itervalues():
        area += data_dict['Area_um2']
    leaf_data["leaf-area_um2"] = area



    leaf_data_csv_data_types = leaf_data.keys()
    csv_path = os.path.join(exp_dir, "leaf_data.csv")
    with open(csv_path, 'wb') as csv_file:
        csvwriter = csv.writer(csv_file, delimiter=',')
        csvwriter.writerow(leaf_data_csv_data_types)
        data_list = []
        for data_type in leaf_data_csv_data_types:
            data_list.append(leaf_data[data_type])
        csvwriter.writerow(data_list)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("exp_dir", help="Experiment Directory")
    args = parser.parse_args()

    main()
