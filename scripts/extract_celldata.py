#!/usr/bin/env python
"""extracts quantitative values and performs a kernel density
estimation from segmented images and writes to .csv file"""

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
    values = None
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


def return_coefficient(t_x, t_y, b_x, b_y):
    # equation of line though tip and base
    m_mv = (t_y - b_y) / (t_x - b_x)
    c_mv = t_y - m_mv * t_x

    # equation of line though tip, perpendicular to tip-base line
    m_tip = -1 / m_mv
    c_tip = t_y - m_tip * t_x

    a = m_tip
    b = -1
    c = c_tip
    return a, b, c


def return_distances(a, b, c, t_x, t_y, b_x, b_y, x_c, y_c):
    distance_mv_pix = np.abs((b_y - t_y) * x_c - (b_x - t_x) * y_c + b_x * t_y - b_y * t_x) / np.sqrt(
        (b_y - t_y) ** 2 + (b_x - t_x) ** 2)
    distance_tip_pix = np.abs(a * x_c + b * y_c + c) / np.sqrt(a ** 2 + b ** 2)
    return distance_mv_pix, distance_tip_pix


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
        print "\"tip.txt\" not found\n"
        return

    print "extracting data from: ", os.path.basename(seg_path)

    tip, base = read_tip(tip_path)

    (t_x, t_y) = (float(tip[1]), float(size[0] - tip[0]))
    (b_x, b_y) = (float(base[1]), float(size[0] - base[0]))

    (a, b, c) = return_coefficient(t_x, t_y, b_x, b_y)

    (t_x_r, t_y_r) = (float(t_x) * voxel[1], float(t_y) * voxel[0])
    (b_x_r, b_y_r) = (float(b_x) * voxel[1], float(b_y) * voxel[0])

    a_r, b_r, c_r = return_coefficient(t_x_r, t_y_r, b_x_r, b_y_r)

    # calculate cell level data

    cell_data_dict = {}
    cell_info = {}

    cell_props = skim.regionprops(id_array)

    for cell in cell_props:
        cell_id = cell['label']
        area_real = cell['area'] * voxel[0] * voxel[1]
        centroid_x_pixels = cell['centroid'][1]
        centroid_y_pixels = cell['centroid'][0]

        (x_c, y_c) = (centroid_x_pixels, size[0] - centroid_y_pixels)
        distance_mv_pix, distance_tip_pix = return_distances(a, b, c, t_x, t_y, b_x, b_y, x_c, y_c)

        (x_c_r, y_c_r) = (centroid_x_pixels * voxel[1], (size[0] - centroid_y_pixels) * voxel[0])
        distance_mv_real, distance_tip_real = return_distances(a_r, b_r, c_r, t_x_r, t_y_r, b_x_r, b_y_r, x_c_r, y_c_r)

        centroid_x_real = centroid_x_pixels * voxel[1]
        centroid_y_real = centroid_y_pixels * voxel[0]

        perimeter = cell['perimeter'] * voxel[0]

        circularity = (4 * np.pi * area_real) / perimeter ** 2

        cell_info = {'Area_um2': area_real,
                     'Centroid-x_pixels': int(centroid_x_pixels),
                     'Centroid-y_pixels': int(centroid_y_pixels),
                     'Centroid-x_um': float(centroid_x_real),
                     'Centroid-y_um': float(centroid_y_real),
                     'Perimeter_um': float(perimeter),
                     'Distance-from-mv_pixels': np.abs(distance_mv_pix),
                     'Distance-from-tip_pixels': np.abs(distance_tip_pix),
                     'Distance-from-mv_um': np.abs(distance_mv_real),
                     'Distance-from-tip_um': np.abs(distance_tip_real),
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
    leaf_data["genotype"] = dir_list
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
