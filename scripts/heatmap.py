"""generates heatmaps of variables, takes directory as argument and uses
the segmented image and output.csv to calculate heatmaps. output to ./heatmaps"""

import os
import json
import argparse

import matplotlib
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import tqdm

import common_functions as cf


def load_json_data(json_path):
    with open(json_path, 'r') as json_handle:
        return json.load(json_handle)


def main():
    exp_dir = args.exp_dir
    seg_path = cf.get_seg_path(exp_dir)
    if not os.path.exists(seg_path):
        print "Segmented Image not found\n"
        return

    print seg_path
    id_array = cf.path2id_array(seg_path)

    data_path = os.path.join(exp_dir, "data.json")
    if not os.path.exists(data_path):
        print "Data json file not found"
        return
    data_dict = load_json_data(data_path)

    outline_array = cf.generate_cell_outline_array(id_array)

    if not os.path.exists(os.path.join(exp_dir, "heatmaps")):
        os.mkdir(os.path.join(exp_dir, "heatmaps"))

    random_data = data_dict.itervalues().next()
    data_to_plot = random_data.keys()

    size, vox = cf.read_dims(exp_dir)

    color_scheme = 'viridis'
    heatmap_shape = [id_array.shape[0], id_array.shape[1], 4]
    # TODO check extent

    extent = [0, heatmap_shape[1] * vox[1], 0, heatmap_shape[0] * vox[0]]

    xs, ys = [], []
    for cell_data in data_dict.itervalues():
        xs.append(cell_data['Centroid-x_um'])
        ys.append(cell_data['Centroid-y_um'])

    im_height = id_array.shape[1] * vox[1]
    border = 0.07 * im_height

    xlims = (min(xs) - border, max(xs) + border)
    ylims = ((im_height - max(ys)) - border, (im_height - min(ys)) + border)

    def units_string(data_type):
        units = data_type.split("_")[1]
        if units == "um":
            return r'$\mu$m'
        if units == "um2":
            return r'$\mu$m$^{2}$'
        if units == "pixels":
            return "Pixels"
        else:
            return ""

    def do_heatmap(data_type):
        data_list = []

        for cell_data in data_dict.itervalues():
            data_list.append(cell_data[data_type])
        print "Painting", data_type, ": ", min(data_list), max(data_list)

        if "area" in data_type:
            norm = colors.LogNorm(vmin=min(data_list), vmax=max(data_list))
            color_map = matplotlib.cm.ScalarMappable(cmap=color_scheme, norm=norm)
        else:
            color_map = matplotlib.cm.ScalarMappable(cmap=color_scheme)

        color_map.set_clim(vmin=min(data_list), vmax=max(data_list))

        heatmap = np.full(shape=heatmap_shape, fill_value=[0, 0, 0, 255], dtype='float32')

        for cell_id, cell_data in tqdm.tqdm(data_dict.iteritems()):
            if 'Density' in data_type and float(cell_data[data_type]) < 0.7:
                color = color_map.to_rgba(cell_data[data_type])
                col_list = list(color)
                col_list[3] = 0.5
                heatmap[id_array == int(cell_id)] = col_list
            else:
                color = color_map.to_rgba(cell_data[data_type])
                heatmap[id_array == int(cell_id)] = color

        plt.figure(figsize=(10, 10))
        ax1 = plt.subplot(111)
        plt.title(data_type.split("_")[0])
        plt.imshow(heatmap, cmap=color_scheme, interpolation='nearest', extent=extent)
        plt.imshow(outline_array, cmap=color_scheme, interpolation='nearest', extent=extent)
        plt.clim(min(data_list), max(data_list))
        cf.add_scale_bar(ax1)
        ax1.axis('on')
        ax1.get_xaxis().set_ticks([])
        ax1.get_yaxis().set_ticks([])
        plt.xlim(xlims)
        plt.ylim(ylims)
        divider1 = make_axes_locatable(ax1)
        cax1 = divider1.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(cax=cax1)
        cbar.set_label(units_string(data_type))
        # plt.show()
        plt.savefig(os.path.join(exp_dir, "heatmaps", data_type + ".png"), format='png', dpi=1000)

    def do_cell_outlines():

        print "Painting Cell Outlines"

        plt.figure(figsize=(10, 10))
        ax1 = plt.subplot(111)
        plt.title("Outline")
        plt.imshow(outline_array, cmap=color_scheme, interpolation='nearest', extent=extent)
        cf.add_scale_bar(ax1)
        ax1.axis('on')
        ax1.get_xaxis().set_ticks([])
        ax1.get_yaxis().set_ticks([])
        plt.xlim(xlims)
        plt.ylim(ylims)

        plt.savefig(os.path.join(exp_dir, "heatmaps", "oultine" + ".png"), format='png', dpi=1000)

    for data_type in data_to_plot:
        do_heatmap(data_type)

    do_cell_outlines()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("exp_dir", help="Experiment Directory")
    args = parser.parse_args()
    main()
