"""generates heatmaps of variables, takes directory as argument and uses
the segmented image and output.csv to calculate heatmaps. output to ./heatmaps"""
import sys
import os
import datetime
import csv
import json

import matplotlib.pyplot as plt
import matplotlib.cm
from matplotlib.colors import LogNorm
import numpy as np

from jicbioimage.illustrate import Canvas
import tensortools.functions as ttf
from celldata import CellData


def read_dims(idir):
    """reads dimensions file in idir and returns image dimensions
    and voxel dimensions"""
    with open(os.path.join(idir, "dims.txt"), 'r') as file_handle:
        dimstr = file_handle.read()
    dimstr = [float(s) for s in dimstr.split(",")]
    size_x, size_y, size_z = dimstr[0], dimstr[1], dimstr[2]
    voxel_x, voxel_y, voxel_z = dimstr[3], dimstr[4], dimstr[5]
    return size_x, size_y, size_z, voxel_x, voxel_y, voxel_z


def norm_var(idir, var):
    """returns dict with variable scaled as (v / max)"""
    flush_message("Reading %s data... " % var)
    vari = {}

    with open(os.path.join(idir, "output.csv"), 'r') as file_handle:
        reader = csv.DictReader(file_handle)
        for row in reader:
            vari[row["cid"]] = row[var]

    max_vari = float(max(vari.itervalues()))
    norm_vari = {k: (float(v) / max_vari) for k, v in vari.iteritems()}
    print "Done"

    return norm_vari


def scaled_var(idir, var):
    """returns dict with variable scaled as ((v - min) / max)"""
    flush_message("Reading %s data... " % var)
    vari = {}

    with open(os.path.join(idir, "output.csv"), 'r') as file_handle:
        reader = csv.DictReader(file_handle)
        for row in reader:
            vari[row["cid"]] = row[var]

    min_vari = float(min(vari.itervalues()))
    shifted_varis = {k: (float(v) - min_vari) for k, v in vari.iteritems()}
    max_shifted = float(max(vari.itervalues()))
    scaled_vari = {k: (float(v) / max_shifted) for k, v in shifted_varis.iteritems()}
    print "Done"

    return scaled_vari


def flush_message(message):
    """flushes message to command line"""
    print message
    sys.stdout.flush()


def paint_cells(var, cdata, idir, posterise=False, scaled=False):
    """creates a canvas and paints the scaled variable data to each cell"""
    if scaled:
        vari_dict = scaled_var(idir, var)
    else:
        vari_dict = norm_var(idir, var)

    size_x, size_y, _, _, _, _ = read_dims(idir)
    canvas = Canvas((size_x, size_y, 4))
    canvas2 = Canvas((size_x, size_y))

    var_list = list(vari_dict.itervalues())
    var_list.remove(max(var_list))

    color_scheme = 'Set1'
    color_map = matplotlib.cm.ScalarMappable(cmap=color_scheme,
                                             norm=LogNorm(vmin=min(var_list),
                                                          vmax=max(var_list)))

    # color_map.set_clim(vmin=min(var_list), vmax=max(var_list))

    flush_message("Painting cells... ")
    for (cid, cell) in cdata:
        if cell.area < 4000:
            cell_mask = np.zeros((size_x, size_y), dtype=bool)
            for point in cell.pl:
                cell_mask[point] = True
            cell_id = vari_dict[str(cid)]
            # canvas.mask_region(cell_mask, (0, cell_id * 255, 0))
            color = color_map.to_rgba(cell_id)
            canvas.mask_region(cell_mask, color)
            canvas2.mask_region(cell_mask, cell_id)
    print "Done"
    with open(os.path.join(idir, "heatmaps/%s.png" % var), "wb") as file_handle:
        file_handle.write(canvas.png())


def load_json(idir):
    json_path = os.path.join(idir, 'output_json.json')
    with open(json_path, 'r') as json_handle:
        data_dict = json.load(json_handle)
    return data_dict


def paint_cells2(var, cid_array, idir):
    """creates a canvas and paints the scaled variable data to each cell"""

    data_dict = load_json(idir)

    shape = [cid_array.shape[0], cid_array.shape[1], 4]
    canvas = np.zeros(shape)
    canvas2 = np.zeros(shape)
    var_list = []

    for data in data_dict.itervalues():
        var_list.append(data[var])

    color_scheme = 'Set1'
    color_map = matplotlib.cm.ScalarMappable(cmap=color_scheme)
    color_map.set_clim(vmin=min(var_list), vmax=max(var_list))

    flush_message("Painting cells... ")

    cells = (np.unique(cid_array))

    for cid, data in data_dict.iteritems():
        canvas[cid_array == int(cid)] = color_map.to_rgba(data[var])

    print "Done"

    cell_outlines = ttf.generate_cell_outline_array(cid_array)

    plt.imshow(canvas, interpolation='nearest')
    plt.imshow(cell_outlines, interpolation='nearest')
    plt.show()

    #with open(os.path.join(idir, "heatmaps/%s_from_json.png" % var), "wb") as file_handle:
    #    file_handle.write(canvas)


def main(args):
    """"reads the data and paints the pictures"""
    start_time = datetime.datetime.now()
    idir = args[1]
    impath = os.path.join(idir, "00000.png")

    mapsdir = os.path.join(idir, "heatmaps")

    cid_array = ttf.path2id_array(impath)

    if not os.path.exists(mapsdir):
        os.makedirs(mapsdir)

    #cdata = CellData(impath)

    # paint_cells2("yrot_re", cid_array, idir)
    # paint_cells("area_pix", cdata, idir)
    # paint_cells("perimeter_pix", cdata, idir)
    # paint_cells("solidity", cdata, idir, scaled=True)
    # paint_cells("circularity", cdata, idir)
    # paint_cells("eccentricity", cdata, idir)
    # paint_cells("axis_ratio", cdata, idir)
    # paint_cells("dist_tip_re", cdata, idir)
    # paint_cells("dist_mv_re", cdata, idir)
    # paint_cells("stomata", cdata, idir)
    # paint_cells_invarea(cdata,idir)

    print "Time: ", (datetime.datetime.now() - start_time)
    print ""

if __name__ == '__main__':
    sys.exit(main(sys.argv))
