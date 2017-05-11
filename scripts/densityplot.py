"""calculates and displays/saves cell density plots for rachel's leaf data"""

import sys
import os
import datetime
import csv
import numpy as np
from scipy import stats
from PIL import Image
import matplotlib.pyplot as plt


def read_dims(idir):
    """reads the dimensions file located in idir and returns
       image and voxel dimensions"""

    with open(os.path.join(idir, "dims.txt"), 'r') as file_handle:
        dimstr = file_handle.read()
        dimstr = [float(s) for s in dimstr.split(",")]
        size_x, size_y, size_z = dimstr[0], dimstr[1], dimstr[2]
        voxel_x, voxel_y, voxel_z = dimstr[3], dimstr[4], dimstr[5]
    return size_x, size_y, size_z, voxel_x, voxel_y, voxel_z


def flush_message(message):
    """flushes a message"""
    print message,
    sys.stdout.flush()


def return_centroids(csvpath, voxel_x, voxel_y):
    """returns cell centroids from csv file at csvpath"""
    with open(csvpath, 'rU') as file_handle:
        reader = csv.DictReader(file_handle)
        x_centroid = []
        y_centroid = []
        for row in reader:
            x_centroid.append(float(row['xcentroid_pix']) * voxel_x)
            y_centroid.append(float(row['ycentroid_pix']) * voxel_y)
            values = np.vstack([x_centroid, y_centroid])
    return x_centroid, y_centroid, values


def main(args):
    """calculates and outputs cell density plots"""
    idir = args[1]

    # print idir

    start_time = datetime.datetime.now()

    print "Start time: %s" % start_time

    # impath = os.path.join(idir, "00000.png")
    # im_file = idir[5:]+'_max-proj.png'
    im_file = os.path.split(idir)[-2][5:]+'_max-proj.png'

    im_path = os.path.join(idir, im_file)

    # print im_path

    csvpath = os.path.join(idir, "output.csv")

    size_x, size_y, _, voxel_x, voxel_y, _ = read_dims(idir)

    grid_x_points, grid_y_points = np.mgrid[0:size_x*voxel_x:100j, 0:size_y*voxel_y:100j]
    positions = np.vstack([grid_x_points.ravel(), grid_y_points.ravel()])

    x_centroid, y_centroid, values = return_centroids(csvpath, voxel_x, voxel_y)

    kernel = stats.gaussian_kde(values)
    density = np.reshape(kernel(positions).T, grid_x_points.shape)
    flipped_density = np.rot90(density,3)
    relative_density = flipped_density/np.max(flipped_density)

    leaf_image = Image.open(im_path)

    print "Time: ", (datetime.datetime.now() - start_time)
    print ""

    # in microns
    bin_size = 15 * 2

    min_x = min(x_centroid)
    min_y = min(y_centroid)
    # max_x = max(x_centroid)
    # max_y = max(y_centroid)

    x_centroid[:] = [v - min_x for v in x_centroid]
    y_centroid[:] = [v - min_y for v in y_centroid]

    border = 10

    x_range_min = min(x_centroid) - border
    x_range_max = max(x_centroid) + border

    y_range_min = min(y_centroid) - border
    y_range_max = max(y_centroid) + border

    crop_image = np.fliplr(np.flipud(np.rot90(leaf_image)))
    image_extent = [size_x*voxel_x - min_x, -min_x,
                    size_y*voxel_y - min_y, -min_y]

    color_map = plt.cm.gist_earth_r

    plt.Figure(figsize=(200, 32), dpi=10000)

    plt.imshow(crop_image, alpha=1, cmap=plt.cm.Greys,
               extent=image_extent)

    plt.imshow(relative_density, cmap=color_map,
               extent=image_extent, alpha=0.5)

    plt.contour(np.flipud(relative_density), [0.7], extent=image_extent)

    # plt.plot(x_centroid, y_centroid, 'k.', markersize=2)
    plt.xlim([x_range_min, x_range_max])
    plt.ylim([y_range_min, y_range_max])
    plt.title(im_file[:-13])

    plt.subplots_adjust(left=0.02, bottom=0.02, right=0.98, top=0.98,
                        wspace=0.1, hspace=0.05)

    _, data_dir = os.path.split(idir)

    out_dir = os.path.join(data_dir, 'cell_densities')

    out_name = im_file + '_density_nocentroids.pdf'
    out_file = os.path.join(idir, out_name)

    plt.savefig(out_file, bbox_inches='tight')


if __name__ == '__main__':
    sys.exit(main(sys.argv))
