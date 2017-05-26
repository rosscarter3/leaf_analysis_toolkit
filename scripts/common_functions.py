"""module containing common functions"""
import platform
import json
import subprocess
import re
import os.path

import numpy as np
from PIL import Image


def read_dims(idir):
    """reads the dimensions file located in idir and returns
       image and voxel dimensions"""

    with open(os.path.join(idir, "dims.txt"), 'r') as file_handle:
        dimstr = file_handle.read()
    dimstr = [float(s) for s in dimstr.split(",")]
    size_x, size_y, size_z = dimstr[0], dimstr[1], dimstr[2]
    voxel_x, voxel_y, voxel_z = dimstr[3], dimstr[4], dimstr[5]
    return (size_x, size_y, size_z), (voxel_x, voxel_y, voxel_z)


def get_seg_path(exp_dir):
    seg_path = ""
    for dirpath, _, filenames in os.walk(exp_dir):
        for f in filenames:
            if "seg" in f:
                seg_path = os.path.abspath(os.path.join(dirpath, f))
    return seg_path


def add_scale_bar(axis):
    xdim = axis.get_xlim()[1]

    def find_scale_bar_length(xdim):
        length = xdim/8
        lengths = [2, 5, 10, 20, 50, 100, 200]
        return min(lengths, key=lambda x: abs(x - length))

    size = find_scale_bar_length(xdim)
    text = str(size)+r'$\mu$m'
    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
    bar = AnchoredSizeBar(axis.transData, size, text,
                          pad=0.5, loc=4, sep=5, borderpad=0.5, frameon=True)
    bar.patch.set(alpha=1, boxstyle='square')
    axis.add_artist(bar)
    axis.axis('off')


def rgb2id_array(rgb_array):
    """Return identifier array from rgb array"""
    b = rgb_array[:, :, 2].astype(np.uint64)
    g = 256 * rgb_array[:, :, 1].astype(np.uint64)
    r = 256 * 256 * rgb_array[:, :, 0].astype(np.uint64)
    return r + g + b


def id_array2rgb(id_array):
    """Return rgb array from identifier array"""
    rgb_array = np.zeros((id_array.shape[0], id_array.shape[1], 3), dtype=np.uint32)

    for row in range(id_array.shape[0]):
        for column in range(id_array.shape[1]):
            c = id_array[row, column]
            r = int(c / (256 * 256)) % 256
            g = int(c / (001 * 256)) % 256
            b = int(c / (001 * 001)) % 256
            rgb_array[row, column] = (r, g, b)
    return rgb_array


def structure_element(i, j, element_order=1):
    """ returns a structuring element of size 'element_order' centred at position i,j"""
    if element_order == 1:
        """
        [0,1,0]
        [1,1,1]
        [0,1,0]
        """
        i += 1
        j += 1

        return [(i, j + 1),
                (i - 1, j), (i, j), (i + 1, j),
                (i, j - 1)]

    if element_order == 2:
        """
        [0,0,1,0,0]
        [0,1,1,1,0]
        [1,1,1,1,1]
        [0,1,1,1,0]
        [0,0,1,0,0]
        """
        i += 2
        j += 2

        return [(i, j + 2),
                (i - 1, j + 1), (i, j + 1), (i + 1, j + 1),
                (i - 2, j), (i - 1, j), (i, j), (i + 1, j), (i + 2, j),
                (i - 1, j - 1), (i, j - 1), (i + 1, j - 1),
                (i, j - 2)]


def path2id_array(image_path):
    """return identifier array from image path"""
    seg_image = Image.open(image_path)
    seg_image_array_rgb = np.asarray(seg_image, dtype="uint64")
    cid_array = rgb2id_array(seg_image_array_rgb)
    return cid_array


def load_neighbours_dictionary_json(json_path):
    """return dictionary of lists of neighbours keyed by cid from json file path"""
    with open(json_path, "r") as input_handle:
        neighbours_dict_raw = json.load(input_handle)

    neighbours_dict = {int(k): v for k, v in neighbours_dict_raw.iteritems()}

    return neighbours_dict


def load_junction_dict(junction_path):
    """return dictionary of junction data keyed by junction id from json file path"""
    with open(junction_path, "r") as input_handle:
        junction_dict_raw = json.load(input_handle)

    junction_dict = {int(k): v for k, v in junction_dict_raw.iteritems()}

    return junction_dict


def path2image_array(image_path):
    """returns a greyscale numpy array from an image path"""
    raw_image = Image.open(image_path)
    raw_image_array = np.asarray(raw_image, dtype="uint64")
    if raw_image_array.ndim == 3:
        return np.amax(raw_image_array, 2)
    else:
        return raw_image_array


def generate_cell_outline_array(cid_array, color='black'):
    """generates a cell outline array for plotting heatmaps"""
    print "Generating cell outlines..."
    ele_ord = 1

    shape = (cid_array.shape[0], cid_array.shape[1], 4)
    outline_array = np.full(shape, [1, 1, 1, 0], dtype='float32')

    for i in range(cid_array.shape[0] - ele_ord * 2):
        for j in range(cid_array.shape[1] - ele_ord * 2):
            struc_element = [cid_array[i][j],
                             cid_array[i + 1][j],
                             cid_array[i][j + 1],
                             cid_array[i + 1][j + 1]]
            if len(set(struc_element)) > 1:
                if color == 'black':
                    outline_array[i][j] = [0, 0, 0, 1]
                elif color == 'white':
                    outline_array[i][j] = [1, 1, 1, 1]
                else:
                    raise ValueError("Color must be black or white")

    return outline_array


def path2outline_array(path):
    outline_array_fname = path + "_outline.npy"
    if os.path.exists(outline_array_fname):
        outline_array = np.load(outline_array_fname)
    else:
        seg_image = Image.open(path)
        seg_image_array_rgb = np.asarray(seg_image, dtype="uint64")
        cid_array = rgb2id_array(seg_image_array_rgb)
        outline_array = generate_cell_outline_array(cid_array)
        np.save(outline_array_fname, outline_array)
    return outline_array


class numpy_JSON_Encoder(json.JSONEncoder):
    """json encoder for encoding numpy data types"""

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(numpy_JSON_Encoder, self).default(obj)


def load_matches(matches_path):
    """loads and returns a matches file as a dictionary"""

    def get_list(raw_list):
        tl = raw_list.translate(None, ' []')
        ids = [int(i) for i in tl.split(',')]
        return ids

    with open(matches_path, 'r') as match_handle:
        lines = [l.strip() for l in match_handle.readlines()]

    splitlines = (l.split(':') for l in lines)
    return dict([(int(cfrom), get_list(cto)) for cfrom, cto in splitlines])


def load_junction_matches(matches_path):
    """return dictionary of junction data keyed by junction id from json file path"""
    with open(matches_path, "r") as input_handle:
        junction_matches_raw = json.load(input_handle)

    junction_matches = {int(k): v for k, v in junction_matches_raw.iteritems()}

    return junction_matches


def load_metadata(metadata_path):
    metadata = []
    with open(metadata_path, 'r') as file_handle:
        for row in file_handle:
            metadata.append(re.findall("[-+]?\d+[.]?\d*[eE]?[-+]?\d*", row))

    metadata_dict = {'vox_x': float(metadata[1][0]), 'vox_y': float(metadata[2][0]), 'vox_z': float(metadata[3][0])}
    return metadata_dict


def speak_complete(phrase='simulation finished'):
    """speak a phrase when called (mac specific)"""
    if platform.system() == 'Darwin':
        command = ['say', phrase]
        subprocess.check_call(command)
    else:
        print phrase


def rand_cmap(nlabels, type='bright', first_color_black=True, last_color_black=False, verbose=True):
    """
    from https://github.com/delestro/rand_cmap
    Creates a random colormap to be used together with matplotlib. Useful for segmentation tasks
    :param nlabels: Number of labels (size of colormap)
    :param type: 'bright' for strong colors, 'soft' for pastel colors
    :param first_color_black: Option to use first color as black, True or False
    :param last_color_black: Option to use last color as black, True or False
    :param verbose: Prints the number of labels and shows the colormap. True or False
    :return: colormap for matplotlib
    """
    from matplotlib.colors import LinearSegmentedColormap
    import colorsys
    import numpy as np

    if type not in ('bright', 'soft'):
        print ('Please choose "bright" or "soft" for type')
        return

    if verbose:
        print('Number of labels: ' + str(nlabels))

    # Generate color map for bright colors, based on hsv
    if type == 'bright':
        randHSVcolors = [(np.random.uniform(low=0.0, high=1),
                          np.random.uniform(low=0.2, high=1),
                          np.random.uniform(low=0.9, high=1)) for i in xrange(nlabels)]

        # Convert HSV list to RGB
        randRGBcolors = []
        for HSVcolor in randHSVcolors:
            randRGBcolors.append(colorsys.hsv_to_rgb(HSVcolor[0], HSVcolor[1], HSVcolor[2]))

        if first_color_black:
            randRGBcolors[0] = [0, 0, 0]

        if last_color_black:
            randRGBcolors[-1] = [0, 0, 0]

        random_colormap = LinearSegmentedColormap.from_list('new_map', randRGBcolors, N=nlabels)

    # Generate soft pastel colors, by limiting the RGB spectrum
    if type == 'soft':
        low = 0.6
        high = 0.95
        randRGBcolors = [(np.random.uniform(low=low, high=high),
                          np.random.uniform(low=low, high=high),
                          np.random.uniform(low=low, high=high)) for i in xrange(nlabels)]

        if first_color_black:
            randRGBcolors[0] = [0, 0, 0]

        if last_color_black:
            randRGBcolors[-1] = [0, 0, 0]
        random_colormap = LinearSegmentedColormap.from_list('new_map', randRGBcolors, N=nlabels)

    # Display colorbar
    if verbose:
        from matplotlib import colors, colorbar
        from matplotlib import pyplot as plt
        fig, ax = plt.subplots(1, 1, figsize=(15, 0.5))

        bounds = np.linspace(0, nlabels, nlabels + 1)
        norm = colors.BoundaryNorm(bounds, nlabels)

        cb = colorbar.ColorbarBase(ax, cmap=random_colormap, norm=norm, spacing='proportional', ticks=None,
                                   boundaries=bounds, format='%1i', orientation=u'horizontal')

    return random_colormap