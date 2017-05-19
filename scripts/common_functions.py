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
    rgb_array = np.zeros((id_array.shape[0], id_array.shape[1], 3))

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


def experimental_distributions():
    mombach_sativum = [0, 0, 0, 0, 0.045024633, 0.2246639, 0.5251844, 0.15445487, 0.045634203, 9.16E-04]
    mombach_cepa = [0, 0, 0, 0, 0.030464832, 0.19987738, 0.5895831, 0.14839482, 0.033524692, 9.16E-04]
    mombach_attenuata = [0, 0, 0, 0, 0.0608693, 0.19823812, 0.46402153, 0.22525957, 0.030464832, 0]
    mombach_arborescens = [0, 0, 0, 0, 0.03385827, 0.23682013, 0.51975787, 0.18229046, 0.02589178, 0]
    mombach_anthurium = [0, 0, 0, 0, 0.06630694, 0.2658155, 0.37409183, 0.21137506, 0.06599349, 0.01186091]
    ross_arabidopsis = [0, 0, 0, 0.069386005, 0.287735939, 0.343747139, 0.201349974, 0.073713303, 0.019390345,
                        0.003937006, 0.00065876, 8.34E-05, 1.03E-05]
    gibson_drosophila = [0, 0, 0, 0, 0.029033862, 0.27950263, 0.45671737, 0.20122315, 0.031333886, 0.001563567]
    gibson_xenopus = [0, 0, 0, 0, 0.029033862, 0.28972605, 0.4290362, 0.18163446, 0.048798237, 0.007954422]
    gibson_cucumber = [0, 0, 0, 0, 0.029033862, 0.25055113, 0.4741881, 0.22423096, 0.030065734, 0.001146958]

    return {
        "mombach_sativum": mombach_sativum,
        "mombach_cepa": mombach_cepa,
        "mombach_attenuata": mombach_attenuata,
        "mombach_arborescens": mombach_arborescens,
        "mombach_anthurium": mombach_anthurium,
        "ross_arabidopsis": ross_arabidopsis,
        "gibson_drosophila": gibson_drosophila,
        "gibson_xenopus": gibson_xenopus,
        "gibson_cucumber": gibson_cucumber
    }

if __name__ == "__main__":
    path = '/Users/carterr/Dropbox/Postdoc_stuff/pcgeometry/data/raw/3002_PD/microscope_metadata/T00.txt'
    print load_metadata(path)
