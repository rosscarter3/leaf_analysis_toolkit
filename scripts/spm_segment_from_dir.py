import os
from sys import platform
import argparse
import shutil
import subprocess


def get_im_paths(exp_dir):
    im_path = []
    for dirpath, _, filenames in os.walk(exp_dir):
        for f in filenames:
            if f.endswith('_proj-g3d.png') or f.endswith('_proj-g3d_rev.png'):
                im_path.append(os.path.abspath(os.path.join(dirpath, f)))
    return im_path


def get_latest_seg_paths(exp_dir):
    fwd_dirs = []
    for dirpath, dirnames, _ in os.walk(os.path.join(exp_dir, "segmented")):
        for directory in dirnames:
            fwd_dirs.append(os.path.abspath(os.path.join(dirpath, directory)))
        break
    rev_dirs = []
    for dirpath, dirnames, _ in os.walk(os.path.join(exp_dir, "segmented_rev")):
        for directory in dirnames:
            rev_dirs.append(os.path.abspath(os.path.join(dirpath, directory)))
        break

    fwd_dirs.sort()
    rev_dirs.sort()

    try:
        fwd_image_path = os.path.join(fwd_dirs[-1], "CellstateFinal", "00000.png")
    except:
        fwd_image_path = None
    try:    
        rev_image_path = os.path.join(rev_dirs[-1], "CellstateFinal", "00000.png")
    except:
        rev_image_path = None

    return fwd_image_path, rev_image_path


def main():
    exp_dir = os.path.abspath(args.exp_dir)

    script_path = os.path.realpath(__file__)

    if platform == "linux" or platform == "linux2":
        bin_path = os.path.join(os.path.dirname(script_path), "../bin/spm2D_1_0_1-reversed-lin_bld")
    elif platform == "darwin":
        bin_path = os.path.join(os.path.dirname(script_path), "../bin/spm2D_1_0_1-reversed-mac_bld")
    else:
        print "system not supported"
        return 0

    par_path = os.path.join(os.path.dirname(script_path), "../parameter_files/spm2d.par")

    im_files = get_im_paths(exp_dir)

    for im_path in im_files:
        if im_path.endswith("_proj-g3d_rev.png"):
            output_dir = os.path.join(exp_dir, "segmented_rev")
        else:
            output_dir = os.path.join(exp_dir, "segmented")
        segment_process = bin_path + " " + output_dir + " " + par_path + " " + im_path
        subprocess.call(segment_process, shell=True)

    fwd_path, rev_path = get_latest_seg_paths(exp_dir)
    fwd_dest_path = os.path.join(exp_dir, "fwd_seg.png")
    rev_dest_path = os.path.join(exp_dir, "rev_seg.png")

    if fwd_path is not None:
        shutil.copy(fwd_path, fwd_dest_path)
    if rev_path is not None:    
        shutil.copy(rev_path, rev_dest_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="submits the two gauss projected images for segmentation \
                                                 (serially), uses par file from ../parameter_files,  \
                                                 detects system and uses relevant build from ../bin")
    parser.add_argument('exp_dir', help="path to the image directory")

    args = parser.parse_args()

    main()
