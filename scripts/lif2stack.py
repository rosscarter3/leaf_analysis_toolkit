import os
import argparse
import subprocess
import ast
import xml.etree.ElementTree
from xml import etree as et


def parse_xml_metadata_2(xml_string, image_id, array_order='xyz'):
    array_order = array_order.upper()
    spatial_array_order = [c for c in array_order if c in 'XYZ']
    size_tags = ['Size' + c for c in array_order]
    res_tags = ['PhysicalSize' + c for c in spatial_array_order]
    #xml_string.encode('utf8', 'ignore')
    metadata_root = et.ElementTree.fromstring(xml_string)

    for child in metadata_root:
        if child.tag.endswith('Image') and str(image_id) in child.attrib['ID']:
            name = child.attrib['Name']
            for grandchild in child:
                if grandchild.tag.endswith('Pixels'):
                    att = grandchild.attrib
                    try:
                        size = tuple([int(att[t]) for t in size_tags])
                    except KeyError:
                        size = tuple([0, 0, 0])
                    try:
                        resolution = tuple([float(att[t]) for t in res_tags])
                    except KeyError:
                        resolution = tuple([0, 0, 0])
    return name, size, resolution


def write_dims(sizes, resolutions, output_dir):
    """writes the image dimensions and voxel dimensions to
    dims.txt"""
    with open(os.path.join(output_dir, 'dims.txt'), 'w') as file_handle:
        file_handle.write("%s,%s,%s,%s,%s,%s" % (sizes[0],
                                                 sizes[1],
                                                 sizes[2],
                                                 resolutions[0],
                                                 resolutions[1],
                                                 resolutions[2]))


def detect_num_images(lif_path):
    """
    returns the number of image stacks contained within a .lif file
    :param lif_path: path to the lif file of interest
    :return: num_ims: number of stacks contained in the provided lif file
    """

    info_process = ["showinf", "-nopix", lif_path]
    info_string = subprocess.check_output(info_process)
    num_ims = 0
    for l in info_string.split("\n"):
        if l.startswith("Series count"):
            num_ims = int(l[-1])
    return num_ims


def main():
    lif_path = args.lif_path
    
    metadata_process = ["showinf", "-nopix", "-omexml-only", lif_path]
    xml_string = subprocess.check_output(metadata_process)
    xml_string.decode('utf8', errors='ignore')

    if args.im_list is not None:
        im_list = ast.literal_eval(args.im_list)
        print "extracting ", im_list
    else:
        num_ims = detect_num_images(lif_path)
        im_list = range(num_ims)
        print "extracting all"
    
    for im_number in im_list:
        name, size, resolution = parse_xml_metadata_2(xml_string, im_number)
        print "extracting", name
        output_dir = "./proj_%s_%s_%s" % (lif_path, im_number, name.replace(" ", "_"))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        stack_output_dir = os.path.join(output_dir, "stack")
        if not os.path.exists(stack_output_dir):
            os.makedirs(stack_output_dir)

        bfconvert_process = ["bfconvert", "-series", str(im_number), lif_path,
                             os.path.join(stack_output_dir, "stack%z.tiff")]
        subprocess.call(bfconvert_process)

        write_dims(size, resolution, output_dir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="extracts the z slices from a particular stack in a "
                                                 " .lif file, requires bioformats command line tools")
    parser.add_argument('lif_path', help="Path to lif file")
    parser.add_argument("-l", '--im_list', help="list of images to extract and project")

    args = parser.parse_args()

    main()
