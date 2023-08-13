## Imports:
import os, sys, glob
import tifftools, json
import platform
import numpy as np
import seaborn as sns
import pandas as pd
import xml.etree.ElementTree as et

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rcParams
from matplotlib.patches import Ellipse
from cycler import cycler
from skimage.io import imread

#from pytrackmate import trackmate_peak_import 

## Change directory to current one
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)

## Set plotting style and print options
sns.set_theme()
sns.set_style("darkgrid")
sns.set_context("paper") #Possible are paper, notebook, talk and poster

d = {'lines.linewidth': 2, 'axes.titlesize': 18, 'axes.labelsize': 18, 'xtick.labelsize': 12, 'ytick.labelsize': 12,\
     'legend.fontsize': 15, 'font.family': 'serif', 'figure.figsize': (9,9)}
d_colors = {'axes.prop_cycle': cycler(color = ['teal', 'navy', 'coral', 'plum', 'purple', 'olivedrab',\
         'black', 'red', 'cyan', 'yellow', 'khaki','lightblue'])}
rcParams.update(d)
rcParams.update(d_colors)
np.set_printoptions(precision = 5, suppress=1e-10)


def get_imlist(path, format = '.jpg'):

    """
    returns a list of filenames for all png images in a directory
    """
    return [os.path.join(path,f) for f in os.listdir(path) if f.endswith(format)]

def merge_tiff(im_dir, img_out_name):
    img_in_names = get_imlist(im_dir, format = '.tif')
    if not img_in_names:
        print("No image!!!")
        return -1
    print("File names: {:}".format(img_in_names))
    t = tifftools.read_tiff(img_in_names[0])
    for img_in_name in img_in_names[1:]:
        t["ifds"].extend(tifftools.read_tiff(img_in_name)["ifds"])

    if os.path.isfile(img_out_name):
        os.unlink(img_out_name)
    tifftools.write_tiff(t, img_out_name)

    imgo = imread(img_out_name)
    imgs = [imread(e) for e in img_in_names]

    print("Output file details: {:}, {:}".format(imgo.shape, imgo.dtype))
    print("Input files details: {:}".format([(img.shape, img.dtype) for img in imgs]))

    if len(img_in_names) in (3, 4):  # @TODO:
        imgo = imgo.transpose((2, 0, 1))

## necessary args:

# if working pipeline:
# path to tif images, recropping argument, trackmate arguments, cellpose arguments, choice of model, save_img = True, 

# if not working pipeline
# path to (possibly) xml and definityely csv spot track edge files

##TODO:
# extend args so as to be able to choose model
# fix hele operating system tingen ift imJ == ekstraher subscript fra imJ-filepath
# omskriv jytoh script, så den dicten som input


class CellSegmentationTracker:
    def __init__(self, imagej_filepath, cellpose_python_filepath, img_folder = None, xml_path = None, output_folder = None,
                  use_model = 'CYTO', custom_model_path = None, show_segmentation = False, cellpose_dict = None, trackmate_dict = None):

        self.img_path = img_folder
        self.xml_path = xml_path
        self.fiji_folder_path = imagej_filepath
        self.cellpose_python_filepath = cellpose_python_filepath
        self.custom_model_path = custom_model_path
        self.output_folder = output_folder
        self.current_path = os.getcwd()

        self.show_segmentation = show_segmentation
        self.cellpose_dict = cellpose_dict
        self.trackmate_dict = trackmate_dict
        self.use_model = use_model
        self.jython_dict = {}

        self.csv = None
        
        if self.xml_path is None:
            im_list = get_imlist(self.img_folder, format = '.tif')
            if len(im_list) == 0:
                print("The image folder contains no .tif images!")
                sys.exit()
            # If folder contains more than 1 tif file, merge them
            elif len(im_list) > 1:
                merge_tiff(self.img_folder, img_out_name = os.path.join(self.img_folder, "merged.tif"))
                self.img_path = os.path.join(self.img_folder, "merged.tif")
            
            # Generate dictionary for jython script
            self.jython_dict.update({"IMG_PATH": self.img_path, "FIJI_FOLDER_PATH": self.fiji_folder_path, \
                                "CELLPOSE_PYTHON_FILEPATH": self.cellpose_python_filepath, "CUSTOM_MODEL_PATH": self.custom_model_path,
                                 "OUTPUT_PATH": self.output_folder, "SHOW_SEGMENTATION": self.show_segmentation})
            if self.cellpose_dict is not None:
                if 'TARGET_CHANNEL' not in self.cellpose_dict.keys():
                    self.cellpose_dict['TARGET_CHANNEL'] = 0
                if 'OPTIONAL_CHANNEL' not in self.cellpose_dict.keys():
                    self.cellpose_dict['OPTIONAL_CHANNEL'] = 0
                if 'CELLPOSE_MODEL' not in self.cellpose_dict.keys():
                    self.cellpose_dict['CELLPOSE_MODEL'] = use_model
                if 'CELL_DIAMETER' not in self.cellpose_dict.keys():
                    self.cellpose_dict['CELL_DIAMETER'] = 30.0
                if 'USE_GPU' not in self.cellpose_dict.keys():
                    self.cellpose_dict['USE_GPU'] = False
                if 'SIMPLIFY_CONTOURS' not in self.cellpose_dict.keys():
                    self.cellpose_dict['SIMPLIFY_CONTOURS'] = True

            if self.trackmate_dict is not None:
                if 'LINKING_MAX_DISTANCE' not in self.trackmate_dict.keys():
                    self.trackmate_dict['LINKING_MAX_DISTANCE'] = 15.0
                if 'GAP_CLOSING_MAX_DISTANCE' not in self.trackmate_dict.keys():
                    self.trackmate_dict['GAP_CLOSING_MAX_DISTANCE'] = 15.0
                if 'MAX_FRAME_GAP' not in self.trackmate_dict.keys():
                    self.trackmate_dict['MAX_FRAME_GAP'] = 2
                if 'ALLOW_TRACK_SPLITTING' not in self.trackmate_dict.keys():
                    self.trackmate_dict['ALLOW_TRACK_SPLITTING'] = True
                if 'ALLOW_TRACK_MERGING' not in self.trackmate_dict.keys():
                    self.trackmate_dict['ALLOW_TRACK_MERGING'] = True

            self.jython_dict.update(self.cellpose_dict)
            self.jython_dict.update(self.trackmate_dict)

            # Save dictionary to json file
            with open(os.path.join(self.output_folder,"jython_dict.json"), 'w') as fp:
                json.dump(self.jython_dict, fp)


            # Run jython script
            os.chdir(self.fiji_folder_path)

            


        
    def print_path1(self):
        print("Path 1:", self.path1)
        
    def print_path2(self):
        print("Path 2:", self.path2)