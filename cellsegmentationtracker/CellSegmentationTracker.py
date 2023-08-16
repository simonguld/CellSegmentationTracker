## Author: Simon Guldager Andersen

## Imports:
import os, sys, glob
import tifftools, json
import platform
import time
from subprocess import Popen, PIPE

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

from utils import trackmate_xml_to_csv, merge_tiff

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




## necessary args:

# if working pipeline:
# path to tif images, recropping argument, trackmate arguments, cellpose arguments, choice of model, save_img = True, 

# if not working pipeline
# path to (possibly) xml and definityely csv spot track edge files

##TODO:

##ESSENTIAL (DO FIRST)


## decide which members should be private

## NONESSENTIAL (IF TIME)
# resize img
# make naming of xml and csv files more flexible
# make roubst under other image formats?


class CellSegmentationTracker:

    def __init__(self, imagej_filepath, cellpose_python_filepath, image_folder = None, xml_path = None, output_folder = None,
                  use_model = 'CYTO', custom_model_path = None, show_segmentation = True, cellpose_dict = dict(), trackmate_dict = dict(),):

        self.img_folder = image_folder
        self.img_path = None
        self.xml_path = xml_path
        self.imagej_filepath = imagej_filepath
        self.fiji_folder_path = os.path.abspath(os.path.join(self.imagej_filepath, os.pardir))
        self.cellpose_python_filepath = cellpose_python_filepath
        self.use_model = use_model
        self.custom_model_path = custom_model_path
        self.output_folder = output_folder
        self.current_path = os.getcwd()
        self.parent_dir = os.path.abspath(os.path.join(self.current_path, os.pardir))

        self.show_segmentation = show_segmentation
        self.cellpose_dict = cellpose_dict
        self.trackmate_dict = trackmate_dict
        self.use_model = use_model
        self.jython_dict = {}

        self.pretrained_models = ['CYTO', 'CYTO2', 'NUCLEI', 'EPI1']
        self.pretrained_models_paths = [os.path.join(self.parent_dir, 'models', 'models', 'cyto_0'), os.path.join(self.parent_dir, 'models', 'models', 'cyto_1'),\
                                        os.path.join(self.parent_dir, 'models', 'models', 'nuclei'), os.path.join(self.parent_dir, 'models', 'epi1')]
        if self.custom_model_path is not None:
            self.use_model = 'CUSTOM'
        # If custom model is not provided, find the path of the pretrained model
        if self.custom_model_path is None and use_model not in self.pretrained_models:
            raise ValueError("No custom model path provided, and use_model not in pretrained models: ", self.pretrained_models)
        elif self.custom_model_path is None and use_model in self.pretrained_models:
            idx = self.pretrained_models.index(use_model)
            self.custom_model_path = self.pretrained_models_paths[idx]
        

        self.cellpose_default_values = {
            'TARGET_CHANNEL' : 0,
            'OPTIONAL_CHANNEL_2': 0,
            'CELLPOSE_PYTHON_FILEPATH': self.cellpose_python_filepath,
            'CELLPOSE_MODEL': self.use_model,
            'CELLPOSE_MODEL_FILEPATH': self.custom_model_path,
            'CELL_DIAMETER': 0.0,
            'USE_GPU': False,
            'SIMPLIFY_CONTOURS': True
            }
        self.trackmate_default_values = {'LINKING_MAX_DISTANCE': 15.0,
                                         'GAP_CLOSING_MAX_DISTANCE': 15.0,
                                         'MAX_FRAME_GAP': 2,
                                         'ALLOW_TRACK_SPLITTING': True,
                                         'ALLOW_TRACK_MERGING': True,
             }

        self.csv = None
        
        # If no xml file is provided (or present in image folder), cell segmentation, tracking and generation of an xml file is performed
        if self.xml_path is None:
            # If no image folder is provided, raise error
            if self.img_folder is None:
                raise ValueError("No image folder nor xml file provided!")

            if self.img_folder[-4:] == ".tif":
                self.img_path = self.img_folder
                self.img_folder = os.path.abspath(os.path.join(self.img_folder, os.pardir))
            # If more than one .tif file is provided, merge them into one
            else:
                # Make so that such a file has not already been created    
                self.img_path = os.path.join(self.img_folder, "merged.tif")
                if os.path.isfile(self.img_path):
                    os.remove(self.img_path)
                merge_tiff(self.img_folder, img_out_name = os.path.join(self.img_folder, "merged.tif"))
                print("Merged tif files into one file: ", self.img_path)

            # Ensure that xml file has not already been created
            if os.path.isfile(self.img_path.strip(".tif") + ".xml"):
                pass
            else:
                # Generate dictionary for jython script, if not already provided
                self.jython_dict.update({"IMG_PATH": self.img_path, "FIJI_FOLDER_PATH": self.fiji_folder_path, \
                                    "CELLPOSE_PYTHON_FILEPATH": self.cellpose_python_filepath, "CUSTOM_MODEL_PATH": self.custom_model_path,
                                    "OUTPUT_PATH": self.output_folder, "SHOW_SEGMENTATION": self.show_segmentation})
                
                # Specify Cellpose and TrackMate parameters not set by user
                if self.cellpose_dict == {}:
                    self.cellpose_dict = self.cellpose_default_values
                else:
                    for key in self.cellpose_default_values.keys():
                        if key not in self.cellpose_dict.keys():
                            self.cellpose_dict[key] = self.cellpose_default_values[key]
                if self.trackmate_dict == {}:
                    self.trackmate_dict = self.trackmate_default_values
                else:
                    for key in self.trackmate_default_values.keys():
                        if key not in self.trackmate_dict.keys():
                            self.trackmate_dict[key] = self.trackmate_default_values[key]

                self.jython_dict.update({'CELLPOSE_DICT': self.cellpose_dict})
                self.jython_dict.update({'TRACKMATE_DICT': self.trackmate_dict})


                # Save dictionary to json file
                with open(os.path.join(self.current_path,"jython_dict.json"), 'w') as fp:
                    json.dump(self.jython_dict, fp)

                print("self.custom_model_path: ", self.custom_model_path)
                print("self.use_model: ", self.use_model)
                # Run jython script
                os.chdir(self.fiji_folder_path)
                executable = list(os.path.split(self.imagej_filepath))[-1]
                jython_path = os.path.join(self.current_path, "jython_cellpose.py")

                if self.show_segmentation: 
                    command = executable + ' --ij2 --console --run "' + jython_path + '"'
                    os.system(command)
                else:
                    self.xml_path = self.img_path.strip(".tif") + ".xml"
                
                    # call jython script as a subprocess
                    pipe = Popen([executable,'--ij2','--headless', '--run', f"{jython_path}"], stdout=PIPE)
                    # wait for xml file to be created
                    try:
                        while not os.path.isfile(self.xml_path):
                            time.sleep(3)
                            continue
                    except: 
                        KeyboardInterrupt
                        sys.exit(0)
                    # kill subprocess
                    pipe.kill()
                    # print output
                    print(pipe.communicate()[0].decode('ascii'))
  
            # Set xml path to image path
            self.xml_path = self.img_path.strip(".tif") + ".xml"

        # Generate csv file from xml file
        self.xml = self.img_folder.strip(".tif") + ".xml"
  
        self.csv = trackmate_xml_to_csv(self.xml_path, get_tracks = True)

    def save_csv(self):
        if self.output_folder is not None:
            self.csv.to_csv(os.path.join(self.output_folder, 'TrackMate.csv'), sep = ',') 
        else:
            self.csv.to_csv(os.path.join(self.img_folder, 'TrackMate.csv'), sep = ',')
        return


        

def main():
    image_path = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\valeriias mdck data for simon\\24.08.22_698x648\\im0.tif"
    input_directory = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\valeriias mdck data for simon\\24.08.22_698x648"
    output_directory = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\valeriias mdck data for simon\\24.08.22_698x648"
    model_directory = 'C:\\Users\\Simon Andersen\\miniconda3\\envs\\cellpose\\lib\\site-packages\\cellpose\\models\\cyto_0'
    cellpose_python_filepath = 'C:\\Users\\Simon Andersen\\miniconda3\\envs\\cellpose\\python.exe'
    imj_path = "C:\\Users\\Simon Andersen\\Fiji.app\\ImageJ-win64.exe"

    show_output = False
    cellpose_dict = {'USE_GPU': True}

    cst = CellSegmentationTracker(imj_path, cellpose_python_filepath, image_path, \
                                  show_segmentation=show_output, cellpose_dict=cellpose_dict, use_model='EPI1')

    cst.save_csv()
    print("you made it bro")

if __name__ == '__main__':
    main()