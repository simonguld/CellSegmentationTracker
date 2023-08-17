## Author: Simon Guldager Andersen

## Imports:
import os, sys, glob
import tifftools
import json
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
from cycler import cycler

from utils import trackmate_xml_to_csv, merge_tiff, get_imlist, prepend_text_to_file, search_and_modify_file

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





##TODO:

##ESSENTIAL (DO FIRST)'

# mod cellpose (forsøg dynamisk) og verificer.
#--> verificer class og readme til at tage de argumenter.
# Træn modeller en for en og: 
# kvalitetsforskel på store/små billedeR?
# gem billeder af gode/dårlige sekventeringer +parametre
# dokumenter løbende
# efterprøv at de virker

# add cellprob, flow_thres til docs, anaconda

# fix cellmask color if possible(se på features/featureutils.java)
# add docstring to class
# #### crop tif img series????


# resize img -- inkl crop tif image series
## FIX utiils --> cellsegmentationtracker.utils.


## NONESSENTIAL (IF TIME)
# make naming of xml and csv files more flexible
# make roubst under other image formats?


class CellSegmentationTracker:

    """
    Documentation for CellSegmentationTracker class.
    """

    def __init__(self, imagej_filepath, cellpose_python_filepath, image_folder_path = None, xml_path = None, output_folder_path = None,
                  use_model = 'CYTO', custom_model_path = None, show_segmentation = True, cellpose_dict = dict(), trackmate_dict = dict(),):

        self.img_folder = image_folder_path
        self.img_path = None
        self.xml_path = xml_path
        self.imagej_filepath = imagej_filepath
        self.cellpose_python_filepath = cellpose_python_filepath   
        self.custom_model_path = custom_model_path
        self.__current_path = os.getcwd()
        self.__parent_dir = os.path.abspath(os.path.join(self.__current_path, os.pardir))
        self.__fiji_folder_path = os.path.dirname(self.imagej_filepath)
        self.__cellpose_folder_path = os.path.join(os.path.dirname(self.cellpose_python_filepath), 'Lib', 'site-packages', 'cellpose')
        
        # Set output folder path to image folder path, if not provided
        self.output_folder = output_folder_path
        if self.output_folder is None:
            if self.img_folder is not None:
                self.output_folder = self.img_folder if os.path.isdir(self.img_folder) else os.path.dirname(self.img_folder)
            else:
                self.output_folder = self.__current_path
        else:
            self.output_folder = output_folder_path

        self.use_model = use_model
        self.show_segmentation = show_segmentation
        self.cellpose_dict = cellpose_dict
        self.trackmate_dict = trackmate_dict
        self.jython_dict = {}
        self.csv = None

        # Extract or set flow and cell probability threshold
        try:
            self.flow_threshold = self.cellpose_dict['FLOW_THRESHOLD']
            del self.cellpose_dict['FLOW_THRESHOLD']
        except:
            self.flow_threshold = 0.4
        try:
            self.cellprob_threshold = self.cellpose_dict['CELLPROB_THRESHOLD']
            del self.cellpose_dict['CELLPROB_THRESHOLD']
        except:
            self.cellprob_threshold = 0.0
        
        self.pretrained_models = ['CYTO', 'CYTO2', 'NUCLEI', 'EPI1']
        self.pretrained_models_paths = [os.path.join(self.__cellpose_folder_path, 'models', 'cyto_0'), \
                                        os.path.join(self.__cellpose_folder_path, 'models', 'cyto_1'),\
                                        os.path.join(self.__cellpose_folder_path, 'models', 'nuclei'), \
                                            os.path.join(self.__parent_dir, 'models', 'epi1')]
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

        
    
        # If no xml file is provided (or present in image folder), cell segmentation, tracking and generation of an xml file is performed
        if self.xml_path is None:
            # Prepare images for cellpose segmentation
            self.__prepare_images()

            # Ensure that xml file has not already been created
            if os.path.isfile(self.img_path.strip(".tif") + ".xml"):
                pass
            else:
                # Generate dictionary for jython script
                self.__generate_jython_dict()

                # If not already done, modify cellpose source code to allow for modification of flow and cell probability threshold
                if not os.path.isfile(os.path.join(self.__cellpose_folder_path, 'params.txt')):
                    self.__modify_cellpose()

                # Save flow and cell probability threshold to txt file
                np.savetxt(os.path.join(self.__cellpose_folder_path, 'params.txt'), \
                           np.array([self.flow_threshold, self.cellprob_threshold]))
                
                
                # Save dictionary to json file
                with open(os.path.join(self.__current_path,"jython_dict.json"), 'w') as fp:
                    json.dump(self.jython_dict, fp)

                # Run jython script
                self.__run_jython_script()

                # Reset flow and cell probability threshold to default values after run
                np.savetxt(os.path.join(self.__cellpose_folder_path, 'params.txt'), \
                           np.array([0.4, 0.0]))
                
            # Set xml path to image path
            self.xml_path = self.img_path.strip(".tif") + ".xml"

        self.csv = trackmate_xml_to_csv(self.xml_path, get_tracks = True)


    # Private methods

    def __prepare_images(self):
        """
        Prepare images for cellpose segmentation.
        """
        # If no image folder is provided, raise error
        if self.img_folder is None:
            raise OSError("No image folder nor xml file provided!")
        # If image folder is provided, check if it is a folder or a .tif file
        elif self.img_folder[-4:] != ".tif" and not os.path.isdir(self.img_folder):
            raise OSError("Image folder path must either be a folder or a .tif file!")
        # If .tif file is provided instead of folder, handle it
        elif self.img_folder[-4:] == ".tif":
            self.img_path = self.img_folder
            self.img_folder = os.path.dirname(self.img_folder)
        else:
            im_list = get_imlist(self.img_folder, '.tif')
            # If more than one .tif file is provided, merge them into one
            if len(im_list) > 1:
                # Make sure that such a file has not already been created    
                self.img_path = os.path.join(self.img_folder, "merged.tif")
                if os.path.isfile(self.img_path):
                    os.remove(self.img_path)
                merge_tiff(self.img_folder, img_out_name = os.path.join(self.img_folder, "merged.tif"))
                print("Merged tif files into one file: ", self.img_path)
            else:
                self.img_path = im_list[0]
                print("Using tif file: ", self.img_path)
        return

    def __modify_cellpose(self):
        """
        Modify cellpose source code to allow for modification of flow and cell probability threshold.
        """
        # Define paths to files to modify
        paths_to_modify = [os.path.join(self.__cellpose_folder_path, name) for name in ['models.py', 'dynamics.py', 'cli.py']]

        # Define search and replace strings
        search_strings = ["'--flow_threshold', default=0.4", "'--cellprob_threshold', default=0", \
                      'flow_threshold=0.4', 'cellprob_threshold=0.0',  \
                      'flow_threshold= 0.4', 'flow_threshold = 0.4', \
                        'cellprob_threshold= 0.0', 'cellprob_threshold = 0.0' ]
        replace_strings = ["'--flow_threshold', default=flow_param", "'--cellprob_threshold', default=cellprop_param", \
                          'flow_threshold=flow_param', 'cellprob_threshold=cellprop_param', \
                            'flow_threshold= flow_param', 'flow_threshold = flow_param', \
                            'cellprob_threshold= cellprop_param', 'cellprob_threshold = cellprop_param' ]

        # Define message to prepend to file to allow for modification of flow and cell probability threshold
        msg = f"\nimport numpy as np\nimport os\ndir = os.path.dirname(os.path.abspath(__file__))\nflow_param, cellprop_param = np.loadtxt(os.path.join(dir, 'params.txt')) \n\n"  

        # Modify files
        for p in paths_to_modify:
            prepend_text_to_file(p, msg)
            for search_string, replace_string in zip(search_strings, replace_strings):
                search_and_modify_file(p, search_string, replace_string)
        print("Cellpose source code files have been modified to allow for modification of flow and cell probability threshold.")
        return

    def __generate_jython_dict(self):
        """
        Generate the dictionary holding the parameters needed to run the jython script.
        """
        # Specify paths and whether to show segmentation
        self.jython_dict.update({"IMG_PATH": self.img_path, "FIJI_FOLDER_PATH": self.__fiji_folder_path, \
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
        return

    def __run_jython_script(self):
        """
        Run the jython script by calling the command line.
        """
        # Change directory to fiji folder
        os.chdir(self.__fiji_folder_path)
        # Get name of executable
        executable = list(os.path.split(self.imagej_filepath))[-1]
        # Set path to jython script
        jython_path = os.path.join(self.__current_path, "jython_cellpose.py")

        if self.show_segmentation: 
            command = executable + ' --ij2 --console --run "' + jython_path + '"'
            os.system(command)
        else:
            self.xml_path = self.img_path.strip(".tif") + ".xml"
        
            # To run headlessly, call jython script as a subprocess
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
            
        # Change directory back to current one
        os.chdir(self.__current_path)    
        return
    

    # Public methods

    def save_csv(self, name = 'TrackMate.csv'):
        path_out = os.path.join(self.output_folder, name)
        self.csv.to_csv(path_out, sep = ',') 
        print("Saved csv file to: ", path_out)
        return

        

def main():
    image_path = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\valeriias mdck data for simon\\24.08.22_698x648\\im0.tif"
    input_directory = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\valeriias mdck data for simon\\24.08.22_698x648"
    output_directory = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\valeriias mdck data for simon\\24.08.22_698x648"
    model_directory = 'C:\\Users\\Simon Andersen\\miniconda3\\envs\\cellpose\\lib\\site-packages\\cellpose\\models\\cyto_0'
    cellpose_python_filepath = 'C:\\Users\\Simon Andersen\\miniconda3\\envs\\cellpose\\python.exe'
    imj_path = "C:\\Users\\Simon Andersen\\Fiji.app\\ImageJ-win64.exe"
    image_path = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\16.06.23_stretch_data_698x648\\im0.tif"

    show_output = True
    cellpose_dict = {'USE_GPU': True, 'FLOW_THRESHOLD': 10, 'CELLPROB_THRESHOLD': -5}

    # xml_path, outfolder, use_model = cyto2,nuclei, 
    xml_path = 'C:\\Users\\Simon Andersen\\Projects\\Projects\\CellSegmentationTracker\\resources\\20.09.xml'

    cst = CellSegmentationTracker(imj_path, cellpose_python_filepath, image_path,
                                  show_segmentation=show_output, cellpose_dict=cellpose_dict, use_model='CYTO2')

    cst.save_csv()
    print("you made it bro")

if __name__ == '__main__':
    main()