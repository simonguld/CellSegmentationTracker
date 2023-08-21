## Author: Simon Guldager Andersen

## Imports:
import os
import sys
import json
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

##ESSENTIAL (DO FIRST)


# Træn sidste model + afslut doku
# TEST: Prøv modeller af på de forskellige datasæt. virker DET?
# gem billeder af dårlige segmenteringer 
# skriv valll tas




## FOR LATER:
# TODO:
# grundig gennemlæsning af class for at checke efter fejl + relevante tests


# SPEND A LITTLE TIME ON:
# make show_tracks opt?
# fix cellmask color if possible(se på features/featureutils.java)


# CONSIDER:
# make imp. faster by 1) not loading all features to csv and/or 2) not calculating everything in trackmate 3) resizing

# BEFORE SENDING TO VALERIIA:
# add docstring to class
# document methods(should csv_dfs be methods? Otherwise throw them in as attributes)
# finish readme and example notebook
# make sure the package works INCLUDING FIX utiils --> cellsegmentationtracker.utils.

## if running full size img are super slow, consider implementing resize after all

## NONESSENTIAL (IF TIME)
# include the resize/splitting stuff (Resizing must preserve rel. proportions. ALSO: alters the observables. Take into account)
# include option to give pixel=length, frame=time 
# train nuclei data. NB: 1 bad img
# make naming of xml and csv files more flexible
# make roubst under other image formats?
# make robust under color, multichannel etc?


class CellSegmentationTracker:

    """
    Documentation for CellSegmentationTracker class.
    """

    def __init__(self, imagej_filepath = None, cellpose_python_filepath = None, image_folder_path = None, xml_path = None, output_folder_path = None,
                  use_model = 'CYTO', custom_model_path = None, show_segmentation = True, cellpose_dict = dict(), trackmate_dict = dict()):

        self.img_folder = image_folder_path
        self.img_path = None
        self.xml_path = xml_path
        self.imagej_filepath = imagej_filepath
        self.cellpose_python_filepath = cellpose_python_filepath   
        self.custom_model_path = custom_model_path
        self.__current_path = os.getcwd()
        self.__parent_dir = os.path.abspath(os.path.join(self.__current_path, os.pardir))
        self.__fiji_folder_path = os.path.dirname(self.imagej_filepath) if self.imagej_filepath is not None else None

        if self.cellpose_python_filepath is not None:
            self.__cellpose_folder_path = os.path.join(os.path.dirname(self.cellpose_python_filepath), 'Lib', 'site-packages', 'cellpose')
            self.pretrained_models_paths = [os.path.join(self.__cellpose_folder_path, 'models', 'cyto_0'), \
                                        os.path.join(self.__cellpose_folder_path, 'models', 'cyto_1'),\
                                        os.path.join(self.__cellpose_folder_path, 'models', 'nuclei'), \
                                        os.path.join(self.__parent_dir, 'models', 'epi500'), \
                                        os.path.join(self.__parent_dir, 'models', 'epi2500')]
        else:
            self.__cellpose_folder_path = None
            self.pretrained_models_paths = None
        
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

        self.flow_threshold = None
        self.cellprob_threshold = None
        self.spots_df = None
        self.tracks_df = None
        self.edges_df = None

        self.pretrained_models = ['CYTO', 'CYTO2', 'NUCLEI', 'EPI500', 'EPI2500']

        self.cellpose_default_values = {
            'TARGET_CHANNEL' : 0,
            'OPTIONAL_CHANNEL_2': 0,
            'CELLPOSE_PYTHON_FILEPATH': self.cellpose_python_filepath,
            'CELLPOSE_MODEL': self.use_model,
            'FLOW_THRESHOLD': 0.4,
            'CELLPROB_THRESHOLD': 0.0,
            'CELL_DIAMETER': 0.0,
            'USE_GPU': False,
            'SIMPLIFY_CONTOURS': True
            }
        self.trackmate_default_values = {'LINKING_MAX_DISTANCE': 15.0,
                                         'GAP_CLOSING_MAX_DISTANCE': 15.0,
                                         'MAX_FRAME_GAP': 2,
                                         'ALLOW_TRACK_SPLITTING': False,
                                         'ALLOW_TRACK_MERGING': False,
             }

        


    #### Private methods ####

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
                print("\nMerged tif files into one file: ", self.img_path)
            else:
                self.img_path = im_list[0]
                print("\nUsing tif file: ", self.img_path)
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
                                    "CELLPOSE_PYTHON_FILEPATH": self.cellpose_python_filepath,
                                    "OUTPUT_PATH": self.output_folder, "SHOW_SEGMENTATION": self.show_segmentation})

        # Specify Cellpose and TrackMate parameters not set by user
        if self.cellpose_dict == {}:
            self.cellpose_dict = self.cellpose_default_values
            del self.cellpose_dict['FLOW_THRESHOLD']
            del self.cellpose_dict['CELLPROB_THRESHOLD']
        else:
            for key in set(self.cellpose_default_values) - set(['FLOW_THRESHOLD', 'CELLPROB_THRESHOLD']):
                if key not in self.cellpose_dict.keys():
                    self.cellpose_dict[key] = self.cellpose_default_values[key]
        if self.trackmate_dict == {}:
            self.trackmate_dict = self.trackmate_default_values
        else:
            for key in self.trackmate_default_values.keys():
                if key not in self.trackmate_dict.keys():
                    self.trackmate_dict[key] = self.trackmate_default_values[key]
        self.cellpose_dict['CELLPOSE_MODEL_FILEPATH'] = self.custom_model_path

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
            try: 
                command = executable + ' --ij2 --console --run "' + jython_path + '"'
                os.system(command)
            except: 
                KeyboardInterrupt
                sys.exit(0)
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
    

    #### Public methods ####

    def run_segmentation_tracking(self):
        """
        Run cellpose segmentation and trackmate tracking.
        """
        # Ensure that imageJ filepath is valid
        if not os.path.isfile(self.imagej_filepath):
            raise OSError("Invalid imageJ filepath!")
        # Ensure that cellpose python filepath is valid
        if not os.path.isfile(self.cellpose_python_filepath):
            raise OSError("Invalid cellpose python filepath!")
        # Ensure that image folder/file path is valid
        if not os.path.isdir(self.img_folder) and not os.path.isfile(self.img_folder):
            raise OSError("Invalid image folder path!")
        
        ## Set relevant attribute values
        # Extract or set flow and cell probability threshold
        try:
            self.flow_threshold = self.cellpose_dict['FLOW_THRESHOLD']
            del self.cellpose_dict['FLOW_THRESHOLD']
        except:
            if self.use_model == 'EPI2500':
                self.flow_threshold = 0.6
            else:
                self.flow_threshold = 0.4
        try:
            self.cellprob_threshold = self.cellpose_dict['CELLPROB_THRESHOLD']
            del self.cellpose_dict['CELLPROB_THRESHOLD']
        except:
            if self.use_model == 'EPI2500':
                self.cellprob_threshold = -1
            elif self.use_model == 'EPI500':
                self.cellprob_threshold = 0.5
            else:
                self.cellprob_threshold = 0.0

        # Set model path
        if self.custom_model_path is not None:
            self.use_model = 'CUSTOM'

        # If custom model is not provided, find the path of the pretrained model
        if self.custom_model_path is None and self.use_model not in self.pretrained_models:
            raise ValueError("No custom model path provided, and use_model not in pretrained models: ", self.pretrained_models)
        elif self.custom_model_path is None and self.use_model in self.pretrained_models:
            idx = self.pretrained_models.index(self.use_model)
            self.custom_model_path = self.pretrained_models_paths[idx]
        

        # Prepare images
        self.__prepare_images()

        # Ensure that xml file has not already been created
        if os.path.isfile(self.img_path.strip(".tif") + ".xml"):
            print("Skipping segmentation and tracking, as xml file already exits at: ", self.img_path.strip(".tif") + ".xml")
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

            # Reinclude flow_threshold and cellprob in cellpose dict 
            self.cellpose_dict.update({'FLOW_THRESHOLD': self.flow_threshold, 'CELLPROB_THRESHOLD': self.cellprob_threshold})

            # Reset flow and cell probability threshold to default values after run
            np.savetxt(os.path.join(self.__cellpose_folder_path, 'params.txt'), \
                        np.array([0.4, 0.0]))
            
        # Set xml path to image path
        self.xml_path = self.img_path.strip(".tif") + ".xml"

        # Warn user about track merging and splitting
        if self.trackmate_dict['ALLOW_TRACK_MERGING'] or self.trackmate_dict['ALLOW_TRACK_SPLITTING']:
            print("\nWARNING: Track merging and/or splitting has been allowed. This may result in inconsistent velocities, \
                  as more spots belonging to the same track can be present at the same time\n")
        return

    def generate_csv_files(self, calculate_velocities = True, get_tracks = True, get_edges = True, save_csv_files = True, name = None):
        """
        Generate csv files from xml file.

        """
        if self.xml_path is not None:
            print("Starting to generate csv files from xml file now. This may take a while... \
                  \nAn XML file with 120.000 spots, 90.000 edges and 20.000 tracks takes about 6-7 minutes to process on a regular labtop.")
                  

            self.spots_df, self.tracks_df, self.edges_df = trackmate_xml_to_csv(self.xml_path, calculate_velocities=True,
                                                            get_track_features = get_tracks, get_edge_features = get_edges)
        else:
            raise OSError("No xml file provided!")

        if save_csv_files:
            self.save_csv(name)
        return


    def save_csv(self, name = None):
        """
        Save csv files to output folder.

        """
        if name is None:
            try:
                name = os.path.basename(self.img_path).strip(".tif")
            except:
                name = 'CellSegmentationTracker'

        if self.spots_df is not None:
            path_out_spots = os.path.join(self.output_folder, name + '_spots.csv')
            self.spots_df.to_csv(path_out_spots, sep = ',') 
            print("Saved spots csv file to: ", path_out_spots)
        if self.tracks_df is not None:
            path_out_tracks = os.path.join(self.output_folder, name + '_tracks.csv')
            self.tracks_df.to_csv(path_out_tracks, sep = ',') 
            print("Saved tracks csv file to: ", path_out_tracks)
        if self.edges_df is not None:
            path_out_edges = os.path.join(self.output_folder, name + '_edges.csv')
            self.edges_df.to_csv(path_out_edges, sep = ',')
            print("Saved edges csv file to: ", path_out_edges)
        return

    def print_settings(self):
        root = et.fromstring(open(self.xml_path).read())

        im_settings = root.find('Settings').find('ImageData')
        basic_settings = root.find('Settings').find('BasicSettings')

        im_settings_list = ['filename', 'folder', 'width', 'height', 'nslices', 'nframes', 'pixelwidth', 'pixelheight', 'voxeldepth','timeinterval']
        basic_settings_list = ['xstart', 'xend', 'ystart', 'yend', 'zstart', 'zend', 'tstart', 'tend']
        settings = [im_settings_list, basic_settings_list]

        print_list = []
        for i, obj in enumerate([im_settings, basic_settings]):
            settings_list = settings[i]      
            for  s in settings_list:
                print_list.append(s + f': {obj.get(s)}') 
        print("\nImage information: ", print_list)
        print("Cellpose settings: ", self.cellpose_dict)
        print("Trackmate settings: ", self.trackmate_dict, "\n")
        return


def main():
    
    model_directory = 'C:\\Users\\Simon Andersen\\miniconda3\\envs\\cellpose\\lib\\site-packages\\cellpose\\models\\cyto_0'
    cellpose_python_filepath = 'C:\\Users\\Simon Andersen\\miniconda3\\envs\\cellpose\\python.exe'
    imj_path = "C:\\Users\\Simon Andersen\\Fiji.app\\ImageJ-win64.exe"

    image_path = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\valeriias mdck data for simon\\24.08.22_698x648\\im0.tif"
    input_directory = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\valeriias mdck data for simon\\24.08.22_698x648"
    image_path = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\16.06.23_stretch_data_698x648\\im0.tif"
    im_p = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\valeriias mdck data for simon\\24.08.22_698x648\\merged.tif"
    show_output = True
    cellpose_dict = {'USE_GPU': True, 'CELL_DIAMETER': 0.0}

    # xml_path, outfolder, use_model = cyto2,nuclei, 
    xml_path = 'C:\\Users\\Simon Andersen\\Projects\\Projects\\CellSegmentationTracker\\resources\\20.09.xml'
    xml_path = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\valeriias mdck data for simon\\24.08.22_698x648\\merged.xml"
    output_directory = "C:\\Users\\Simon Andersen\\Projects\\Projects\\CellSegmentationTracker\\resources"
    
    path = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\NucleiCorrected.tif"
    new_path = path.strip(".tif") + "_resized.tif"
    np = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\merged.tif"
    path_stretch = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\16.06.23_stretch_data_split"
    pn = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\t_164_428 - 2023_03_03_TL_MDCK_2000cells_mm2_10FNh_10min_int_frame1_200_cropped_split"
    xml_path = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\t_164_428 - 2023_03_03_TL_MDCK_2000cells_mm2_10FNh_10min_int_frame1_200_cropped_split\\merged.xml"
    t1= time.time()
    cst = CellSegmentationTracker(imj_path, cellpose_python_filepath, pn, output_folder_path=output_directory, \
                                  show_segmentation=show_output, cellpose_dict=cellpose_dict, use_model='EPI2500',)
    cst.run_segmentation_tracking()
    t2= time.time()
    
    print("RUNTIME: ", (t2-t1)/60)

    cst.generate_csv_files()
    t3 = time.time()
    print("CSV gen time: ", (t3-t2)/60)
    

    #cst.run_segmentation_tracking()
    
    if 0:
        print(cst.flow_threshold, cst.cellprob_threshold)

        print(cst.spots_df.info()  )
        print(cst.spots_df.describe())
        print(cst.tracks_df.info())
        print(cst.edges_df.info())



        cst.save_csv()
        cst.print_settings()
        print("you made it bro")

if __name__ == '__main__':
    main()
