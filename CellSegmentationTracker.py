## Imports:
import os, sys, glob
import tifftools, json
import platform
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
    t = tifftools.read_tiff(img_in_names[0])
    for img_in_name in img_in_names[1:]:
        t["ifds"].extend(tifftools.read_tiff(img_in_name)["ifds"])

    if os.path.isfile(img_out_name):
        os.unlink(img_out_name)
    tifftools.write_tiff(t, img_out_name)

    imgo = imread(img_out_name)
    imgs = [imread(e) for e in img_in_names]

    if len(img_in_names) in (3, 4):  # @TODO:
        imgo = imgo.transpose((2, 0, 1))

def trackmate_xml_to_csv(trackmate_xml_path, include_features_list = None, get_tracks=True):
    """Import detected peaks with TrackMate Fiji plugin.

    Parameters
    ----------
    trackmate_xml_path : str
        TrackMate XML file path.
    get_tracks : boolean
        Add tracks to label
    """

    root = et.fromstring(open(trackmate_xml_path).read())

    objects = []
    minimal_labels = {'FRAME': 'Frame',
                      'ID': 'Spot ID',}
    object_labels = {'FRAME': 'Frame',
                     'POSITION_T': 'T',
                     'POSITION_X': 'X',
                     'POSITION_Y': 'Y',
                     'POSITION_Z': 'Z',
                     'RADIUS': 'Radius',
                     'VISIBILITY': 'Visibility',
                     'MANUAL_SPOT_COLOR': 'Manual spot color',
                     'ELLIPSE_X0': 'Ellipse center x0',
                        'ELLIPSE_Y0': 'Ellipse center y0',
                     'ELLIPSE_MINOR': 'Ellipse short axis',
                        'ELLIPSE_MAJOR': 'Ellipse long axis',
                        'ELLIPSE_THETA': 'Ellipse angle',
                        'ELLIPSE_ASPECTRATIO': 'Ellipse aspect ratio',
                        'AREA': 'Area',
                        'PERIMETER': 'Perimeter',
                        'CIRCULARITY': 'Circularity',
                        'SOLIDITY': 'Solidity',
                        'SHAPE_INDEX': 'Shape index',
                     'QUALITY': 'Quality',      
                     'MEAN_INTENSITY_CH1': 'Mean intensity ch1',
                     'MEDIAN_INTENSITY_CH1': 'Median intensity ch1',
                     'MIN_INTENSITY_CH1': 'Min intensity ch1',
                     'MAX_INTENSITY_CH1': 'Max intensity ch1',
                     'TOTAL_INTENSITY_CH1': 'Sum intensity ch1',
                    'STD_INTENSITY_CH1': 'Std intensity ch1',
                    'CONTRAST_CH1': 'Contrast ch1',
                    'SNR_CH1': 'Signal/Noise ratio ch1',
                    'ID': 'Spot ID',
                    }
    
    # Include only features in include_features_list
    if include_features_list is not None:
        del_list = []
        for el in object_labels.keys():
            if object_labels[el] not in include_features_list:
                del_list.append(el)
        for el in del_list:
            del object_labels[el]
    object_labels.update(minimal_labels)

    features = root.find('Model').find('FeatureDeclarations').find('SpotFeatures')
    #features = [c.get('feature') for c in features.getchildren()] + ['ID']
    features = [c.get('feature') for c in list(features)] + ['ID']

    spots = root.find('Model').find('AllSpots')
    trajs = pd.DataFrame([])
    objects = []
    for frame in spots.findall('SpotsInFrame'):
        for spot in frame.findall('Spot'):
            single_object = []
            for label in features:
                single_object.append(spot.get(label))
            objects.append(single_object)

    trajs = pd.DataFrame(objects, columns=features)
    trajs = trajs.astype(float)

    # Apply initial filtering
    initial_filter = root.find("Settings").find("InitialSpotFilter")

    trajs = filter_spots(trajs,
                         name=initial_filter.get('feature'),
                         value=float(initial_filter.get('value')),
                         isabove=True if initial_filter.get('isabove') == 'true' else False)

    # Apply filters
    spot_filters = root.find("Settings").find("SpotFilterCollection")

    for spot_filter in spot_filters.findall('Filter'):

        trajs = filter_spots(trajs,
                             name=spot_filter.get('feature'),
                             value=float(spot_filter.get('value')),
                             isabove=True if spot_filter.get('isabove') == 'true' else False)

    trajs = trajs.loc[:, object_labels.keys()]
    trajs.columns = [object_labels[k] for k in object_labels.keys()]
    trajs['TRACK_ID'] = np.arange(trajs.shape[0])

    # Get tracks
    if get_tracks:
        filtered_track_ids = [int(track.get('TRACK_ID')) for track in root.find('Model').find('FilteredTracks').findall('TrackID')]

        label_id = 0
        trajs['TRACK_ID'] = np.nan

        tracks = root.find('Model').find('AllTracks')
        for track in tracks.findall('Track'):

            track_id = int(track.get("TRACK_ID"))
            if track_id in filtered_track_ids:

                spot_ids = [(edge.get('SPOT_SOURCE_ID'), edge.get('SPOT_TARGET_ID'), edge.get('EDGE_TIME')) for edge in track.findall('Edge')]
                spot_ids = np.array(spot_ids).astype('float')[:, :2]
                spot_ids = set(spot_ids.flatten())

                trajs.loc[trajs["Spot ID"].isin(spot_ids), "TRACK_ID"] = label_id
                label_id += 1

        # Label remaining columns
        single_track = trajs.loc[trajs["TRACK_ID"].isnull()]
        trajs.loc[trajs["TRACK_ID"].isnull(), "TRACK_ID"] = label_id + np.arange(0, len(single_track))

    return trajs


def filter_spots(spots, name, value, isabove):
    if isabove:
        spots = spots[spots[name] > value]
    else:
        spots = spots[spots[name] < value]

    return spots

def print_all_possible_spot_features():

    object_labels = {'FRAME': 'Frame',
                     'POSITION_T': 'T',
                     'POSITION_X': 'X',
                     'POSITION_Y': 'Y',
                     'POSITION_Z': 'Z',
                     'RADIUS': 'Radius',
                     'VISIBILITY': 'Visibility',
                     'MANUAL_SPOT_COLOR': 'Manual spot color',
                     'ELLIPSE_X0': 'Ellipse center x0',
                        'ELLIPSE_Y0': 'Ellipse center y0',
                     'ELLIPSE_MINOR': 'Ellipse short axis',
                        'ELLIPSE_MAJOR': 'Ellipse long axis',
                        'ELLIPSE_THETA': 'Ellipse angle',
                        'ELLIPSE_ASPECTRATIO': 'Ellipse aspect ratio',
                        'AREA': 'Area',
                        'PERIMETER': 'Perimeter',
                        'CIRCULARITY': 'Circularity',
                        'SOLIDITY': 'Solidity',
                        'SHAPE_INDEX': 'Shape index',
                     'QUALITY': 'Quality',      
                     'MEAN_INTENSITY_CH1': 'Mean intensity ch1',
                     'MEDIAN_INTENSITY_CH1': 'Median intensity ch1',
                     'MIN_INTENSITY_CH1': 'Min intensity ch1',
                     'MAX_INTENSITY_CH1': 'Max intensity ch1',
                     'TOTAL_INTENSITY_CH1': 'Sum intensity ch1',
                    'STD_INTENSITY_CH1': 'Std intensity ch1',
                    'CONTRAST_CH1': 'Contrast ch1',
                    'SNR_CH1': 'Signal/Noise ratio ch1',
                    'ID': 'Spot ID',
                    }
    
    print(object_labels.values())




## necessary args:

# if working pipeline:
# path to tif images, recropping argument, trackmate arguments, cellpose arguments, choice of model, save_img = True, 

# if not working pipeline
# path to (possibly) xml and definityely csv spot track edge files

##TODO:

##ESSENTIAL (DO FIRST)
# make os.system work or find another way to run jython script
# extend model path to be relative you know
# fig out a good way to name trackmate xmls
# extend args so as to be able to choose model(place it in models folder, such that ..\models\model kan be used to access model from this directory)

# evt omskriv model_path til use_model, hvis du får relativ ting path til at # udvid xml-loader

## NONESSENTIAL (IF TIME)
# resize img


class CellSegmentationTracker:

    def __init__(self, imagej_filepath, cellpose_python_filepath, image_folder = None, xml_path = None, output_folder = None,
                  use_model = 'CYTO', custom_model_path = None, show_segmentation = False, cellpose_dict = dict(), trackmate_dict = dict(),):

        self.img_folder = image_folder
        self.img_path = None
        self.xml_path = xml_path
        self.imagej_filepath = imagej_filepath
        self.fiji_folder_path = os.path.abspath(os.path.join(self.imagej_filepath, os.pardir))
        self.cellpose_python_filepath = cellpose_python_filepath
        self.custom_model_path = custom_model_path
        self.output_folder = output_folder
        self.current_path = os.getcwd()

        self.show_segmentation = show_segmentation
        self.cellpose_dict = cellpose_dict
        self.trackmate_dict = trackmate_dict
        self.use_model = use_model
        self.jython_dict = {}

        self.cellpose_default_values = {
            'TARGET_CHANNEL' : 0,
            'OPTIONAL_CHANNEL_2': 0,
            'CELLPOSE_PYTHON_FILEPATH': self.cellpose_python_filepath,
            'CELLPOSE_MODEL': 'CYTO',
            'CELLPOSE_MODEL_FILEPATH': self.cellpose_python_filepath,
            'CELL_DIAMETER': 30.0,
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
        
        # If no xml file is provided, cell segmentation, tracking and generation of an xml file is performed
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


            if 0:
                if 'TARGET_CHANNEL' not in self.cellpose_dict.keys():
                    self.cellpose_dict['TARGET_CHANNEL'] = 0
                if 'OPTIONAL_CHANNEL_2' not in self.cellpose_dict.keys():
                    self.cellpose_dict['OPTIONAL_CHANNEL_2'] = 0
                if 'CELLPOSE_MODEL' not in self.cellpose_dict.keys():
                    self.cellpose_dict['CELLPOSE_MODEL'] = use_model
                # When cell diameter is set to 0, Cellpose will estimate the cell radius
                if 'CELL_DIAMETER' not in self.cellpose_dict.keys():
                    self.cellpose_dict['CELL_DIAMETER'] = 0 
                if 'USE_GPU' not in self.cellpose_dict.keys():
                    self.cellpose_dict['USE_GPU'] = False
                if 'SIMPLIFY_CONTOURS' not in self.cellpose_dict.keys():
                    self.cellpose_dict['SIMPLIFY_CONTOURS'] = True

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

            self.jython_dict.update({'CELLPOSE_DICT': self.cellpose_dict})
            self.jython_dict.update({'TRACKMATE_DICT': self.trackmate_dict})


            # Save dictionary to json file
            with open(os.path.join(self.current_path,"jython_dict.json"), 'w') as fp:
                json.dump(self.jython_dict, fp)


            # Run jython script
            os.chdir(self.fiji_folder_path)
            executable = list(os.path.split(self.imagej_filepath))[-1]
            jython_path = os.path.join(self.current_path, "jython_cellpose.py")

            if self.show_segmentation: 
                command = executable + ' --ij2 --run "' + jython_path + '"'
            else:
                command = executable + ' --ij2 --headless --run "' + jython_path + '"'

            os.system(command)

            # Set xml path to image path
            self.xml_path = self.img_path.strip(".tif") + ".xml"

 
            print("XMLfile saved to: ", self.xml_path)
        # Generate csv file from xml file
        self.csv = trackmate_xml_to_csv(self.xml_path, get_tracks = True)

    def save_csv(self):
        if self.output_folder is not None:
            self.csv.to_csv(os.path.join(self.output_folder, 'TrackMate.csv'), sep = ',') 
        else:
            self.csv.to_csv(os.path.join(self.img_folder, 'TrackMate.csv'), sep = ',')
        return


        


# tests:
# Get it working without gui
# include merging tiffs
# get it working using models in ./models[also, fix model_path for user who simply can choose a must]
# get it working with all args() DOES DICT ARG WORK?
# fix naming
#fix xml output path


def main():
    image_path = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\valeriias mdck data for simon\\24.08.22_698x648\\im0.tif"
    input_directory = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\valeriias mdck data for simon\\24.08.22_698x648"
    output_directory = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\valeriias mdck data for simon\\24.08.22_698x648"
    model_directory = 'C:\\Users\\Simon Andersen\\miniconda3\\envs\\cellpose\\lib\\site-packages\\cellpose\\models\\cyto_0'
    cellpose_python_filepath = 'C:\\Users\\Simon Andersen\\miniconda3\\envs\\cellpose\\python.exe'
    imj_path = "C:\\Users\\Simon Andersen\\Fiji.app\\ImageJ-win64.exe"

    show_output = True
    cellpose_dict = {'USE_GPU': True}
    trackmate_dict = {'MAX_FRAME_GAP': 3}
  
    cst = CellSegmentationTracker(imj_path, cellpose_python_filepath, image_path, \
                                  show_segmentation=show_output, cellpose_dict=cellpose_dict, trackmate_dict=trackmate_dict)
    print(cst.cellpose_dict)
    print(cst.trackmate_dict)
    cst.save_csv()
    print("you made it bro")

if __name__ == '__main__':
    main()