import os
import pandas as pd
import numpy as np
from matplotlib import rcParams

from CellSegmentationTracker_dev import CellSegmentationTracker

d = {'figure.figsize': (7,6)}
rcParams.update(d)

# Set whether to run segmentation tracking, whether to generate csv files or whether to simply load csv files
run_segmentation_tracking = False
generate_csv_files = False
load_csv_files = False

#  If you want to do segmentation and tracking, you need to provide paths to the
# cellpose folder , the ImageJ program and cellpose python executable.

# Set the path to the cellpose folder. It can be found in the virtual environment created for
#running cellpose. 
# On windows, it typically found in
#  ./path_to_anaconda_folder/envs/cellpose/Lib/site-packages/cellpose
# On MAC, it is typically found in 
# ./path_to_anaconda_folder/envs/cellpose/lib/python3.[insert version here]/site-packages/cellpose
cellpose_folder_path = \
  'C:\\Users\\Simon Andersen\\miniconda3\\envs\\cellpose\\Lib\\site-packages\\cellpose'

# The ImageJ program is found in the Fiji.app folder, and has the extension .exe on Windows
#  and -macosx on MAC.
imj_path = "C:\\Users\\Simon Andersen\\Fiji.app\\ImageJ-win64.exe"

# the python executable is found in the virtual environment created for running cellpose.
cellpose_python_filepath = \
  'C:\\Users\\Simon Andersen\\miniconda3\\envs\\cellpose\\python.exe'

# Set path to .tif image file (or folder with images). If you already have an xml file, 
# you don't have to provide it
image_path = "C:\\Users\\Simon Andersen\\Projects\\Projects\\CellSegmentationTracker\\resources\\epi2500.tif"
image_path = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\data\\valeriias mdck data for simon\\24.08.22_698x648\\merged1.tif"

# If you have already done segmentation and tracking, you can simply provide the path
#  to the xml file.
xml_path = None

# Set path to output folder. If None, results will be outputted in the same folder
#  as the image.
output_folder_path = None

# Set whether to use a pretrained model or not. If not, you need to provide the path
#  to a custom model
use_model = 'EPI500'
custom_model_path = None

# Set whether to open Fiji and show the segmentation and tracking results. 
# If you choose to show the results, you must close the window before the program
#  can continue. If true, the overlays will be saved as .avi file in the image folder
show_segmentation = True

# Set cellpose and trackmate settings. If you don't provide any, the default settings
#  will be used
cellpose_dict = {
          'TARGET_CHANNEL' : 0,
         # 'OPTIONAL_CHANNEL_2': 0,
          'FLOW_THRESHOLD': 0.5, 
        #  'CELLPROB_THRESHOLD': -1.0,
          'CELL_DIAMETER': 0.0, # If 0.0, the diameter will be estimated by Cellpose
          'USE_GPU': True,
          'SIMPLIFY_CONTOURS': True
          }

# Beware that if you set ALLOW_TRACK_SPLITTING and/or ALLOW_TRACK_MERGING to True, 
# the calculated velocities might be incorrect, as several cells will be merged into
#  one track and can be present at the same time
trackmate_dict = {'LINKING_MAX_DISTANCE': 15.0,
                                        'GAP_CLOSING_MAX_DISTANCE': 15.0,
                                        'MAX_FRAME_GAP': 2,
                                        'ALLOW_TRACK_SPLITTING': False, 
                                        'ALLOW_TRACK_MERGING': False,
                                        'ALLOW_GAP_CLOSING': True,
            }

# Finally, you can set the names of the physical units used in the image, as well
#  as the frame interval in physical units. If you don't provide any, the default
#  values will be used (pixels and frames)
# For the epi2500 images, the frame interval is 600 seconds, but since this information
# is already stored in the .tif file, we don't need to provide it here.
# NB: The length and time scales will be automatically converted to physical units
# by TrackMate, if the physical units are set correctly in the .tif file.
unit_conversion_dict= {
                                                  'frame_interval_in_physical_units': 600,
                                                  'physical_length_unit_name': r'$\mu$m',
                                                    'physical_time_unit_name': 's',
        }


# Now having set all parameters, we are ready to initialise the CellSegmentationTracker
#  object:
cst = CellSegmentationTracker(cellpose_folder_path, imagej_filepath = imj_path, \
    cellpose_python_filepath = cellpose_python_filepath, image_folder_path = image_path, \
      xml_path = xml_path, output_folder_path = output_folder_path,\
      use_model = use_model, custom_model_path = custom_model_path,\
      show_segmentation = show_segmentation, cellpose_dict = cellpose_dict, \
        trackmate_dict = trackmate_dict, unit_conversion_dict = unit_conversion_dict)


images_indices = [0, 7]

cst.do_sample_segmentation(image_indices_to_segment=images_indices, \
                           cellpose_model_to_estimate_diameter='CYTO2', plot_overlays=True, plot_as_you_go=False, save_overlays=True)