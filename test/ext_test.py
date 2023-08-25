import os
import time
import pandas as pd
from matplotlib import rcParams

from cellsegmentationtracker import CellSegmentationTracker


d = {'figure.figsize': (6,6)}
rcParams.update(d)



## Set paths to executables:
#  If you want to do segmentation and tracking, you need to provide paths to ImageJ and Cellpose. 
# If you already have an xml file and just want to use CellSegmentationTracker to generate csv files, 
# calculate velocities or analyse the results, you don't have to provide these
cellpose_python_filepath = 'C:\\Users\\Simon Andersen\\miniconda3\\envs\\cellpose\\python.exe'
imj_path = "C:\\Users\\Simon Andersen\\Fiji.app\\ImageJ-win64.exe"

# Set path to .tif image file (or folder with images). If you already have an xml file, 
# you don't have to provide it
image_path = "C:\\Users\\Simon Andersen\\Projects\\Projects\\CellSegmentationTracker\\resources\\epi2500.tif"
image_path = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\16.06.23_stretch_data_split_orig\\im9.tif"
# If you have already done segmentation and tracking, you can simply provide the path to the xml file.
xml_path = None 

# Set path to output folder. If None, results will be outputted in the same folder as the image.
output_folder_path = None

# Set whether to use a pretrained model or not. If not, you need to provide the path to a custom model
use_model = 'CYTO'
custom_model_path = None

# Set whether to open Fiji and show the segmentation and tracking results. 
# If you choose to show the results, you must close the window before the program can continue
show_segmentation = False
# Set cellpose and trackmate settings. If you don't provide any, the default settings will be used
cellpose_dict = {
          'TARGET_CHANNEL' : 0,
          'OPTIONAL_CHANNEL_2': 0,
          'FLOW_THRESHOLD': 0.5, 
          'CELLPROB_THRESHOLD': -1.0,
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
            }


# Now having set all parameters, we are ready to initialise the CellSegmentationTracker object
cst = CellSegmentationTracker(imagej_filepath = imj_path, \
    cellpose_python_filepath = cellpose_python_filepath, image_folder_path = image_path, \
      xml_path = xml_path, output_folder_path = output_folder_path,\
      use_model = use_model, custom_model_path = custom_model_path,\
      show_segmentation = show_segmentation, cellpose_dict = cellpose_dict, trackmate_dict = trackmate_dict,)


# If you are using the epi2500 sample images (without GPU), it will take 5-6 minutes to run
# If you are showing the results, you must close the window before the program can continue
t1 = time.time()
cst.run_segmentation_tracking()
t2 = time.time()
print('Segmentation and tracking took {} seconds'.format(t2-t1))
