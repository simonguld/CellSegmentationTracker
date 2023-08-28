import os
import time
import pandas as pd
from matplotlib import rcParams

from cellsegmentationtracker import CellSegmentationTracker


d = {'figure.figsize': (6,6)}
rcParams.update(d)



## Set paths to executables:
# Set the path to the cellpose folder. It can be found in the virtual environment created for
#running cellpose. 
# On windows, it typically found in
#  ./path_to_anaconda_folder/envs/cellpose/Lib/site-packages/cellpose
# On MAC, it is typically found in 
# ./path_to_anaconda_folder/envs/cellpose/lib/python3.[insert version here]/site-packages/cellpose
cellpose_folder_path = \
  'C:\\Users\\Simon Andersen\\miniconda3\\envs\\cellpose\\Lib\\site-packages\\cellpose'


#  If you want to do segmentation and tracking, you need to provide paths to 
# ImageJ and Cellpose. 
# If you already have an xml file and just want to use CellSegmentationTracker to generate
# csv files, calculate velocities or analyse the results, you don't have to provide these
cellpose_python_filepath = \
  'C:\\Users\\Simon Andersen\\miniconda3\\envs\\cellpose\\python.exe'
imj_path = "C:\\Users\\Simon Andersen\\Fiji.app\\ImageJ-win64.exe"

# Set path to .tif image file (or folder with images). If you already have an xml file, 
# you don't have to provide it
image_path = "C:\\Users\\Simon Andersen\\Projects\\Projects\\CellSegmentationTracker\\resources\\epi2500.tif"

# If you have already done segmentation and tracking, you can simply provide the path
#  to the xml file.
xml_path = None 

# Set path to output folder. If None, results will be outputted in the same folder
#  as the image.
output_folder_path = None

# Set whether to use a pretrained model or not. If not, you need to provide the path
#  to a custom model
use_model = 'EPI2500'
custom_model_path = None

# Set whether to open Fiji and show the segmentation and tracking results. 
# If you choose to show the results, you must close the window before the program
#  can continue
show_segmentation = True

# Set cellpose and trackmate settings. If you don't provide any, the default settings
#  will be used
cellpose_dict = {'USE_GPU': True, 'CELL_DIAMETER': 0.0}


# Finally, you can set the names of the physical units used in the image, as well
#  as the frame interval in physical units. If you don't provide any, the default
#  values will be used (pixels and frames)
# For the epi2500 images, the frame interval is 5 minutes.
# NB: The length scales will be automatically converted to physical units
# by TrackMate, if the physical units are set correctly in the .tif file.
unit_conversion_dict= {
                                                  'frame_interval_in_physical_units': 5,
                                                  'physical_length_unit_name': r'$\mu$m',
                                                    'physical_time_unit_name': 'min',
        }


# Now having set all parameters, we are ready to initialise the CellSegmentationTracker
#  object:
cst = CellSegmentationTracker(cellpose_folder_path, imagej_filepath = imj_path, \
    cellpose_python_filepath = cellpose_python_filepath, image_folder_path = image_path, \
      xml_path = xml_path, output_folder_path = output_folder_path,\
      use_model = use_model, custom_model_path = custom_model_path,\
      show_segmentation = show_segmentation, cellpose_dict = cellpose_dict, \
      unit_conversion_dict = unit_conversion_dict)



# If you are using the epi2500 sample images (without GPU), it will take 5-6 minutes to run
# If you are showing the results, you must close the window before the program can continue
t1 = time.time()
cst.run_segmentation_tracking()
t2 = time.time()
print('Segmentation and tracking took {} seconds'.format(t2-t1))

