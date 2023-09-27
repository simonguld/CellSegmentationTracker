import os
import pandas as pd
import numpy as np
from matplotlib import rcParams

from CellSegmentationTracker_dev import CellSegmentationTracker

d = {'figure.figsize': (7,6)}
rcParams.update(d)


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
cellpose_dict = {
          'TARGET_CHANNEL' : 0,
          'OPTIONAL_CHANNEL_2': 0,
          'FLOW_THRESHOLD': 0.5, 
          'CELLPROB_THRESHOLD': -1.0,
          'CELL_DIAMETER': 0.0, # If 0.0, the diameter will be estimated by Cellpose
          'USE_GPU': False,
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



## It is useful to note that if you have already generated the csv files, 
# you can simply load them.
# To do this, you must leave the image path and xml_path empty when initialising 
# the CellSegmentationTracker object, after which you can provide the path to 
# the csv files as follows:
cst.spots_df = pd.read_csv('C:\\Users\\Simon Andersen\\Projects\\Projects\\CellSegmentationTracker\\resources\\epi2500_spots.csv')
cst.tracks_df = pd.read_csv('C:\\Users\\Simon Andersen\\Projects\\Projects\\CellSegmentationTracker\\resources\\epi2500_tracks.csv')
cst.edges_df = pd.read_csv('C:\\Users\\Simon Andersen\\Projects\\Projects\\CellSegmentationTracker\\resources\\epi2500_edges.csv')
    


Ngrid = 50
# Choose any additional spot features to include in the grid statistics
include_features = []

xmin, xmax = 300, 900
ymin, ymax = 300, 800
boundaries = [(xmin, xmax), (ymin, ymax)]

name = None
# Decide whether to save csv. It will be named after the image name if provided,
# otherwise as CellSegmentationTracker_grid_stats.csv. It will be saved in the 
# output folder, if provided, otherwise in the image folder.

grid_df = cst.calculate_grid_statistics(boundaries, Ngrid = Ngrid, include_features=include_features, \
    save_csv=True, name = name, save_mat_file=True, x_increases_left_to_right=False)


# Finally, to plot the velocity flow field or streamlines, 
# call the plot_velocity_field() method.

# Choose whether to plot the velocity flow field or streamlines
mode = 'field' # 'streamlines' or 'field'


# Choose which frames to plot. If [0,0], all frames will be plotted
frame_range = [0,1] # left inclusive, right exclusive

# Choose whether to calculate the average over the entire frame range
calculate_average = False

# Choose whether to animate the plot (doesn't work in jupyter notebook)
animate = False
# Choose the frame interval for the animation)
frame_interval = 1100
# Choose whether to show the plot
show = True


## Now that we have calculated the grid statistics, we can plot them.
# To plot scalar fields, call the visualize_grid_statistics() method.

# Choose which feature to plot. If None, number_density will be plotted
# NB: If you choose another feature than number_density, you must have
# included it in the grid statistics calculation
feature = 'cell_number'

# Set the unit of the feature
feature_unit = None

# Choose which frames to plot. If [0,0], all frames will be plotted
frame_range = [0,1] # left inclusive, right exclusive

# Choose whether to calculate the average over the entire frame range
calculate_average = False

# Choose whether to animate the plot (doesn't work in jupyter notebook)
animate = False
# Choose the frame interval for the animation)

if 0:
  cst.visualize_grid_statistics(feature = 'number_density', feature_unit = feature_unit,  frame_range = frame_range, \
      calculate_average = calculate_average, animate = animate, \
      frame_interval = frame_interval, show = False)
  cst.visualize_grid_statistics(feature = 'interp_number_density', feature_unit = feature_unit,  frame_range = frame_range, \
      calculate_average = calculate_average, animate = animate, \
      frame_interval = frame_interval, show = True)


  # All the remaining arguments are the same as for visualize_grid_statistics()

if 1:
  cst.plot_velocity_field(mode = mode, frame_range = frame_range, \
      calculate_average = calculate_average, animate = animate, \
      frame_interval = frame_interval, use_interpolated_velocities=True, show = show)