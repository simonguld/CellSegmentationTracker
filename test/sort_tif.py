import os
import sys   

import pandas as pd
from matplotlib import rcParams

from CellSegmentationTracker_dev import CellSegmentationTracker


d = {'figure.figsize': (6,6)}
rcParams.update(d)



def get_imlist(path, format = '.jpg'):

    """
    returns a list of filenames for all png images in a directory
    """
    if os.path:
        return [os.path.join(path,f) for f in os.listdir(path) if f.endswith(format)]
    else:
        raise OSError("No such directory: " + path)
    



img_folder = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\data\\16.06.23_stretch_data"


im_list = get_imlist(img_folder, format = '.tif')



cellpose_folder_path = \
  'C:\\Users\\Simon Andersen\\miniconda3\\envs\\cellpose\\Lib\\site-packages\\cellpose'
imj_path = "C:\\Users\\Simon Andersen\\Fiji.app\\ImageJ-win64.exe"

cellpose_python_filepath = \
  'C:\\Users\\Simon Andersen\\miniconda3\\envs\\cellpose\\python.exe'

image_path = "C:\\Users\\Simon Andersen\\Projects\\Projects\\CellSegmentationTracker\\resources\\epi2500.tif"


print(im_list)

# Set whether to use a pretrained model or not. If not, you need to provide the path
#  to a custom model
use_model = 'EPI6000'
custom_model_path = None

# Set whether to open Fiji and show the segmentation and tracking results. 
# If you choose to show the results, you must close the window before the program
#  can continue
show_segmentation = False
output_folder = img_folder
# Set cellpose and trackmate settings. If you don't provide any, the default settings
#  will be used
cellpose_dict = {
          'TARGET_CHANNEL' : 0,
          'OPTIONAL_CHANNEL_2': 0,
          'CELL_DIAMETER': 0.0, # If 0.0, the diameter will be estimated by Cellpose
          'USE_GPU': True,
          'SIMPLIFY_CONTOURS': True
          }

# Beware that if you set ALLOW_TRACK_SPLITTING and/or ALLOW_TRACK_MERGING to True, 
# the calculated velocities might be incorrect, as several cells will be merged into
#  one track and can be present at the same time
trackmate_dict = {'LINKING_MAX_DISTANCE': 15.0,
                                        'GAP_CLOSING_MAX_DISTANCE': 30.0,
                                        'MAX_FRAME_GAP': 3,
                                        'ALLOW_TRACK_SPLITTING': False, 
                                        'ALLOW_TRACK_MERGING': False,
                                        'ALLOW_GAP_CLOSING': True
            }

# Now having set all parameters, we are ready to initialise the CellSegmentationTracker
#  object:
if 0:

  cst = CellSegmentationTracker(cellpose_folder_path, imagej_filepath = imj_path, \
      cellpose_python_filepath = cellpose_python_filepath, image_folder_path=image_path, output_folder_path=output_folder, \
        use_model = use_model, custom_model_path = custom_model_path,\
        show_segmentation = show_segmentation, cellpose_dict = cellpose_dict, \
          trackmate_dict = trackmate_dict,)

  cst.run_segmentation_tracking()
if 1:
  for im_path in [im_list[2]]:
    cst = CellSegmentationTracker(cellpose_folder_path, imagej_filepath = imj_path, \
    cellpose_python_filepath = cellpose_python_filepath, image_folder_path = im_path, \
      use_model = use_model, custom_model_path = custom_model_path,\
      show_segmentation = show_segmentation, cellpose_dict = cellpose_dict, \
        trackmate_dict = trackmate_dict,)
    
    cst.run_segmentation_tracking()

    df_spots, df_tracks, df_edges = \
    cst.generate_csv_files()

