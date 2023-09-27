
import os
import pandas as pd
from matplotlib import rcParams

from cellsegmentationtracker import CellSegmentationTracker


def main():

    cellpose_folder_path = \
    'C:\\Users\\Simon Andersen\\miniconda3\\envs\\cellpose\\Lib\\site-packages\\cellpose'

    unit_conversion_dict= {
                                                    'frame_interval_in_physical_units': 600.213562,
                                                    'physical_length_unit_name': r'$\mu$m',
                                                        'physical_time_unit_name': 's',
            }


    cst = CellSegmentationTracker(cellpose_folder_path, unit_conversion_dict = unit_conversion_dict)

    cst.spots_df = pd.read_csv('C:\\Users\\Simon Andersen\\Projects\\Projects\\CellSegmentationTracker\\resources\\epi2500_spots.csv')
    cst.tracks_df = pd.read_csv('C:\\Users\\Simon Andersen\\Projects\\Projects\\CellSegmentationTracker\\resources\\epi2500_tracks.csv')
    cst.edges_df = pd.read_csv('C:\\Users\\Simon Andersen\\Projects\\Projects\\CellSegmentationTracker\\resources\\epi2500_edges.csv')
    

    Ngrid = 50
    # Choose any additional spot features to include in the grid statistics
    include_features = ['Area',]

    cst.calculate_grid_statistics(Ngrid = Ngrid, include_features=include_features, \
        save_csv=True, save_mat_file=True)


    feature = 'Area'

    # Set the unit of the feature
    feature_unit = r'$\mu$m$^2$'

    # Choose whether to calculate the average over the entire frame range
    calculate_average = False

    # Choose whether to animate the plot (doesn't work in jupyter notebook)
    animate = True
    # Choose the frame interval for the animation)
    frame_interval = 1200

    cst.visualize_grid_statistics(feature = 'interp_area', feature_unit = feature_unit, \
        calculate_average = calculate_average, animate = animate)

    mode = 'streamlines' # 'streamlines' or 'field'

    cst.plot_velocity_field(mode = 'streamlines',\
        calculate_average = calculate_average, animate = animate, \
        frame_interval = frame_interval)


if __name__ == '__main__':
    main()