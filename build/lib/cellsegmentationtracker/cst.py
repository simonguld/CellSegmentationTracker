import pandas as pd
from CellSegmentationTracker_dev import CellSegmentationTracker
import scipy.io
import numpy as np

def main():

    # Load matlab file
    path32 = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\gHZIcyBQnp\\MDCK_well2_piv_settings_32_16_16.mat"
    path64 = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\gHZIcyBQnp\\MDCK_well2_piv_settings_64_32_32_16.mat"
    spots_path = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\gHZIcyBQnp\\MDCK_well2_crop_xyCorrected_spots.csv"

    
    
    mat32 = scipy.io.loadmat(path32)
    math64 = scipy.io.loadmat(path64)

    x64 = math64['X'][0][0]
    y64 = math64['Y'][0][0]

    x32 = mat32['X'][0][0]
    y32 = mat32['Y'][0][0]

    xmin, xmax = x32.min(), x32.max()
    ymin, ymax = y32.min(), y32.max()
    mesh_size_x = x32[1][0] - x32[0][0]
    mesh_size_y = y32[0][1] - y32[0][0]
    assert(mesh_size_x == mesh_size_y)
    Lx = xmax - xmin + mesh_size_x
    Ly = ymax - ymin + mesh_size_y

    Ngrid = int(min(Lx, Ly) / mesh_size_x)

    print(min(Lx, Ly) / mesh_size_x)
    print(xmin, xmax, ymin, ymax, mesh_size_x, mesh_size_y, Ngrid)

    # Initialize CellSegmentationTracker object
    output_folder= "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\gHZIcyBQnp"
    cst = CellSegmentationTracker(output_folder_path=output_folder)

    # Set path to spots csv file
    spots_df = pd.read_csv(spots_path)

    # Crop csv
    cst.spots_df = spots_df[(spots_df['X'] >= xmin) & (spots_df['X'] <= xmax) & (spots_df['Y'] >= ymin) & (spots_df['Y'] <= ymax)]

    cst.grid_df = pd.read_csv("C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\gHZIcyBQnp\\MDCK_well2_grid_stats.csv")

  #  cst.calculate_grid_statistics(Ngrid = Ngrid, name = "MDCK_well2_grid_stats",)

    print(cst.grid_df.loc[:,['x_center', 'y_center']].describe())
   

if __name__ == '__main__':
    main()