import os
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata


if 0:
    def __interpolate_grid_velocities(self, interpolation_method = 'cubic'):
        """
        Interpolates the missing velocities in the grid dataframe using the grid centers as interpolation points.

        Parameters:
        -----------
        interpolation_method : (string, default = 'cubic') - interpolation method to use. Possible values are 'linear', 'nearest' and 'cubic'.
        """
 
        # Etract grid csv as array
        grid_data = self.grid_df.values

        # The columns of the grid dataframe for reference
        grid_columns = ['Frame', 'T', 'Ngrid','x_center', 'y_center', 'cell_number', 'number_density','mean_velocity_X','mean_velocity_Y']

        # Find grid centers
        x_centers = np.unique(grid_data[:, grid_columns.index('x_center')])
        y_centers = np.unique(grid_data[:, grid_columns.index('y_center')])

        # Create grid stack, ie. array of all grid centers with shape (Ngrid**2, 2)
        grid_stack = np.array(np.meshgrid(x_centers, y_centers)).T.reshape(-1,2)

        # Initialize array for interpolated velocities
        v_stack_full = np.zeros([0,2])

        for frame in np.arange(self.__Nframes):
            # get mask for current frame
            frame_mask = (grid_data[:, grid_columns.index('Frame')] == frame)

            # Extract velocities for current frame
            vx = grid_data[frame_mask, grid_columns.index('mean_velocity_X')]
            vy = grid_data[frame_mask, grid_columns.index('mean_velocity_Y')]

            # Stack velocities, into array with shape (Ngrid**2, 2)
            v_stack = np.vstack([vx, vy]).T
            
            # Find nan values
            nan_mask = np.isnan(v_stack[:,0])

            # Interpolate nan values for each component
            for i in [0,1]:
                v_stack[nan_mask, i] = griddata(grid_stack[~nan_mask], v_stack[~nan_mask,i], \
                                                grid_stack[nan_mask], method = interpolation_method)
                
                # Find remaining nan values
                remaining_nans_mask = np.isnan(v_stack[:,i])

                # Estimate remaning nan points at the boundary by nearest interpolation
                if remaining_nans_mask.any():
                    v_stack[remaining_nans_mask, i] = griddata(grid_stack[~remaining_nans_mask], v_stack[~remaining_nans_mask,i], \
                                                            grid_stack[remaining_nans_mask], method='nearest')

            v_stack_full = np.vstack([v_stack_full, v_stack])

        self.grid_df['interp_velocity_X'] = v_stack_full[:,0]
        self.grid_df['interp_velocity_Y'] = v_stack_full[:,1]
        return


def main():
    path = 'C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\gHZIcyBQnp\\MDCK_well2_grid_stats_mesh_size_20.csv'
    dir_name = os.path.dirname(path)

    grid_csv = pd.read_csv(path, sep=',')
    Nframes = grid_csv['Frame'].max() + 1
    
    x_centers = np.unique(grid_csv['x_center'].values)
    y_centers = np.unique(grid_csv['y_center'].values)

    grid_stack = np.array(np.meshgrid(x_centers, y_centers)).T.reshape(-1,2)

    v_stack_full = np.zeros([0,2])

    print("Npoints :",grid_stack.shape)

    for frame in np.arange(Nframes):
        vx = grid_csv[grid_csv['Frame'] == frame]['mean_velocity_X'].values
        vy = grid_csv[grid_csv['Frame'] == frame]['mean_velocity_Y'].values

        v_stack = np.vstack([vx, vy]).T
        nan_idx = [1, 5, 9, 11, 19, 21,]
        v_stack[nan_idx] = [np.nan, np.nan]

        nan_mask = np.isnan(v_stack[:,0])

        print("Fraction of nans: ", nan_mask.sum() / len(nan_mask))

        print(grid_stack.shape, v_stack.shape)

        for i in [0,1]:
            v_stack[nan_mask, i] = griddata(grid_stack[~nan_mask], v_stack[~nan_mask,i], grid_stack[nan_mask], method='cubic')
            remaining_nans_mask = np.isnan(v_stack[:,i])

            # Estimate remaning nan points at the boundary by nearest interpolation
            if remaining_nans_mask.any():
                v_stack[remaining_nans_mask, i] = griddata(grid_stack[~remaining_nans_mask], v_stack[~remaining_nans_mask,i], \
                                                        grid_stack[remaining_nans_mask], method='nearest')

        v_stack_full = np.vstack([v_stack_full, v_stack])
   
    # save new grid csv
    grid_csv['interp_velocity_X'] = v_stack_full[:,0]
    grid_csv['interp_velocity_Y'] = v_stack_full[:,1]
    grid_csv.to_csv(dir_name + '\\grid_statistics_interpolated_mesh20.csv', index=False)



if __name__ == '__main__':
    main()


if 0:
    
    xmin, xmax = x_centers.min(), x_centers.max()
    ymin, ymax = y_centers.min(), y_centers.max()

    Nxgrid = len(x_centers)
    Nygrid = len(y_centers)
    mesh_size = 1
    
    x_new, y_new = np.mgrid[xmin:xmax:mesh_size, ymin:ymax:mesh_size]
    grid_new = np.array(np.meshgrid(x_new, y_new)).T.reshape(-1,2)

    print(grid_new)
    print(x_new)
    print(y_new)

    

    for idx in nan_idx:
        v_stack[idx] = [np.nan, np.nan]

    v_res = v_stack.reshape(Nxgrid, Nygrid, 2)

    print(v_res[:,:,0])


    x_centers = np.linspace(0, 5, 6)
    y_centers = 1* np.linspace(0, 5, 6)

    grid_stack = np.array(np.meshgrid(x_centers, y_centers)).T.reshape(-1,2)
  
    vx = 2 * np.linspace(0, 5, 6)
    vy = -1 * np.linspace(0, 5, 6)
    v_stack = np.array(np.meshgrid(vx, vy)).T.reshape(-1,2)

    nan_idx = [1, 5, 9, 11, 19, 21,]
    v_stack[nan_idx] = [np.nan, np.nan]
    
    nan_mask = np.isnan(v_stack[:,0])


    for i in [0,1]:
        v_stack[nan_mask, i] = griddata(grid_stack[~nan_mask], v_stack[~nan_mask,i], grid_stack[nan_mask], method='cubic')
        remaining_nans_mask = np.isnan(v_stack[:,i])
        
        print(v_stack[:,i])
        # Estimate remaning nan points at the boundary by nearest interpolation
        if remaining_nans_mask.any():
            v_stack[remaining_nans_mask, i] = griddata(grid_stack[~remaining_nans_mask], v_stack[~remaining_nans_mask,i], \
                                                       grid_stack[remaining_nans_mask], method='nearest')

        print(v_stack[:,i])