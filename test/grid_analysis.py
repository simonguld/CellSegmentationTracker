# Author: Simon Guldager Andersen
# Date (latest update): 

### SETUP ------------------------------------------------------------------------------------

## Imports:
import os
import sys
import time
import warnings
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats, integrate, interpolate, optimize
from scipy.special import sici, factorial
from iminuit import Minuit
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.animation as ani
from matplotlib import rcParams
from cycler import cycler

sys.path.append('..\\Appstat2022\\External_Functions')
from ExternalFunctions import Chi2Regression, BinnedLH, UnbinnedLH
from ExternalFunctions import nice_string_output, add_text_to_ax    # Useful functions to print fit results on figure


## Change directory to current one
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)

## Set plotting style and print options
sns.set_theme()
sns.set_style("darkgrid")
sns.set_context("paper") #Possible are paper, notebook, talk and poster

d = {'lines.linewidth': 2, 'axes.titlesize': 18, 'axes.labelsize': 18, 'xtick.labelsize': 12, 'ytick.labelsize': 12,\
     'legend.fontsize': 15, 'font.family': 'serif', 'figure.figsize': (9,6)}
d_colors = {'axes.prop_cycle': cycler(color = ['teal', 'navy', 'coral', 'plum', 'purple', 'olivedrab',\
         'black', 'red', 'cyan', 'yellow', 'khaki','lightblue'])}
rcParams.update(d)
rcParams.update(d_colors)
np.set_printoptions(precision = 5, suppress=1e-10)


### FUNCTIONS ----------------------------------------------------------------------------------

### MAIN ---------------------------------------------------------------------------------------


def calculate_grid_statistics(dataframe, Ngrid, include_features = [], return_absolute_cell_counts = False,\
                               save_csv = True, name_path = './grid_statistics'):
    """
    Calculates the mean value of a given feature in each grid square for each frame and returns a dataframe with the results.

    Parameters:
    ----------
    dataframe : pandas dataframe generated by the function 'generate_csv_files'
    Ngrid : int - number of grid squares in the smallest dimension. The number of grid squares in the other dimension 
            is determined by the aspect ratio of the image, with the restriction that the grid squares are square.
    include_features : list of strings - list of features to include in the grid dataframe, in addition to the standard features
                       number_density, mean_velocity_X and mean_velocity_Y. The possible feature keys are the columns of the 
                       spots dataframe generated by the function 'generate_csv_files'.
    return_absolute_cell_counts : bool - if True, the number of cells in each grid square is returned instead of the density.
    save_csv : bool - if True, the grid dataframe is saved as a csv file.
    name : string - name of the csv file to be saved. It will be saved in the output_folder, if provided, otherwise in the image folder

    Returns:
    -------
    grid_df : pandas dataframe - dataframe containing the grid statistics for each frame.
    """

    cols = ['Frame', 'X', 'Y', 'T', 'Velocity_X', 'Velocity_Y']
    add_to_grid = []
    for feature in include_features:
        if feature not in cols:
            cols.append(feature)    
            add_to_grid.append("mean_" + feature.lower())

    # Extract relevant data from dataframe as array
    data = dataframe.loc[:, cols].values

    grid_columns = ['Frame', 'T', 'Ngrid','x_center', 'y_center', 'number_density','mean_velocity_X','mean_velocity_Y']
    grid_columns.extend(add_to_grid)

    # Create empty array to store grid data
    Nframes = int(data[:,0].max() + 1)
    
    xmin, xmax = data[:,1].min(), data[:,1].max()
    ymin, ymax = data[:,2].min(), data[:,2].max()
    Lx, Ly = xmax - xmin, ymax - ymin

    # Calculate number of grid squares in x and y direction
    if Ly < Lx:
        Nygrid = Ngrid
        Nxgrid = int(np.floor(Ngrid * Lx / Ly))     
        grid_len = Ly / Nygrid
        Residual_x = Lx % grid_len
        xmin += Residual_x / 2
        xmax -= Residual_x / 2
    elif Lx < Ly:
        Nxgrid = Ngrid
        Nygrid = int(np.floor(Ngrid * Ly / Lx))
        grid_len = Lx / Nxgrid
        Residual_y = Ly % grid_len
        ymin += Residual_y / 2
        ymax -= Residual_y / 2

    Nsquares = Nxgrid * Nygrid
    dA = grid_len**2

    grid_arr = np.zeros([Nframes * Nsquares, len(grid_columns)])

    # Loop over all frames
    for i, frame in enumerate(np.arange(Nframes)):
        arr_T = data[data[:, 0] == frame]
        time = arr_T[0, 3]
        Ngrid = 0
        for y in np.linspace(ymax - grid_len, ymin, Nygrid):
            mask = (arr_T[:, 2] >= y) & (arr_T[:, 2] < y + grid_len)
            arr_Y = arr_T[mask]

            for x in np.linspace(xmin, xmax - grid_len, Nxgrid):
                mask = (arr_Y[:, 1] >= x) & (arr_Y[:, 1] < x + grid_len)
          
                if mask.sum() == 0:
                    grid_arr[frame * Nsquares + Ngrid, 0:6] = [frame, time, Ngrid, x + grid_len / 2, y + grid_len / 2, 0]
                    grid_arr[frame * Nsquares + Ngrid, 6:] = np.nan * np.zeros(len(grid_columns) - 6)
                else:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", category=RuntimeWarning)
                        vx = np.nanmean(arr_Y[mask, 4])
                        vy = np.nanmean(arr_Y[mask, 5])

                    # If true, return the number of cells in each grid as opposed to the density
                    if return_absolute_cell_counts:
                        density = mask.sum() 
                    else:
                        denisty = mask.sum() / dA      

                    grid_arr[frame * Nsquares + Ngrid, 0:8] = [frame, time, Ngrid, x + grid_len / 2, \
                                                               y + grid_len / 2, density, vx, vy]

                    # Add additional features to grid (if any)
                    for j, feature in enumerate(add_to_grid):
                        grid_arr[frame * Nsquares + Ngrid, 8 + j] = np.nanmean(arr_Y[mask, 6 + j])

                Ngrid += 1

    grid_df = pd.DataFrame(grid_arr, columns = grid_columns)

    if save_csv:
        grid_df.to_csv(name_path, sep = ',')

    return grid_df

def visualize_grid_statistics(grid_dataframe, feature = 'number_density', frame_range = [0,0], calculate_average = False, \
                             animate = True, frame_interval = 800, show = True):
    
    xmin, xmax = grid_dataframe['x_center'].min() , grid_dataframe['x_center'].max()
    ymin, ymax = grid_dataframe['y_center'].min(), grid_dataframe['y_center'].max()
    vmin, vmax = grid_dataframe[feature].min(), grid_dataframe[feature].max()
    Nsquares = int(np.round(grid_dataframe['Ngrid'].max()) + 1)
    Nframes = int(np.round(grid_dataframe['Frame'].max() + 1)) if frame_range == [0,0] else frame_range[1] - frame_range[0]
    frame_range[1] = frame_range[1] if frame_range[1] != 0 else Nframes

    Lx, Ly = xmax - xmin, ymax - ymin

    # Calculate number of grid squares in x and y direction
    if Ly < Lx:
        y_arr = grid_dataframe['y_center'].diff().values
        idx = int(np.argwhere(y_arr > 0.5)[0])

        grid_len = grid_dataframe['x_center'].loc[idx + 1] - grid_dataframe['x_center'].loc[0]
        Lx += grid_len
        Ly += grid_len
        Nygrid = int(np.round(Ly / grid_len))
        Nxgrid = int(np.round(Nsquares / Nygrid) )
    elif Lx <= Ly:
        grid_len = grid_dataframe['x_center'].loc[1] - grid_dataframe['x_center'].loc[0]
        Lx += grid_len
        Ly += grid_len
        Nxgrid = int(np.round(Lx / grid_len))
        Nygrid = int(np.round(Nsquares / Nxgrid) )

    data = grid_dataframe.loc[:, ['Frame', 'T', 'Ngrid', feature]].values
    data = data[frame_range[0] * Nsquares:frame_range[1] * Nsquares,:]

    if calculate_average:
        arr_T = data[:,-1].reshape(Nframes, Nsquares)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            arr_T = np.nanmean(arr_T, axis = 0)
        title = title=f'{feature.capitalize()} heatmap averaged over frames {frame_range[0]}:{frame_range[1]}'
        heatmap_plotter(0, arr_T, Nxgrid, Nygrid, xmin, xmax, ymin, ymax, grid_len, feature, vmin, vmax,\
                 title=title)
    else:
        heatmap_plotter_res = lambda i, frame : heatmap_plotter(i, frame, Nxgrid, Nygrid, xmin, xmax, ymin, ymax, grid_len, feature, vmin, vmax)
        if animate:
            animator(data, heatmap_plotter_res, (frame_range[0],frame_range[1]), inter=frame_interval, show=show)
        else:
            for i in range(frame_range[0], frame_range[1]):
                frame = data[data[:, 0] == i, -1]
                heatmap_plotter_res(i, frame)
    if show:
        plt.show()
    return

def heatmap_plotter(i, frame, Nxgrid, Nygrid, xmin, xmax, ymin, ymax, grid_len, feature, vmin, vmax, title=None):
    """
    Plots a heatmap of the given feature for a given frame.
    """
    fig, ax = plt.subplots()
    # reshape to grid
    #arr_grid = np.flip(frame.reshape(Nxgrid, Nygrid).T, axis = 0)
    arr_grid = frame.reshape(Nygrid, Nxgrid)
    Lx, Ly = xmax - xmin, ymax - ymin

    xticklabels = np.round(np.linspace(xmin, xmax, 4), 2) 
    yticklabels = np.round(np.linspace(ymax, ymin, 4), 2)

    ax = sns.heatmap(arr_grid, cmap = 'viridis', vmin=vmin, vmax=vmax)

    title = f'{feature.capitalize()} heatmap for frame = {i}' if title is None else title
    ax.set(xlabel = 'x', ylabel = 'y', title = title)
    ax.set_xticks(ticks = np.linspace(0.5,Nxgrid-0.5,4), labels=xticklabels)
    ax.set_yticks(ticks = np.linspace(0.5,Nygrid-0.5,4), labels=yticklabels)
    return

def animator(arr, fn, rng, animate_wrapper = None, inter=200, show=True):
    """Show a frame-by-frame animation.

    Parameters:
    arr -- the data array
    fn -- the plot function (argument: fram number, frame)
    rng -- range of the frames to be plotted
    interval -- time between frames (ms)
    """

    # create the figure
    fig = plt.figure()

    if animate_wrapper is None:
        # the local animation function
        def animate_fn(i):
            # we want a fresh figure everytime
            fig.clf()
            # add subplot, aka axis
            #ax = fig.add_subplot(111)
            # load the frame
            frame = arr[arr[:, 0] == i][:, -1]
            # call the global function
            fn(i, frame)
    else:
        animate_fn = lambda i: animate_wrapper(i, arr, fig, fn)
            
    anim = ani.FuncAnimation(fig, animate_fn,
                             frames=np.arange(rng[0], rng[1]),
                             interval=inter, blit=False)
    if show==True:
      plt.show()
      return

    return anim

def plot_velocity_field(grid_dataframe, frame_range = [0,0], calculate_average = False, \
                                animate = True, frame_interval = 800, show = True):
    

    xmin, xmax = grid_dataframe['x_center'].min() , grid_dataframe['x_center'].max()
    ymin, ymax = grid_dataframe['y_center'].min(), grid_dataframe['y_center'].max()

    Nsquares = int(np.round(grid_dataframe['Ngrid'].max()) + 1)
    Nframes = int(np.round(grid_dataframe['Frame'].max() + 1)) if frame_range == [0,0] else frame_range[1] - frame_range[0]
    frame_range[1] = frame_range[1] if frame_range[1] != 0 else Nframes

    Lx, Ly = xmax - xmin, ymax - ymin

    # Calculate number of grid squares in x and y direction
    if Ly < Lx:
        y_arr = grid_dataframe['y_center'].diff().values
        idx = int(np.argwhere(y_arr > 0.5)[0])

        grid_len = grid_dataframe['x_center'].loc[idx + 1] - grid_dataframe['x_center'].loc[0]
        Lx += grid_len
        Ly += grid_len
        Nygrid = int(np.round(Ly / grid_len))
        Nxgrid = int(np.round(Nsquares / Nygrid) )
    elif Lx <= Ly:
        grid_len = grid_dataframe['x_center'].loc[1] - grid_dataframe['x_center'].loc[0]
        Lx += grid_len
        Ly += grid_len
        Nxgrid = int(np.round(Lx / grid_len))
        Nygrid = int(np.round(Nsquares / Nxgrid) )

    data= grid_dataframe.loc[:, ['Frame', 'T', 'x_center', 'y_center', 'mean_velocity_X', 'mean_velocity_Y']].values
    data = data[frame_range[0] * Nsquares:frame_range[1] * Nsquares, :]

    x_center = data[:Nsquares, 2].reshape(Nygrid, Nxgrid)
    y_center = data[:Nsquares, 3].reshape(Nygrid, Nxgrid)

    if calculate_average:
        arr_vx = data[:, -2].reshape(Nframes, Nsquares)
        arr_vy = data[:, -1].reshape(Nframes, Nsquares)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            arr_vx = np.nanmean(arr_vx, axis = 0).reshape(Nygrid, Nxgrid)
            arr_vy = np.nanmean(arr_vy, axis = 0).reshape(Nygrid, Nxgrid)
        velocity_field_plotter(x_center, y_center, arr_vx, arr_vy, \
                               title=f'Velocity field averaged over frames {frame_range[0]}:{frame_range[1]}')
        
    else:
    
        if animate:
            anim_wrapper = lambda i, frame, fig, fn: animate_flow_field(i, x_center, y_center, data, fig, fn, Nxgrid, Nygrid)
            animator(data, velocity_field_plotter, (frame_range[0],frame_range[1]), anim_wrapper, inter=frame_interval, show=show)
        else:
            for i in range(frame_range[0], frame_range[1]):
                arr_vx = data[data[:, 0] == i][:, -2].reshape(Nygrid, Nxgrid)
                arr_vy = data[data[:, 0] == i][:, -1].reshape(Nygrid, Nxgrid)
                velocity_field_plotter(x_center, y_center, arr_vx, arr_vy, i=i)
    if show:
        plt.show()

def animate_flow_field(i, X, Y, arr, fig, fn, Nx, Ny):
    # we want a fresh figure everytime
    fig.clf()
    # add subplot, aka axis
    #ax = fig.add_subplot(111)
    # load the frame
    arr_vx = arr[arr[:, 0] == i][:, -2].reshape(Ny, Nx)
    arr_vy = arr[arr[:, 0] == i][:, -1].reshape(Ny, Nx)
    # call the global function
    velocity_field_plotter(X, Y, arr_vx, arr_vy, i=i)
    return
    
def velocity_field_plotter(X,Y, VX, VY, i = 0, title = None):

    plt.quiver(X, Y, VX, VY, units='dots', scale_units='dots' )
    #ax.streamplot(np.unique(x_center), np.unique(y_center), vx, vy, density = 1.5, color = 'black')
    title = f'Velocity field for frame = {i}' if title is None else title
    plt.xlabel(xlabel = 'x')
    plt.ylabel(ylabel = 'y')
    plt.title(title)
    return



def main():

    df = pd.read_csv('../resources/CellSegmentationTracker_spots.csv')
    grid_df = calculate_grid_statistics(df, Ngrid = 15, return_absolute_cell_counts=True, include_features=[], save_csv = False,)

    plot_velocity_field(grid_df, frame_range = [0,5], calculate_average = True, \
                                animate = True, frame_interval = 2000, show = True)

if __name__ == '__main__':
    main()






