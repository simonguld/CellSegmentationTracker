# Author: Simon Guldager Andersen
# Date (latest update): 

### SETUP ------------------------------------------------------------------------------------

## Imports:
import os
import sys
import time
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats, integrate, interpolate, optimize
from scipy.special import sici, factorial
from iminuit import Minuit
from mpl_toolkits.axes_grid1 import make_axes_locatable
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



def main():
    df = pd.read_csv('../resources/epi500_sample_images_spots.csv')

    data = df.loc[:, ['Frame','X', 'Y', 'VELOCITY_X', 'VELOCITY_Y']].values

 
    Nframes = int(data[:,0].max() + 1)
    #Nframes = data

    xmin, xmax = data[:,1].min(), data[:,1].max()
    ymin, ymax = data[:,2].min(), data[:,2].max()

    Lx = xmax - xmin
    Ly = ymax - ymin

    Nygrid = 50

    dy = Ly / Nygrid
    dA = dy**2

    Residual_x = Lx % dy
    Nxgrid = int(np.floor(Lx / dy))
    Nsquares = Nxgrid * Nygrid

    xmin += Residual_x / 2
    xmax -= Residual_x / 2
    Lx = xmax - xmin


   
    columns = ['Frame','Ngrid','xmin', 'xmax', 'ymin', 'ymax', 'density','vx','vy']
    grid_arr = np.zeros([Nframes * Nsquares, len(columns)])

    ### HANDLE NANS!!!

    t1 = time.time()
    for i, frame in enumerate(np.arange(Nframes)):
        arr_T = data[data[:, 0] == frame]
        Ngrid = 0
        
        for j, x in enumerate(np.linspace(xmin, xmax - dy, Nxgrid)):
            mask = (arr_T[:, 1] >= x) & (arr_T[:, 1] < x + dy)
            arr_X = arr_T[mask]
      
            for k, y in enumerate(np.linspace(ymin, ymax - dy, Nygrid)):
                mask = (arr_X[:, 2] >= y) & (arr_X[:, 2] < y + dy)

                density = mask.sum()  / dA
                vx = np.nanmean(arr_X[mask, 3])
                vy = np.nanmean(arr_X[mask, 4])

                grid_arr[frame * Nsquares + Ngrid] = [frame, Ngrid, x, x + dy, y, y + dy, density, vx, vy]
                Ngrid += 1

    t2 = time.time()
    print(f'Time elapsed: {t2 - t1:.2f} s')
    
           
    print(grid_arr[:20])
    print(grid_arr.shape)
    print(grid_arr[-20:]    )
    grid_df = pd.DataFrame(grid_arr, columns = columns)

    print(grid_df.info())
    print(grid_df.head())
    print(grid_df.describe())
    
    

if __name__ == '__main__':
    main()




"""
for j, x in enumerate(np.linspace(xmin, xmax, Nxgrid - 1)):
            interval = pd.IntervalIndex([pd.Interval(x, x + dy, closed='both')])
            mask = df_T['X'].apply(lambda x: x in interval)
            df_X = df_T.loc[mask]
            for k, y in enumerate(np.linspace(ymin, ymax, Nygrid - 1)):
                
                interval = pd.IntervalIndex([pd.Interval(y, y + dy, closed='both')])
                mask = df_X['Y'].apply(lambda x: x in interval)
                df_XY = df_X.loc[mask]
                Nspots = len(df_XY.index)
                density = Nspots / dA
                vx = x + dy / 2
                vy = y + dy / 2
                grid_df.loc[j*Nygrid + k] = [frame, Nspots, x, x + dy, y, y + dy, density, vx, vy]


           
"""





