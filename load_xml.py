
## Imports:
import os, sys
import numpy as np
import seaborn as sns
import pandas as pd
import xml.etree.ElementTree as et

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rcParams
from matplotlib.patches import Ellipse
from cycler import cycler


#from pytrackmate import trackmate_peak_import 

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

path_mock = 'FakeTracks1.xml' 
path_celss = '20.09_tracks.xml' 

### FUNCTIONS -----------------------------------------------------------------



def trackmate_peak_import(trackmate_xml_path, include_features_list = None, get_tracks=True):
    """Import detected peaks with TrackMate Fiji plugin.

    Parameters
    ----------
    trackmate_xml_path : str
        TrackMate XML file path.
    get_tracks : boolean
        Add tracks to label
    """

    root = et.fromstring(open(trackmate_xml_path).read())

    objects = []
    minimal_labels = {'FRAME': 'Frame',
                      'ID': 'Spot ID',}
    object_labels = {'FRAME': 'Frame',
                     'POSITION_T': 'T',
                     'POSITION_X': 'X',
                     'POSITION_Y': 'Y',
                     'POSITION_Z': 'Z',
                     'RADIUS': 'Radius',
                     'VISIBILITY': 'Visibility',
                     'MANUAL_SPOT_COLOR': 'Manual spot color',
                     'ELLIPSE_X0': 'Ellipse center x0',
                        'ELLIPSE_Y0': 'Ellipse center y0',
                     'ELLIPSE_MINOR': 'Ellipse short axis',
                        'ELLIPSE_MAJOR': 'Ellipse long axis',
                        'ELLIPSE_THETA': 'Ellipse angle',
                        'ELLIPSE_ASPECTRATIO': 'Ellipse aspect ratio',
                        'AREA': 'Area',
                        'PERIMETER': 'Perimeter',
                        'CIRCULARITY': 'Circularity',
                        'SOLIDITY': 'Solidity',
                        'SHAPE_INDEX': 'Shape index',
                     'QUALITY': 'Quality',      
                     'MEAN_INTENSITY_CH1': 'Mean intensity ch1',
                     'MEDIAN_INTENSITY_CH1': 'Median intensity ch1',
                     'MIN_INTENSITY_CH1': 'Min intensity ch1',
                     'MAX_INTENSITY_CH1': 'Max intensity ch1',
                     'TOTAL_INTENSITY_CH1': 'Sum intensity ch1',
                    'STD_INTENSITY_CH1': 'Std intensity ch1',
                    'CONTRAST_CH1': 'Contrast ch1',
                    'SNR_CH1': 'Signal/Noise ratio ch1',
                    'ID': 'Spot ID',
                    }
    
    # Include only features in include_features_list
    if include_features_list is not None:
        del_list = []
        for el in object_labels.keys():
            if object_labels[el] not in include_features_list:
                del_list.append(el)
        for el in del_list:
            del object_labels[el]
    object_labels.update(minimal_labels)

    features = root.find('Model').find('FeatureDeclarations').find('SpotFeatures')
    #features = [c.get('feature') for c in features.getchildren()] + ['ID']
    features = [c.get('feature') for c in list(features)] + ['ID']

    spots = root.find('Model').find('AllSpots')
    trajs = pd.DataFrame([])
    objects = []
    for frame in spots.findall('SpotsInFrame'):
        for spot in frame.findall('Spot'):
            single_object = []
            for label in features:
                single_object.append(spot.get(label))
            objects.append(single_object)

    trajs = pd.DataFrame(objects, columns=features)
    trajs = trajs.astype(float)

    # Apply initial filtering
    initial_filter = root.find("Settings").find("InitialSpotFilter")

    trajs = filter_spots(trajs,
                         name=initial_filter.get('feature'),
                         value=float(initial_filter.get('value')),
                         isabove=True if initial_filter.get('isabove') == 'true' else False)

    # Apply filters
    spot_filters = root.find("Settings").find("SpotFilterCollection")

    for spot_filter in spot_filters.findall('Filter'):

        trajs = filter_spots(trajs,
                             name=spot_filter.get('feature'),
                             value=float(spot_filter.get('value')),
                             isabove=True if spot_filter.get('isabove') == 'true' else False)

    trajs = trajs.loc[:, object_labels.keys()]
    trajs.columns = [object_labels[k] for k in object_labels.keys()]
    trajs['label'] = np.arange(trajs.shape[0])

    # Get tracks
    if get_tracks:
        filtered_track_ids = [int(track.get('TRACK_ID')) for track in root.find('Model').find('FilteredTracks').findall('TrackID')]

        label_id = 0
        trajs['label'] = np.nan

        tracks = root.find('Model').find('AllTracks')
        for track in tracks.findall('Track'):

            track_id = int(track.get("TRACK_ID"))
            if track_id in filtered_track_ids:

                spot_ids = [(edge.get('SPOT_SOURCE_ID'), edge.get('SPOT_TARGET_ID'), edge.get('EDGE_TIME')) for edge in track.findall('Edge')]
                spot_ids = np.array(spot_ids).astype('float')[:, :2]
                spot_ids = set(spot_ids.flatten())

                trajs.loc[trajs["Spot ID"].isin(spot_ids), "label"] = label_id
                label_id += 1

        # Label remaining columns
        single_track = trajs.loc[trajs["label"].isnull()]
        trajs.loc[trajs["label"].isnull(), "label"] = label_id + np.arange(0, len(single_track))

    return trajs


def filter_spots(spots, name, value, isabove):
    if isabove:
        spots = spots[spots[name] > value]
    else:
        spots = spots[spots[name] < value]

    return spots

def print_all_possible_spot_features():

    object_labels = {'FRAME': 'Frame',
                     'POSITION_T': 'T',
                     'POSITION_X': 'X',
                     'POSITION_Y': 'Y',
                     'POSITION_Z': 'Z',
                     'RADIUS': 'Radius',
                     'VISIBILITY': 'Visibility',
                     'MANUAL_SPOT_COLOR': 'Manual spot color',
                     'ELLIPSE_X0': 'Ellipse center x0',
                        'ELLIPSE_Y0': 'Ellipse center y0',
                     'ELLIPSE_MINOR': 'Ellipse short axis',
                        'ELLIPSE_MAJOR': 'Ellipse long axis',
                        'ELLIPSE_THETA': 'Ellipse angle',
                        'ELLIPSE_ASPECTRATIO': 'Ellipse aspect ratio',
                        'AREA': 'Area',
                        'PERIMETER': 'Perimeter',
                        'CIRCULARITY': 'Circularity',
                        'SOLIDITY': 'Solidity',
                        'SHAPE_INDEX': 'Shape index',
                     'QUALITY': 'Quality',      
                     'MEAN_INTENSITY_CH1': 'Mean intensity ch1',
                     'MEDIAN_INTENSITY_CH1': 'Median intensity ch1',
                     'MIN_INTENSITY_CH1': 'Min intensity ch1',
                     'MAX_INTENSITY_CH1': 'Max intensity ch1',
                     'TOTAL_INTENSITY_CH1': 'Sum intensity ch1',
                    'STD_INTENSITY_CH1': 'Std intensity ch1',
                    'CONTRAST_CH1': 'Contrast ch1',
                    'SNR_CH1': 'Signal/Noise ratio ch1',
                    'ID': 'Spot ID',
                    }
    
    print(object_labels.values())

def draw_ellipses(ax, dataframe):
    for index, row in dataframe.iterrows():

        ellipse = Ellipse(xy=(row['X'] - row['Ellipse center x0'], row['Y'] - row['Ellipse center y0']),
                          width=row['Ellipse long axis'],
                          height=row['Ellipse short axis'],
                          angle=row['Ellipse angle'],
                          edgecolor='r',
                          fc='None',
                          lw=1)
        ax.add_patch(ellipse)

def animate( plotter, dataframe, t_range, inter=400, plot_masks=True, show=True, repeat=False):
    """Show a frame-by-frame animation.
    """
   
    if len(t_range) == 0:
        t_range = (0, dataframe['T'].max())

    fig, ax = plt.subplots()

    # the local animation function
    def plotter_wrapper(i):
        # we want a fresh figure everytime
       # fig.clf()
        ax.clear()
        # add subplot, aka axis
        #ax = fig.add_subplot(111)
        # load the frame
        plotter(ax, dataframe, i, plot_masks=True)

    anim = FuncAnimation(fig, plotter_wrapper,
                             frames=t_range,
                             interval=inter, blit=False, repeat=repeat)
    plt.show()
    return

    return anim

def plot_spots(ax, spots, t, xlim = None, ylim = None, label=None, plot_masks = False):
    spots_t = spots[spots['T'] == t]
    ax.plot(spots_t['X'], spots_t['Y'], 'o', label=label)
    ax.set(xlabel='X', ylabel='Y', title='Spots for frame = {}'.format(t+1))
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    if plot_masks:
        draw_ellipses(ax, spots_t)
    plt.gca().invert_yaxis()



### MAIN -----------------------------------------------------------------------


def main():
    include_features = ['Frame', 'T', 'X', 'Y', 'Ellipse long axis', 'Ellipse short axis', 'Ellipse angle', 'Ellipse center x0', 'Ellipse center y0', 'Spot ID']

    #tree = ET.parse('FakeTracks1.xml')
    #spots = tree.getroot()

    spots = trackmate_peak_import(path_mock, include_features_list=include_features, get_tracks=False)
    

    # spots = trackmate_peak_import(path_mock)
    print(spots.head())
    print(spots.info())
    print(spots.describe())

    time_range = np.arange(spots['T'].min(), spots['T'].max() + 1)
    xlim, ylim = (0,120), (0,120)
    plotter = lambda ax, spots, t, plot_masks: plot_spots(ax, spots, t, xlim, ylim, plot_masks=plot_masks)

    anim = animate(plotter, spots, t_range=time_range[:25], inter=300, plot_masks=True,)
    # fig, ax = plt.subplots()
    # for t in time_range:
    #     fig, ax = plt.subplots()
    #     plot_spots(ax, spots, t, xlim, ylim, label='t = {}'.format(time_range[0]), plot_masks=True)
    #     plt.show()
    # for t in [0,]:
    #     spots_t = spots[spots['T'] == t]
    #     ax.plot(spots_t['X'], spots_t['Y'], 'o', label='t = {}'.format(t),)
    #     draw_ellipses(ax, spots_t)
    #     ax.legend()
    #     ax.set(xlabel='X', ylabel='Y')
    #     plt.gca().invert_yaxis()
    

if __name__ == '__main__':
    main()
