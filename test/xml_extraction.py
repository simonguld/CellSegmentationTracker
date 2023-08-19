# Author: Simon Guldager Andersen
# Date (latest update): 

# This script contains all the helper functions needed in CellSegmentationTracker.py


## IMPORTS:
import os
import tifftools
import numpy as np
import pandas as pd
import xml.etree.ElementTree as et
from skimage.io import imread
from skimage.transform import resize
from PIL import Image

pd.set_option('mode.chained_assignment', None)

## FUNCTIONS:


def trackmate_xml_to_csv(trackmate_xml_path, include_features_list = None, get_track_features=False, get_edge_features=False):
    """
    This function is a modification of Hadrien Mary's function from the pytrackmate python package
    ...
    Copyright (c) 2018, Hadrien Mary
    All rights reserved.    
    ... and is used here in modified form in accordance with the license.

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
    df_spots = pd.DataFrame([])
    objects = []
    for frame in spots.findall('SpotsInFrame'):
        for spot in frame.findall('Spot'):
            single_object = []
            for label in features:
                single_object.append(spot.get(label))
            objects.append(single_object)

    df_spots = pd.DataFrame(objects, columns=features)
    df_spots = df_spots.astype(float)

    # Apply initial filtering
    initial_filter = root.find("Settings").find("InitialSpotFilter")

    df_spots = filter_spots(df_spots,
                         name=initial_filter.get('feature'),
                         value=float(initial_filter.get('value')),
                         isabove=True if initial_filter.get('isabove') == 'true' else False)

    # Apply filters
    spot_filters = root.find("Settings").find("SpotFilterCollection")

    for spot_filter in spot_filters.findall('Filter'):

        df_spots = filter_spots(df_spots,
                             name=spot_filter.get('feature'),
                             value=float(spot_filter.get('value')),
                             isabove=True if spot_filter.get('isabove') == 'true' else False)

    df_spots = df_spots.loc[:, object_labels.keys()]
    df_spots.columns = [object_labels[k] for k in object_labels.keys()]
    df_spots['TRACK_ID'] = np.arange(df_spots.shape[0])

    # Get track IDs
    filtered_track_ids = [int(track.get('TRACK_ID')) for track in root.find('Model').find('FilteredTracks').findall('TrackID')]

    label_id = 0
    df_spots['TRACK_ID'] = np.nan

    tracks = root.find('Model').find('AllTracks')

    if get_track_features:
        track_features = root.find('Model').find('FeatureDeclarations').find('TrackFeatures')
        track_features = [c.get('feature') for c in list(track_features)]
  
    objects = []
    for track in tracks.findall('Track'):
        if get_track_features:
            single_object = []
            for feature in track_features:
                single_object.append(track.get(feature))
            objects.append(single_object)
            
        track_id = int(track.get("TRACK_ID"))
        if track_id in filtered_track_ids:

            spot_ids = [(edge.get('SPOT_SOURCE_ID'), edge.get('SPOT_TARGET_ID'), edge.get('EDGE_TIME')) for edge in track.findall('Edge')]
            spot_ids = np.array(spot_ids).astype('float')[:, :2]
            spot_ids = set(spot_ids.flatten())

            df_spots.loc[df_spots["Spot ID"].isin(spot_ids), "TRACK_ID"] = label_id
            label_id += 1

    # Label remaining columns
    single_track = df_spots.loc[df_spots["TRACK_ID"].isnull()]
    df_spots.loc[df_spots["TRACK_ID"].isnull(), "TRACK_ID"] = label_id + np.arange(0, len(single_track))

    # Extract velocities
    df_spots['VELOCITY_X'] = np.nan
    df_spots['VELOCITY_Y'] = np.nan

    for track in np.arange(df_spots['TRACK_ID'].max()):
        idx = df_spots.index[df_spots['TRACK_ID'] == track]
        df_spots_res = df_spots.loc[idx].copy()

        frame_gap = np.hstack([df_spots_res['Frame'].diff().values[1:], df_spots_res['Frame'].diff().values[-1:]])
        velocity_x = np.hstack([df_spots_res['X'].diff().values[1:], df_spots_res['X'].diff().values[-1:]]) / frame_gap
        velocity_y = np.hstack([df_spots_res['Y'].diff().values[1:], df_spots_res['Y'].diff().values[-1:]]) / frame_gap

        df_spots['VELOCITY_X'].iloc[idx] = velocity_x
        df_spots['VELOCITY_Y'].iloc[idx] = velocity_y

    if get_track_features:
        df_tracks = pd.DataFrame(objects, columns=track_features).astype(float)
        return df_spots, df_tracks
    else:
        return df_spots

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
    return


def main():
    xml_path = "C:\\Users\\Simon Andersen\\Projects\\Projects\\CellSegmentationTracker\\resources\\FakeTracks.xml"

    csv_path = "resources\\Trackmate.csv"


    # Calculate track features:



    root = et.fromstring(open(xml_path).read())

    track_features = root.find('Model').find('FeatureDeclarations').find('TrackFeatures')
    track_features = [c.get('feature') for c in list(track_features)]
    tracks = root.find('Model').find('AllTracks')

    objects = []
    for track in tracks.findall('Track'):
        print(track.get('TRACK_ID'))
        single_object = []
        for feature in features:
            single_object.append(track.get(feature))
        objects.append(single_object)
    df_tracks = pd.DataFrame(objects, columns=features).astype(float)

    print(df_tracks.info())
    print(df_tracks.head(10))


    if 0:
        df_tracks = pd.DataFrame([])
        objects = []
        for track in tracks.findall('Track'):
            single_object = []
            for feature in track:
                single_object.append(track.get(feature))
            objects.append(single_object)
            print(objects)

        df_tracks = pd.DataFrame(objects, columns=features).astype(float)


        

  

    if 0:
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
        trajs['TRACK_ID'] = np.arange(trajs.shape[0])

        # Get tracks
        if get_tracks:
            filtered_track_ids = [int(track.get('TRACK_ID')) for track in root.find('Model').find('FilteredTracks').findall('TrackID')]

            label_id = 0
            trajs['TRACK_ID'] = np.nan

            tracks = root.find('Model').find('AllTracks')
            for track in tracks.findall('Track'):

                track_id = int(track.get("TRACK_ID"))
                if track_id in filtered_track_ids:

                    spot_ids = [(edge.get('SPOT_SOURCE_ID'), edge.get('SPOT_TARGET_ID'), edge.get('EDGE_TIME')) for edge in track.findall('Edge')]
                    spot_ids = np.array(spot_ids).astype('float')[:, :2]
                    spot_ids = set(spot_ids.flatten())

                    trajs.loc[trajs["Spot ID"].isin(spot_ids), "TRACK_ID"] = label_id
                    label_id += 1

          

     







   

    if 0:

        print("......",df['TRACK_ID'].max())
        

        idx = df.index[df['TRACK_ID'] == track]
        df_res = df.loc[idx].copy()

        print(df_res['Frame'].diff().values[1:])
        print(df_res['Frame'].diff().values[1:].shape)
    
        frame_gap = np.hstack([df_res['Frame'].diff().values[1:], df_res['Frame'].diff().values[-1:]])
        velocity_x = np.hstack([df_res['X'].diff().values[1:], df_res['X'].diff().values[-1:]]) / frame_gap
        velocity_y = np.hstack([df_res['Y'].diff().values[1:], df_res['Y'].diff().values[-1:]]) / frame_gap

        print(frame_gap)


  #  df_res['VELOCITY_X'].iloc[:-1]= df_res['X'].diff()[1:] / frame_gap[1:]
   # df_res['VELOCITY_X'].iloc[-1:] = df_res['X'].diff()[-1:] / frame_gap[-1:]
 
    #df_res['VELOCITY_Y'].iloc[:-1]= df_res['Y'].diff()[1:] / frame_gap[1:]
    #df_res['VELOCITY_Y'].iloc[-1:] = df_res['Y'].diff()[-1:] / frame_gap[-1:]

        df_res['VELOCITY_X'] = velocity_x
        df_res['VELOCITY_Y'] = velocity_y


        print(df_res)
        print(df_res['Frame'].diff())

   # print(df_res.loc[-1])
    if 0:
        print(df_res['VELOCITY_X'].loc[:-1])
        df_res['VELOCITY_X'].iloc[:-1]= df_res['X'].diff()[1:]  # np.array([df_res['X'].diff().iloc[1:], df_res['X'].diff().iloc[-1]])
        df_res['VELOCITY_X'].loc[-1] = df_res['X'].diff()[-1]
    
        df_track = df.loc[idx, ['X','Y','Frame']]
        
        
        print(df_track)
        print(df_track['X'].diff())


    # generate track info:
    if 0:
        root = et.fromstring(open(xml_path).read())
        track_features = root.find('Model').find('FeatureDeclarations').find('TrackFeatures')
        features = [c.get('feature') for c in list(track_features)]

        spots = root.find('Model').find('AllTracks')
        trajs = pd.DataFrame([])
        objects = []
        for frame in spots.findall('Track'):
            print(frame)
            single_object = []
            for label in features:
                single_object.append(frame.get(label))
            objects.append(single_object)

    if 0:
        root = et.fromstring(open(xml_path).read())
        im_settings = root.find('Settings').find('ImageData')
        basic_settings = root.find('Settings').find('BasicSettings')

        im_settings_list = ['filename', 'folder', 'width', 'height', 'nslices', 'nframes', 'pixelwidth', 'pixelheight', 'voxeldepth','timeinterval']
        basic_settings_list = ['xstart', 'xend', 'ystart', 'yend', 'zstart', 'zend', 'tstart', 'tend']
        settings = [im_settings_list, basic_settings_list]

        print_list = []
        for i, obj in enumerate([im_settings, basic_settings]):
            settings_list = settings[i]      
            for  s in settings_list:
                print_list.append(s + f': {obj.get(s)}') 
        print("Image information: \n", print_list)

    """
 

    """


    if 0:
        for spot in frame.findall('Edge'):
            print(spot)
            single_object = []
            for label in features:
                single_object.append(spot.get(label))
            objects.append(single_object)


        trajs = pd.DataFrame(objects, columns=features)
        trajs = trajs.astype(float)

        print(trajs.info())
        print(trajs.head(10))
        print(trajs.describe())


   
if __name__ == "__main__":
    main()

