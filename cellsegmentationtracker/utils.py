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


## FUNCTIONS:

def prepend_text_to_file(file_path, text):
    """
    Prepend text to a file.
    """
    with open(file_path, 'r+') as file:
        content = file.read()
        file.seek(0, 0)
        file.write(text.rstrip('\r\n') + '\n' + content)
        file.close()
    return

def search_and_modify_file(file_path, search_string, modify_string):
    """
    Search for a string in a file and replace it with another string.
    """
    with open(file_path, 'r') as file:
        filedata = file.read()
    filedata = filedata.replace(search_string, modify_string)
    with open(file_path, 'w') as file:
        file.write(filedata)
        file.close()
    return

def get_imlist(path, format = '.jpg'):

    """
    returns a list of filenames for all png images in a directory
    """
    if os.path:
        return [os.path.join(path,f) for f in os.listdir(path) if f.endswith(format)]
    else:
        raise OSError("No such directory: " + path)

def get_target_dimensions(im_path, division_factor = 1):
    """
    Get target dimensions for a single tif image (possibly a stack) to target_width x target_height
    """
    # Find image dimensions
    img = Image.open(im_path)
    if len(img.size) == 2:
        target_width, target_height = int(img.width / division_factor) , int(img.height / division_factor)
    elif img.ndim == 3:
        if img.width > img.shape[2]:
            target_width, target_height = int(img.width / division_factor) , int(img.height / division_factor)
        else:
            target_width, target_height = int(img.height / division_factor) , int(img.shape[2] / division_factor)
    img.close()
    return target_width, target_height

def resize_imlist(im_dir, target_width, target_height, output_dir = None):
    """
    Resize all images in a directory to target_width x target_height
    """
    if output_dir is None:
        output_dir = im_dir  + f"_resized_{target_width}x{target_height}"
        os.mkdir(output_dir)
        print("Images will be saved in: " + output_dir)
    im_list = get_imlist(im_dir, format = '.tif')
    
    for i, f in enumerate(im_list):
        im = Image.open(f).resize((target_width,target_height))
        file_name = f"im{i}"
        im.save(os.path.join(output_dir , f"{file_name}.tif"))
        im.close()
    return

def resize_tif(im_path, division_factor = 1, Nframes = None):
    """
    Resize a single tif image (possibly a stack) to target_width x target_height
    """

    # Find image dimensions
    target_width, target_height = get_target_dimensions(im_path, division_factor = division_factor)

    im_dir = os.path.dirname(im_path)
    im_path_new = im_path.strip(".tif") + f"_resized_{target_width}x{target_height}.tif"
    img = imread(im_path)

    # If only 1 frame is present, simply resize and save
    if img.ndim == 2:
        im = Image.open(im_path).resize((target_width,target_height))
        im.save(im_path_new)
    else:
        # Create directory for resized images
        output_dir = os.path.join(im_dir, f"_resized_{target_width}x{target_height}")
        os.mkdir(output_dir)
        # Split tif file
        info = tifftools.read_tiff(im_path)
        for ifd in info['ifds']:
            tifftools.write_tiff(ifd, im_path_new)
        # Resize all images in directory
        resize_imlist(output_dir, target_width, target_height, Nframes)
        # Merge resized images
        merge_tiff(output_dir, im_path_new, Nframes)
        # Delete resized images
        os.rmdir(output_dir)
        os.rmdir(output_dir + f"_resized_{target_width}x{target_height}")
    return

def split_tiff(im_path, name = None):
    """
    Split a tif file containing multiple frames into multiple tif files
    """
    im_dir = os.path.dirname(im_path)
    if name is None:
        name = f"_split"
    output_dir = im_path.strip(".tif") + name
    os.mkdir(output_dir)
    # Split tif file
    info = tifftools.read_tiff(im_path)
    for i,ifd in enumerate(info['ifds']):
        tifftools.write_tiff(ifd, os.path.join(output_dir, f'im{i}.tif'))
    return

def merge_tiff(im_dir, img_out_path, Nframes = None):
    im_list = get_imlist(im_dir, format = '.tif')
    Nframes = len(im_list) if Nframes is None else Nframes
    if len(im_list) == 0:
        print("No image!!!")
        return -1
    t = tifftools.read_tiff(im_list[0])

    for img_in_name in im_list[1:Nframes]:
        t["ifds"].extend(tifftools.read_tiff(img_in_name)["ifds"])

    if os.path.isfile(img_out_path):
        os.unlink(img_out_path)
    tifftools.write_tiff(t, img_out_path)

    imgo = imread(img_out_path)
    imgs = [imread(e) for e in im_list]


    if len(im_list) in (3, 4):  # @TODO:
        imgo = imgo.transpose((2, 0, 1))

    return

def trackmate_xml_to_csv(trackmate_xml_path, include_spot_features_list = None, calculate_velocities = True,\
                          get_track_features=False, get_edge_features=False):
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

    df_spots = None
    df_tracks = None
    df_edges = None

    root = et.fromstring(open(trackmate_xml_path).read())

    minimal_labels = {'FRAME': 'Frame',
                    'ID': 'Spot ID',
                    'POSITION_X': 'X',
                    'POSITION_Y': 'Y',}
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
    
    # Include only features in include_spot_features_list
    if include_spot_features_list is not None:
        del_list = []
        for el in object_labels.keys():
            if object_labels[el] not in include_spot_features_list:
                del_list.append(el)
        for el in del_list:
            del object_labels[el]
    object_labels.update(minimal_labels)

    spot_features = root.find('Model').find('FeatureDeclarations').find('SpotFeatures')
    spot_features = [c.get('feature') for c in list(spot_features)] + ['ID']

    spots = root.find('Model').find('AllSpots')
    spot_objects = []
    for frame in spots.findall('SpotsInFrame'):
        for spot in frame.findall('Spot'):
            single_object = []
            for label in spot_features:
                single_object.append(spot.get(label))
            spot_objects.append(single_object)

    df_spots = pd.DataFrame(spot_objects, columns=spot_features).astype('float')
 
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
        track_objects = []
    if get_edge_features:
        edge_features = root.find('Model').find('FeatureDeclarations').find('EdgeFeatures')
        edge_features = [c.get('feature') for c in list(edge_features)]
        edge_objects = []


    for track in tracks.findall('Track'):
        if get_track_features:
            single_object = []
            for feature in track_features:
                single_object.append(track.get(feature))
            track_objects.append(single_object)

        if get_edge_features:
            for edge in track.findall('Edge'):
                single_object = []
                for label in edge_features:
                    single_object.append(edge.get(label))
                edge_objects.append(single_object)
            
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
    

    if calculate_velocities:
        # Extract velocities
        df_spots['Velocity_X'] = np.nan
        df_spots['Velocity_Y'] = np.nan

        for track in np.arange(df_spots['TRACK_ID'].max()):
            idx = df_spots.index[df_spots['TRACK_ID'] == track]
            df_spots_res = df_spots.loc[idx].copy()

            frame_gap = np.hstack([df_spots_res['Frame'].diff().values[1:], df_spots_res['Frame'].diff().values[-1:]])
            velocity_x = np.hstack([df_spots_res['X'].diff().values[1:], df_spots_res['X'].diff().values[-1:]]) / frame_gap
            velocity_y = np.hstack([df_spots_res['Y'].diff().values[1:], df_spots_res['Y'].diff().values[-1:]]) / frame_gap

            df_spots['Velocity_X'].iloc[idx] = velocity_x
            df_spots['Velocity_Y'].iloc[idx] = velocity_y

    if get_edge_features:
        df_edges = pd.DataFrame(edge_objects, columns=edge_features).astype(float)
    if get_track_features:
        df_tracks = pd.DataFrame(track_objects, columns=track_features).astype(float)

    return df_spots, df_tracks, df_edges

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

    from PIL import Image
    path = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject"
    new_path = path.strip(".tif") + "_resized.tif"
    path2 = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\t_164_428 - 2023_03_03_TL_MDCK_2000cells_mm2_10FNh_10min_int_frame1_200_cropped.tif"
    path3 = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\16.06.23_stretch_data\\im0.tif"
 


    dir = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\16.06.23_stretch_data\\im0_split"
    dir2 = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\t_164_428 - 2023_03_03_TL_MDCK_2000cells_mm2_10FNh_10min_int_frame1_200_cropped_split"
    dir3 = os.path.join(path3.strip(".tif") + "_split")
    target_width = 698
    target_height = 648
    names = ['epi6000.tif', 'epi2500.tif']
    
    im_p = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\tt"

    merge_tiff(im_p, os.path.join(im_p, "merged.tif"), Nframes = 4)

    if 0:
    
        for i, dirr in enumerate([dir2]):
            im0 = get_imlist(dirr, format = '.tif')[0]
            target_width, target_height = get_target_dimensions(im0, division_factor = 1)
            resize_imlist(dirr, target_width, target_height)
            dirrr = dirr + f"_resized_{target_width}x{target_height}"

            merge_tiff(dirrr, os.path.join(path, f"{names[i]}"), Nframes = 4)


                #split_tiff(im_path)

            #split_tif(im_path)

    if 0:
        img = Image.open(path2)
        img_dir = os.path.dirname(path2)
        # Create directory for resized images
        im_path = path2
        output_dir = im_path.strip(".tif") + f"_split"
        os.mkdir(output_dir)
        # Split tif file
        info = tifftools.read_tiff(im_path)
        for i,ifd in enumerate(info['ifds']):
            tifftools.write_tiff(ifd, os.path.join(output_dir, f'im{i}.tif'))


    if 0:
        for i in range(4):
            try:
                img.seek(i)
                img.save('page_%s.tif'%(i,))
            except EOFError:
                break

    if 0:

        # Create directory for resized images
        im_path = dir4 + ".tif"
        output_dir = dir4 + f"_split_v2"
        os.mkdir(output_dir)
        # Split tif file
        info = tifftools.read_tiff(im_path)
        for i,ifd in enumerate(info['ifds']):
            tifftools.write_tiff(ifd, output_dir + f'im{i}.tif')

        for d in [dir, dir2, dir3]:
            im_path = get_imlist(d, format = '.tif')[0]
            target_width, target_height = get_target_dimensions(im_path, division_factor = 1)
            resize_imlist(d, target_width, target_height)


    if 0:
        #STEP 1: does imp work for Nframes,W,H?

        par_dir = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject"
        dir = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\nuclei"
        target_width = 698
        target_height = 648
        new_dir = dir  + f"_resized_{target_width}x{target_height}"
        os.mkdir(new_dir)
        im_list = get_imlist(dir, format = '.tif')
        
    
        for i, f in enumerate(im_list):
            im = Image.open(f).resize((target_width,target_height))
            file_name = f"im{i}"
            print(os.path.join(new_dir,f"{file_name}.tif"))
            im.save(os.path.join(new_dir , f"{file_name}.tif"))

        merge_tiff(new_dir, par_dir + "\\merged.tif", Nframes = 4)


    if 0:

        imgo = imread(path4)
        print(imgo.shape, type(imgo))

        if 1:
            for p in [path]: #, path2, path3, path4]:
                data = tifffile.imread(p)
                print(data.shape, type(data))
                print("Assuming W x H x Nchannels x Nframes")
                shape = data.shape
                new_shape = (100,100, shape[0])
                resized_data = resize(data, new_shape, anti_aliasing=True)
                print(resized_data.shape)
                tifffile.imwrite(new_path, resized_data, planarconfig='CONTIG')

        #    new_d = tifffile.imread(path5)
        #   print(new_d.shape, type(new_d))
            
    #     resized_data = resize(data, new_shape, anti_aliasing=True)


            new_width, new_height = 100, 120
            #data = tifftools.read_tiff(path)
            #data = imread(path)
        # data = tifftools.read_tiff(path)
        #   print(data.shape, type(data))
        if 0:
            print(np.swapaxes(data, 0, 2).shape)
            new_shape = (new_height, new_width, data.shape[2])
            resized_data = resize(data, new_shape, anti_aliasing=True)
            tifftools.write_tiff(resized_data, path.strip('.tif') + '_resized.tif')
            print(resized_data.shape)


if __name__ == "__main__":
    main()

