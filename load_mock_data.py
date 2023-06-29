# Author: Simon Guldager Andersen
# Date (latest update): 

### SETUP ------------------------------------------------------------------------------------

## Imports:
import os, sys
import zipfile
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats, integrate, interpolate, optimize
from iminuit import Minuit
from PIL import Image
from matplotlib import rcParams
from cycler import cycler


#sys.path.append('C:\\Users\\Simon\\miniconda3\\envs\\cellpose')

#from cellpose import plot, utils

## Change directory to current one
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)
sys.path.append(os.getcwd())

## 
im_path1 = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\valeriias mdck data for simon\\20.09.22"
im_path2 = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\valeriias mdck data for simon\\24.08.22"
im_path3 = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\16.06.23 stretch data for simon (2)"

### FUNCTIONS ----------------------------------------------------------------------------------

def get_imlist(path, format = '.jpg'):

    """
    returns a list of filenames for all png images in a directory
    """
    return [os.path.join(path,f) for f in os.listdir(path) if f.endswith(format)]
def get_imlist_from_zip(path):
    """
    Construct a str list of all images from a zip.file. Each entry corresponds to a picture
    To open a picture, say picture k, with pil, write Image.open(im_list[k]).show()

    path: is the full path of the zip-file.
    """
    imzip = zipfile.ZipFile(path)
    infolist = imzip.infolist()

    #initialize image list
    im_list = [] 
    
    #add pictures to list
    for f in infolist:
        ifile = imzip.open(f)
        im_list.append(ifile)

    return im_list

def generate_info_dict(image, supress_filename = False):
    if supress_filename:
        info_dict = {}
    else:
        info_dict = {"Filename": image.filename}

    info_dict0 = {
    "Image Size": image.size,
    "Image Height": image.height,
    "Image Width": image.width,
    "Image Format": image.format,
    "Image Mode": image.mode,
    }
    info_dict.update(info_dict0)
    return info_dict


def images_to_data_list (image_list):
    """
    Takes a list of images, potentially of different sizes, and returns a list of pixel data, in which each entry is a flattened array containing
    all pixel values. It also returns a list in which each entry contains the shape of a picture in the format (pixel height, pixel width,
    number of integers used to describe 1 pixel) [provided that more than 1 byte is used]
    If the picture has the pixel area NxN, each entry array of the returned list will have length (pixel area x bytes/pixel) ,where bytes/pixel is number of
    integers used to decribe a pixel (1 for greyscale, 3 for 24 bit pictures). 
    
    """
    #Find the number of images in the list
    no_images = len(image_list)
    #Store the shape of each image 
    im_shape_list = []
    #Initialize picture list
    pic_list = []

    for k in range(no_images):
        im = np.array(Image.open(image_list[k]))
        im_shape_list.append(im.shape)
        pic_list.append(im.flatten())
    
    return pic_list, im_shape_list
def images_to_data_arr (image_list, pixel_shape, color = True):

    """
    Takes a list of images, potentially of different sizes, reshapes them to pixel_shape, 
    and returns an array of pixel data, in which each row is a flattened array containing
    all pixel values. 
    The returned array has dmensions (No. of images) x (pixel_shape * bytes/pixel)
    params:
        image_list: a list of images
        pixel_shape: list holding dimensions that all pictures will be resized to in the format 
                    [pixel height, pixel width]
        color: Boolean. If true, the pictures will be assumed to be 24 bit RGB images. 
                        If False, the pictures will be assumed to be 8 bit grey images

    """
    #Find no. of pictures
    no_pictures = len(image_list)

    # Find length of each row
    if color:
        bytes_per_pixel = 3
    else:
        bytes_per_pixel = 1
    
    arr = np.empty([no_pictures, bytes_per_pixel * pixel_shape[0] * pixel_shape[1]])

    #Build matrix arr, with all pixels from one picture per row
    for image in range(no_pictures):
        im_resized = Image.open(image_list[image]).resize((pixel_shape[0], pixel_shape[1]))
        arr[image] = np.array(im_resized).flatten()
    return arr
def convert_to_grey_scale(image_list, pixel_shape): 
    """
    Takes a list of images, potentially of different sizes, reshapes them to pixel_shape, converts them to greyscale,
    and returns an array of pixel data, in which each row is a flattened array containing
    all pixel values. 
    The returned array has dmensions (No. of images) x (pixel_shape * bytes/pixel)
    params:
        image_list: a list of images
        pixel_shape: list holding dimensions that all pictures will be resized to in the format 
                    [pixel height, pixel width]
    """
    n_pictures = len(image_list)

    grey_arr = np.empty([n_pictures, pixel_shape[0] * pixel_shape[1] ])

    for image in range(n_pictures):
        im = Image.open(image_list[image]).convert('L').resize((pixel_shape[0], pixel_shape[1]))
        grey_arr [image] = np.array(im).flatten()
    return grey_arr



### MAIN ---------------------------------------------------------------------------------------

# Set plotting style
sns.set_theme()
sns.set_style("darkgrid")
sns.set_context("paper") #Possible are paper, notebook, talk and poster
rcParams['lines.linewidth'] = 2 
rcParams['axes.titlesize'] =  18
rcParams['axes.labelsize'] =  18
rcParams['xtick.labelsize'] = 12
rcParams['ytick.labelsize'] = 12
rcParams['legend.fontsize'] = 15
rcParams['font.family'] = 'serif'
rcParams['figure.figsize'] = (9,6)
rcParams['axes.prop_cycle'] = cycler(color = ['teal', 'navy', 'coral', 'plum', 'purple', 'olivedrab',\
         'black', 'red', 'cyan', 'yellow', 'khaki','lightblue'])
np.set_printoptions(precision = 5, suppress=1e-10)

def main():
    ## OBJECTIVE: Load mock data images, show them, and convert them to arrays
    show_ims1, show_ims2, show_ims3 = False, False, True

    im_list1 = get_imlist(im_path1, format = '.tif')
    im_list2 = get_imlist(im_path2, format = '.tif')
    im_list3 = get_imlist(im_path3, format = '.tif')

    im0 = Image.open(im_list1[0])
    im1 = Image.open(im_list1[1])

    im0_dict = generate_info_dict(im0, supress_filename = True)
    im0_arr = np.array(im0)
   
    print(im0_dict['Image Size'][0])
    ## imsize
    target_width, target_height = 698, 648
    im0 = im0.resize((target_width, target_height))
    #im0 = im0.resize((1024, 1024))#, Image.Resampling.LANCZOS)
    print (im0.size)
    im0.save(im_path1 + "\\im0_698.tif")

    im1_list = [None] * len(im_list1)
    im2_list = [None] * len(im_list2)
    im3_list = [None] * len(im_list3)

    if show_ims1:
        for i, im in enumerate(im_list1):
            im1_list[i] = Image.open(im).resize((target_width, target_height))
            file_name = f"im{i}"
       #     print(im_path1 + f"_698\\{file_name}.tif")
            im1_list[i].save(im_path1  + f"_{target_width}x{target_height}\\{file_name}.tif")
            im1_list[i].show()
    if show_ims2:
        for i, im in enumerate(im_list2):
            im2_list[i] = Image.open(im).resize((target_width, target_height))
            file_name = f"im{i}"
            im2_list[i].save(im_path2  + f"_{target_width}x{target_height}\\{file_name}.tif")
            im2_list[i].show()
    if show_ims3:
        for i, im in enumerate(im_list3):
            im3_list[i] = Image.open(im).resize((target_width, target_height))
            file_name = f"im{i}"
            im3_list[i].save(im_path3  + f"_{target_width}x{target_height}\\{file_name}.tif")
            im3_list[i].show()

    if 0:
        ## Try loading cellpose output
        output_path = "C:\\Users\\Simon\\Documents\\Uni\\Speciale\\Mock data\\valeriias mdck data for simon\\valeriias mdck data for simon\\20.09.22\\im0_seg.npy"

        dat = np.load(output_path, allow_pickle=True).item()

        #fig = plt.subplots()
    #  mask_RGB = plot.mask_overlay(dat['img'], dat['masks'], colors=np.array(dat['colors']))
        plt.imshow(dat['img'])
    # print(dat['masks'].shape)
        mask = np.argwhere(dat['masks'] != 0)

        outlines = utils.outlines_list(dat['masks'])

        #print(type(outlines))
    # print(outlines[0].shape)
        #print(dat['masks'][600:620,500:510])
        numbers = 0
        #print(dat['img'].shape)

        if 1:
            # plot image with outlines overlaid in red
            outlines = utils.outlines_list(dat['masks'])
            plt.imshow(dat['img'])
            for o in outlines[0:10]:
                Nskip = 4
                range = np.arange(0, o.shape[0], Nskip)
                numbers += len(range) * 2
                plt.plot(o[range,0], o[range,1],)# color='r')

        print("Numbers in image and outlines: ", dat['img'].shape[0] * dat['img'].shape[1], numbers)

        plt.show()

if __name__ == '__main__':
    main()
