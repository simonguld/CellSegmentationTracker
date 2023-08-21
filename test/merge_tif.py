#!/usr/bin/env python


im_dir = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\valeriias mdck data for simon\\24.08.22_698x648"
img_out_name = im_dir + "\\merged.tif"

def get_imlist(path, format = '.jpg'):

    """
    returns a list of filenames for all png images in a directory
    """
    return [os.path.join(path,f) for f in os.listdir(path) if f.endswith(format)]


def merge_tiff(im_dir, img_out_name):
    img_in_names = get_imlist(im_dir, format = '.tif')
    if not img_in_names:
        print("No image!!!")
        return -1
    print("File names: {:}".format(img_in_names))
    t = tifftools.read_tiff(img_in_names[0])
    for img_in_name in img_in_names[1:]:
        t["ifds"].extend(tifftools.read_tiff(img_in_name)["ifds"])

    if os.path.isfile(img_out_name):
        os.unlink(img_out_name)
    tifftools.write_tiff(t, img_out_name)

    imgo = imread(img_out_name)
    imgs = [imread(e) for e in img_in_names]

    print("Output file details: {:}, {:}".format(imgo.shape, imgo.dtype))
    print("Input files details: {:}".format([(img.shape, img.dtype) for img in imgs]))

    if len(img_in_names) in (3, 4):  # @TODO:
        imgo = imgo.transpose((2, 0, 1))


def main():
    merge_tiff(im_dir, img_out_name)


if __name__ == "__main__":    main()

