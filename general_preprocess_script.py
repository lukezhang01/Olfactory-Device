# =================================
# Ethen Sun 29 June 2024
# =================================
# Instructions: prepare the config file and run the script.

# Changes:
# 29 June 2024 Performance improvement
#   config file can now optionally be specified as -c parameter in terminal
#   recentering code separated out into separate function
#   stack function now saves entire stack in an npy file for analysis

import numpy as np
from datetime import datetime
from scipy.signal import fftconvolve
import os
import argparse
import json
from tqdm import tqdm
from astropy.io import fits
from astropy.modeling.functional_models import Lorentz1D
import cv2
#import imageio.v3 as iio
import matplotlib.pyplot as plt

#todo: job log, bitpix

# Interpret config ------------------------------------
parser = argparse.ArgumentParser(description = "Eat that data")
parser.add_argument('-c',
                    dest    = 'configfile',
                    type    = str,
                    default = 'hama-processing-config.json',
                    help    = 'Name of configuration file.'
                    )
args = parser.parse_args()
conf = json.load(open(args.configfile))
if conf['job_id'] is None:
    conf['job_id'] = datetime.utcnow().strftime('%Y%m%d-%H%M%S')
images_location = [conf['images_location'] + path for path in conf['dataset']]
flats_location = [conf['images_location'] + path for path in conf['flats_location']]
darks_location = [conf['images_location'] + path for path in conf['darks_location']]
biases_location = [conf['images_location'] + path for path in conf['biases_location']]

# Define gaussian and lorentzian profiles for convolution use
recenter_filter = conf['alignment']
if isinstance(conf['alignment'], list):
    x, y = np.meshgrid(np.linspace(-1, 1, conf['alignment']["box_size"]+1), np.linspace(-1, 1, conf['alignment']["box_size"]+1))
    d = np.sqrt(x*x + y*y)
    if conf['alignment']['function_type'] == "gaussian":
        recenter_filter = np.exp(-((d)**2 / (2.0 * conf['alignment']['width']**2)))
    elif conf['alignment']['function_type'] == "lorentzian":
        recenter_filter = Lorentz1D(fwhm=conf['alignment']['width'])(d)
    elif conf['alignment']['function_type'] == "first_image":
        recenter_filter = conf['alignment']['box_size']

# functions-----------------------------------------------------------------
def get_file_list(location, file_type=lambda x: x==x):
    # Inputs:
    # folder location or list of folder locations
    # function for matching filename, e.g. lambda x: x.split("_")[0] == "Light"
    # Returns list of matching FITS images
    files = []
    if isinstance(location, str):
        location = [location]
    for sub_location in location:
        files += [sub_location + y for y in sorted(list(filter(lambda x: x.lower().endswith(".fits") \
                                                                      and not x.lower().startswith(".") \
                                                                      and not x.lower().startswith("hide_") \
                                                                      and not x.lower().startswith("output_") \
                                                                           and file_type(x),
                                   next(os.walk(sub_location))[2])))]
    #print(files)
    return files
def debayer(image, debayer_matrix=None, depth=65535):
    # Take a bayered image and return an rgb array
    # converts to uint16 because debayer requires it, and converts back to float
    if not debayer_matrix or len(image.shape) > 2:
        return image
    elif debayer_matrix == "RGGB":
        debayer_matrix = cv2.COLOR_BAYER_BGGR2RGB
    else:
        print("huh")
    scalefactor = depth/np.max(image)
    debayered_image = cv2.cvtColor(np.uint16(image*scalefactor), debayer_matrix)
    #normalised_image = cv2.normalize(debayered_image, None, alpha=0, beta=1, norm_type=cv2.NORM_MINMAX,
    #                                 dtype=cv2.CV_32F)
    #imgplot = plt.imshow(normalised_image)
    return debayered_image/scalefactor
def identify_center(grey_map, recenter_filter):

    g = fftconvolve(grey_map - np.median(grey_map), recenter_filter, mode="same")
    # coordinates of detected target
    x, y = np.unravel_index(g.argmax(), g.shape)
    x -= len(recenter_filter) // 2
    y -= len(recenter_filter) // 2

    return x, y
def recenter(stack, centers, align="maximum"):
    if align == "minimum":
        print("Aligning stack (minimum size)")
        # right edge needs to be moved in by (maximum x - midpoint x)
        # lower edge needs to be moved in by (maximum y - midpoint y)
        crop_dimensions = [np.min([dimensions[0] - np.max(centers[:, 0]), np.min(centers[:, 0])]) * 2 - 2,
                           np.min([dimensions[1] - np.max(centers[:, 1]), np.min(centers[:, 1])]) * 2 - 2]
        lefts = centers[:, 0] - crop_dimensions[0] // 2
        tops = centers[:, 1] - crop_dimensions[1] // 2
        rights = lefts + crop_dimensions[0] + 1
        bottoms = tops + crop_dimensions[1] + 1

        aligned_stack = np.zeros_like(stack[:, :crop_dimensions[0] + 1, :crop_dimensions[1] + 1])
        for n in tqdm(range(len(aligned_stack))):
            # print(aligned_stack[n].shape, stack[n].shape, lefts[n], rights[n], tops[n], bottoms[n], centers[n], aligned_stack[n].shape, stack[n][lefts[n]:rights[n], tops[n]:bottoms[n]].shape)
            aligned_stack[n] = stack[n][lefts[n]:rights[n], tops[n]:bottoms[n]]
    elif align == "maximum":
        crop_dimensions = [dimensions[0] + np.max(centers[:, 0]) - np.min(centers[:, 0]),
                           dimensions[1] + np.max(centers[:, 1]) - np.min(centers[:, 1])]
        lefts = np.max(centers[:, 0]) - centers[:, 0]
        tops = np.max(centers[:, 1]) - centers[:, 1]
        rights = lefts + dimensions[0]
        bottoms = tops + dimensions[1]

        print("Aligning stack (maximum size)", crop_dimensions)
        align_shape = stack.shape
        align_shape = np.array(align_shape)
        align_shape[1:3] = crop_dimensions
        # offset of 128 bytes is room to store the header later, to save the aligned stack as npy
        aligned_stack = np.memmap("stack_overflow.npy", mode="w+", shape=tuple(align_shape), dtype=stack.dtype, offset=128)
        for n in tqdm(range(len(stack))):
            # print(aligned_stack[n].shape, stack[n].shape, lefts[n], rights[n], tops[n], bottoms[n], centers[n], dimensions, aligned_stack[n, lefts[n]:rights[n], tops[n]:bottoms[n]].shape)
            aligned_layer = np.nan * np.ones(align_shape[1:])
            aligned_layer[lefts[n]:rights[n], tops[n]:bottoms[n]] = stack[n]
            aligned_stack[n] = aligned_layer
        print("aligned stack:", aligned_stack.shape)
        return aligned_stack, align_shape
def simplestack(images, mode="median", flat=1, dark=0, bias=0):
    # reads every file in a provided list of images and stacks them
    # kwargs:
    #    mode:
    #    dark, bias:    if provided, subtract these before stacking
    # Does not debayer
    stack = False
    for image in images:
        path = image
        with fits.open(path) as hdul:
            hdu = hdul[0]
            header = hdu.header
            image_data = hdu.data.astype("float64")
            image_data = image_data - dark - bias
            image_data = image_data / flat
            if mode == "premedian":
                # Normalize the image before stacking: useful for flats with different brightnesses
                image_data = image_data / np.median(image_data)
            if isinstance(stack, bool):
                stack = image_data
            else:
                stack = np.dstack((stack, image_data))
    if len(images) == 1:
        return stack
    elif mode == "median" or mode == "premedian":
        return np.median(stack, axis=2)
    else:
        print("Unknown stacking mode. Defaulting to median.")
        return np.median(stack, axis=2)
def stack(conf, images, recenter_filter=None, mode="median", flat=1, dark=0, bias=0):
    """reads every file in a provided list of images and stacks them.
    args:
        images: list of file paths
        recenter_filter: a representative 2D image of the target that should be centered in the stack
            False = stack as is and do not align
            n = use central nxn px of 1st image as the reference image
    kwargs:
        mode:
        flat, dark, bias:    if provided, normalize/subtract these before stacking
    returns aligned stack and recenter filter.
    saves a stack.npy file in the same directory that contains the entire stack.
    debayers if image is bayered, in which case debayer_matrix must be supplied.
    first axis is stack axis, not last.
    maximum or minimum flag inside function sets if output should be overlap area only or contain blank areas"""
    centers = []
    lightcurves = []
    dimensions = [1, 1]
    # force memory saving measures if analysing more than 30 images.
    if (len(images) > 30):
        conf['memorysaving'] = True

    for n in tqdm(range(len(images))):
        path = images[n]
        with fits.open(path) as hdul:
            hdu = hdul[0]
            header = hdu.header
            # this is same as image_data.shape
            dimensions[1] = header['NAXIS1']
            dimensions[0] = header['NAXIS2']
            x,y = dimensions[0]//2, dimensions[1]//2
            # float32 prevents negative overflow errors from negative dark, even if flat division is not used
            image_data = hdu.data.astype("float32")
            image_data = image_data - dark - bias
            image_data = image_data/flat
            grey_map = image_data
            if conf['debayer_matrix']:
                image_data = debayer(image_data, debayer_matrix=conf['debayer_matrix'])
                grey_map = np.sum(image_data, axis=2)

            # things to do on first image read
            if n==0:
                # grab the correct recenter filter if necessary
                if isinstance(recenter_filter, int):
                    # Use the central recenter_filter x recenter_filter pixels of the first image as the kernel
                    recenter_filter = grey_map[dimensions[0] // 2 - recenter_filter // 2:dimensions[
                                                                                             0] // 2 + recenter_filter // 2 + 1,
                                      dimensions[1] // 2 - recenter_filter // 2: dimensions[
                                                                                     1] // 2 + recenter_filter // 2 + 1]

                # setup empty stack
                stack_shape = [len(images)] + list(image_data.shape)
                if conf['memorysaving']:
                    # use disk for stack; make memmap and leave 128 bytes for storing header
                    if recenter_filter:
                        stack = np.memmap("overflow.npy", mode="w+", shape=tuple(stack_shape), dtype=image_data.dtype,
                                      offset=128)
                    else:
                        stack = np.memmap("stack_overflow.npy", mode="w+", shape=tuple(stack_shape), dtype=image_data.dtype,
                                      offset=128)
                else:
                    # use RAM for stack
                    stack = np.empty(stack_shape, dtype=image_data.dtype)

            # Store the brightness of central 200 px by 200 px and the overall brightness for plotting
            lightcurves += [[np.median(grey_map[x - 100:x + 101, y - 100:y + 101]), np.median(grey_map)]]
            if recenter_filter:
                x, y = identify_center(image_data, recenter_filter)
                if conf['diagnostics']:
                    if (not os.path.exists(images_location + '/diagnostics')):
                        os.mkdir(images_location + '/diagnostics')
                    newpath = path[:path.rfind('/')] + '/diagnostics' + path[path.rfind('/'):]
                    plt.imshow(fftconvolve(grey_map, recenter_filter, mode="valid"))
                    plt.plot(y - (len(recenter_filter) // 2), x - (len(recenter_filter) // 2), "bo", markerfacecolor="none")
                    plt.savefig(newpath + "_convolution.png")
                    plt.close()
                    plt.imshow(grey_map)
                    plt.plot(y, x, "bo", markerfacecolor="none")
                    plt.savefig(newpath + "_detection.png")
                    plt.close()

            stack[n] = image_data
            if conf['memorysaving'] and (n % 10) == 0:
                    # flush RAM every 10 files
                    stack.flush()
            centers.append(np.array([x, y]))
    print("Done reading images")
    #--------------------------------------------------------------
    if recenter_filter:
        aligned_stack, align_shape = recenter(stack, centers, align="maximum")
    else:
        aligned_stack = stack
        align_shape = stack.shape
        align_shape = np.array(align_shape)
        align_shape[1:3] = dimensions
    aligned_stack.flush()
    npy_header = np.lib.format.header_data_from_array_1_0(aligned_stack)
    npy_filename = aligned_stack.filename
    with open(npy_filename, 'r+b') as f:
        np.lib.format.write_array_header_1_0(f, npy_header)
    #---------------------------------------------------------------
    print("Collapsing stack")
    if len(images) == 1:
        output = aligned_stack[0]
    elif conf['memorysaving']:
        print("permuting axes")
        output = np.zeros(align_shape[1:])
        # Now first axis is rows
        aligned_stack = np.moveaxis(aligned_stack, 0, 1)
        # print("Collapsing aligned stack...")
        # aligned_stack = np.nanmedian(aligned_stack, axis=0, overwrite_input=True)
        print("Collapsing aligned stack by row...")
        for n in tqdm(range(len(aligned_stack))):
            output[n] = np.nanmedian(aligned_stack[n], axis=0)
            output[n] = np.nan_to_num(output[n])
        print("done nanmedian")
        # output = np.nan_to_num(aligned_stack)
    else:
        output = np.median(stack, axis=0)
    del aligned_stack
    lightcurves = np.array(lightcurves)
    stacked = {"data": output,
               "recenter_filter": recenter_filter,
               "lightcurves": lightcurves}
    #os.rename(npy_filename, 'stack.npy')
    return stacked

# Get files and find calibrations ----------------------------------------------------------
light_files = get_file_list(images_location)
# get dimensions of images
dimensions = [1,1]
with fits.open(light_files[0]) as hdul:
    hdu = hdul[0]
    header = hdu.header
    dimensions[1] = header['NAXIS1']
    dimensions[0] = header['NAXIS2']

print("Getting biases")
bias_files = get_file_list(biases_location, lambda x: x.split("_")[0] == "MasterBias")
if len(bias_files) == 0:
    print("    No master bias found. Attempting to stack biases.")
    bias_files = get_file_list(biases_location)
    if len(bias_files) == 0:
        print("    No bias found")
        master_bias = np.zeros(dimensions)
    else:
        print("    Stacking biases")
        master_bias = simplestack(bias_files)
        hdu = fits.PrimaryHDU(master_bias)
        hdu.writeto(biases_location[0] + "/MasterBias_" + conf['job_id'] + ".fits", overwrite=True)
else:
    print("Master bias found.")
    master_bias = simplestack(bias_files)

print("Getting darks")
dark_files = get_file_list(darks_location, lambda x: x.split("_")[0] == "MasterDark")
if len(dark_files) == 0:
    print("    No master dark found. Attempting to stack darks.")
    dark_files = get_file_list(darks_location)
    if len(dark_files) == 0:
        print("    No darks found")
        master_dark = np.zeros(dimensions)
    else:
        print("    Stacking darks")
        master_dark = simplestack(dark_files, bias=master_bias)
        hdu = fits.PrimaryHDU(master_dark)
        hdu.writeto(darks_location[0] + "/MasterDark_" + conf['job_id'] + ".fits", overwrite=True)
else:
    print("Master dark found.")
    master_dark = simplestack(dark_files, bias=master_bias)

print("Getting flats")
flat_files = get_file_list(flats_location, lambda x: x.split("_")[0] == "MasterFlat")
if len(flat_files) == 0:
    print("    No master flat found. Attempting to stack flats.")
    flat_files = get_file_list(flats_location)
    if len(flat_files) == 0:
        print("    Warning: no flats found")
        master_flat = np.ones(dimensions)
    else:
        print("    Stacking flats")
        master_flat = simplestack(flat_files, mode="premedian", dark=master_dark, bias=master_bias)
        master_flat /= np.percentile(master_flat,95)
        hdu = fits.PrimaryHDU(master_flat)
        hdu.writeto(flats_location[0] + "/MasterFlat_" + conf['job_id'] + ".fits", overwrite=True)
else:
    print("Master flat found.")
    master_flat = simplestack(flat_files)

# Stack and save =============================================================
print(f"Reading and stacking from {images_location} ({len(light_files)} images)")
stacked = stack(conf, light_files, recenter_filter=recenter_filter, flat=master_flat, dark=master_dark, bias=master_bias)

print("Job complete. Writing outputs.")
#output = debayer(stacked['data'], debayer_matrix=colour).astype("float32")
output = stacked['data']
np.savetxt(images_location[0] + "/lightcurves_" + conf['job_id'] + ".csv", stacked['lightcurves'])

#print(np.min(stacked), np.min(dark))
#stacked -= np.percentile(stacked, 30)
#stacked -= np.min(stacked)-1


hdu = fits.PrimaryHDU(output)
hdu.writeto(images_location[0] + "output_" + conf['job_id'] + ".fits", overwrite=True)
# if this doesnt work you'll have to catch a FileNotFoundError
os.replace('stack_overflow.npy', images_location[0] + "stack_" + conf['job_id'] + ".npy")
print("Dynamic Range:", np.min(output), np.max(output))


if not conf['debayer_matrix']:
    hdu = fits.PrimaryHDU(output)
    hdu.writeto(images_location[0] + "/output_" + conf['job_id'] + ".fits", overwrite=True)
else:
    hdu = fits.PrimaryHDU(np.sum(output, axis=2))
    hdu.writeto(images_location[0] + "/output_grey_" + conf['job_id'] + ".fits", overwrite=True)
    hdu = fits.PrimaryHDU(output[:,:,0])
    hdu.writeto(images_location[0] + "/output_red_" + conf['job_id'] + ".fits", overwrite=True)
    hdu = fits.PrimaryHDU(output[:, :, 1])
    hdu.writeto(images_location[0] + "/output_green_" + conf['job_id'] + ".fits", overwrite=True)
    hdu = fits.PrimaryHDU(output[:, :, 2])
    hdu.writeto(images_location[0] + "/output_blue_" + conf['job_id'] + ".fits", overwrite=True)