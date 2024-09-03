from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from photutils.detection import DAOStarFinder
from scipy.optimize import curve_fit

plt.rcParams.update({'font.size': 12})


def open_fits_file(file_path):
    with fits.open(file_path, do_not_scale_image_data=True) as hdul:
        primary_hdu = hdul[0]
        image_data = primary_hdu.data

        bscale = primary_hdu.header.get('BSCALE', 1.0)
        bzero = primary_hdu.header.get('BZERO', 0.0)
        image_data = bscale * image_data + bzero
    return image_data


def show_data(data, vmin: int = 0.0001, vmax: int = 1000):
    plt.figure()
    plt.imshow(data, cmap='inferno', aspect='auto', norm=LogNorm(vmin=vmin, vmax=vmax), origin='lower')
    cb = plt.colorbar()
    cb.set_label(label='Intensity [ADU]', size=20)


def cell_finder(data, threshold: int = 110, fwhm: int = 5, x_min: int = 0, x_max: int = 9576, y_min: int = 0,
                y_max: int = 6388):
    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold)
    cells = daofind(data[x_min:x_max, y_min:y_max])

    coords = []
    for i in range(len(cells['xcentroid'])):
        coords.append([int(cells['xcentroid'][i] + x_min), int(cells['ycentroid'][i] + y_min)])

    print(f'Number of cells detected = {len(coords)}')

    return coords


def show_cells(data, coords, vmin: int = 0.0001, vmax: int = 1000):
    x_coords = [i[0] for i in coords]
    y_coords = [i[1] for i in coords]

    plt.figure()
    plt.imshow(data, cmap='inferno', aspect='auto', norm=LogNorm(vmin=vmin, vmax=vmax), origin='lower')
    cb = plt.colorbar()
    cb.set_label(label='Intensity [ADU]', size=20)
    plt.scatter(x_coords, y_coords, s=50, color='b', facecolor='none')


def sub_image(x_coord, y_coord, data, boxsize: int = 100):
    sub_im = data[x_coord - boxsize: x_coord + boxsize + 1, y_coord - boxsize: y_coord + boxsize + 1]
    return sub_im


def find_exact_center(sub_img):
    stack_x = np.zeros(sub_img.shape[1])
    stack_y = np.zeros(sub_img.shape[0])
    for j in range(sub_img.shape[0])[-10 + sub_img.shape[0] // 2:10 + sub_img.shape[0] // 2]:
        stack_x += sub_img[j, :]
        stack_y += sub_img[:, j]

    x_cent = np.where(stack_x == max(stack_x))[0][0]
    y_cent = np.where(stack_y == max(stack_y))[0][0]
    center_coords = np.array([x_cent, y_cent])

    return [stack_x, stack_y], center_coords


def std_radius(stacks, center_coords, boxsize: int = 100):
    def gaus(x, a, x0, sigma):
        return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

    axis = np.arange(0, 2 * boxsize + 1)

    try:
        popt_x, _ = curve_fit(gaus, axis, stacks[0], p0=[max(stacks[0]), center_coords[0], 10], maxfev=2000)
        popt_y, _ = curve_fit(gaus, axis, stacks[1], p0=[max(stacks[1]), center_coords[1], 10], maxfev=2000)
    except RuntimeError as e:
        print(f"Error fitting Gaussian: {e}")
        return None

    std_rad = np.sqrt((popt_x[2] ** 2) + (popt_y[2] ** 2))

    return np.round(std_rad)


def aperture(sub_image, std: float):
    i, j = np.indices(sub_image.shape)
    radius = (i - len(sub_image[0]) / 2) ** 2 + (j - len(sub_image[1]) / 2) ** 2
    mask = radius <= (1.5 * std) ** 2
    masked = np.zeros_like(sub_image)
    masked[mask] = 1

    return masked


def annulus(sub_image, std: float):
    i, j = np.indices(sub_image.shape)
    radius = (i - len(sub_image[0]) / 2) ** 2 + (j - len(sub_image[1]) / 2) ** 2
    mask = (radius >= (2.5 * std) ** 2) & (radius <= (4 * std) ** 2)
    masked = np.zeros_like(sub_image)
    masked[mask] = 1

    return masked


def flux_calc(centered_sub_img, std_rad: float):
    N = np.count_nonzero(annulus(centered_sub_img, std_rad))
    bg = np.sum(centered_sub_img * annulus(centered_sub_img, std_rad)) / N if N != 0 else 0
    flux = np.sum((centered_sub_img - bg) * aperture(centered_sub_img, std_rad))

    return flux


def flux_uncertainty(sub_im_cent, std: float, flux: float):
    app = sub_im_cent * aperture(sub_im_cent, std)
    ann = sub_im_cent * annulus(sub_im_cent, std)
    aperture_pixels = np.count_nonzero(app)
    ann_1d = ann[ann != 0]
    ann_std = np.std(ann_1d)
    bg_unc = np.sqrt(aperture_pixels * ann_std)
    source_unc = np.sqrt(flux / 4)
    flux_unc = np.sqrt(bg_unc ** 2 + source_unc ** 2)

    return flux_unc


def signal_to_noise(flux: float, flux_uncertainty: float):
    return flux / flux_uncertainty


def view_aperture_annulus(sub_im_cent, std_rad: float, vmin: int = 0.0001, vmax: int = 1000):
    plt.figure()
    plt.imshow(sub_im_cent * aperture(sub_im_cent, std_rad) + sub_im_cent * annulus(sub_im_cent, std_rad),
               cmap='inferno', aspect='auto', norm=LogNorm(vmin=vmin, vmax=vmax), origin='lower')
    plt.title('Masked Target Star')
    cb = plt.colorbar()
    cb.set_label(label='Intensity [ADU]', size=20)


if __name__ == '__main__':
    # Define image path
    image_path = '../images/experiment 9/515LP_LANO/with_signal_1.fits'

    # Load in the data
    data = open_fits_file(image_path)
    print(data)

    # View image if needed
    show_data(data, vmin=0.0001, vmax=1000)

    # Find the approximate positions of the cells
    coords = cell_finder(data, threshold=110, fwhm=5, x_min=0, x_max=9576, y_min=0, y_max=6388)

    # View cell locations to confirm the code works as expected
    show_cells(data, coords, vmin=0.0001, vmax=1000)

    # Initialize empty lists to store values from each cell
    fluxes = []
    uncertainties = []
    s2n_ratios = []

    # This essentially determines the size of the sub-image, feel free to change if needed
    boxsize = 100

    # Start a for-loop to do the analysis for each cell individually
    for i in coords:
        x_coord, y_coord = i

        # Make the sub-image centered at the approximate center coords found by the algorithm
        sub_im = sub_image(x_coord, y_coord, data, boxsize)

        # Find the actual coordinates of the center within the sub-image
        stacks, center_coords = find_exact_center(sub_im)

        # Recalculate the coords of the center to make them be at the center of the sub-image
        cent_x_coord, cent_y_coord = center_coords - boxsize + i

        # Remake the sub-image with the new center coords
        cent_sub_img = sub_image(cent_x_coord, cent_y_coord, data, boxsize)

        # Calculate the radius factor needed for the aperture and annulus
        std_rad = std_radius(stacks, center_coords, boxsize)
        if std_rad is None:
            continue

        # View the aperture and annulus if needed
        view_aperture_annulus(cent_sub_img, std_rad, vmin=0.0001, vmax=1000)

        # Calculate the background subtracted flux of the cell
        flux = flux_calc(cent_sub_img, std_rad)

        # Calculate the uncertainty in the flux
        flux_unc = flux_uncertainty(cent_sub_img, std_rad, flux)

        # Calculate the signal-to-noise ratio of the cell
        s2n = signal_to_noise(flux, flux_unc)

        # Append all the calculated values to their respected lists
        fluxes.append(flux)
        uncertainties.append(flux_unc)
        s2n_ratios.append(s2n)
