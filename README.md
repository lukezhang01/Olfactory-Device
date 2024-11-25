# Experiment 12 Flourescene Levels 
Outline of experiment can be found here: https://docs.google.com/document/d/1BuckM3dgzxQlzf2iB9w16xnlOVBJWG-d9uebHRNtUQo/edit#heading=h.mf4r9x54241h



# Signal-to-Noise Ratio and Flux Detection Script

## Overview
This project implements a procedure to detect and analyze objects (e.g., cells) in images based on their flux and signal-to-noise ratio (SNR). The script processes image data, calculates flux for regions of interest, and evaluates the SNR for better object characterization.

### Camera Models:
Small Ximea 

<img width="351" alt="image" src="https://github.com/user-attachments/assets/c6b52fec-d618-4495-aa50-8ae680097ddd">

Large Ximea

<img width="355" alt="image" src="https://github.com/user-attachments/assets/43585229-512a-4615-b1e3-e1a710c90427">


Notably the large Ximea camera has a higher average pixel value.
Manually using DS9, I noticed a larger gap between the background and fluorescent pixel values. (*) That is, 
For the large Ximea, the fluorescence pixels are usually between 1000-1500 value while the background pixels are usually between 100-200 value 
For the small Ximea, the fluorescence pixels are usually between 40-60 value while the background pixels are usually between 30-50 value

---

## Requirements
- **Dependencies:**
  - Python 3.x
  - NumPy, Matplotlib, DAOStarFinder, OpenCV
- **Environment:**
  - Works best with dark-field images or data with minimal noise.

---

## Features
- **Image Subsetting:** Extracts sub-images (regions of interest) from a larger dataset for detailed analysis.
- **Flux Calculation:** Determines the total flux from objects using predefined apertures and statistical measures.
- **Signal-to-Noise Ratio Calculation:** Computes SNR using physical parameters like exposure time, read noise, and pixel counts.
- **Visualization:** Displays sub-images with adjustable visualization parameters for better understanding.

---

## Key Functions
### `sub_image(x_coord, y_coord, data, boxsize)`
  1. Extracts a subsection of the image centered at (x_coord, y_coord) with dimensions defined by boxsize.
    For each pixel (x, y) that is assumed to be the center of a light source, create a subimage centred around (x, y). The subimage is 2*boxsize in width and height, where boxsize is a hyperparameter. See example below:

<img width="450" alt="image" src="https://github.com/user-attachments/assets/31a643e1-11aa-47e6-a490-946be800918e">

### `find_exact_center(sub_im)`
Determines the precise center of an object within a sub-image using intensity-weighted statistics:

1. Calculating Weighted Centroids:
   - The center coordinates $(x_{\text{center}}, y_{\text{center}})$ are computed as:
   
$$
x_{\text{center}} = \frac{\sum_{x, y} x \cdot I(x, y)}{\sum_{x, y} I(x, y)}
$$
$$
y_{\text{center}} = \frac{\sum_{x, y} y \cdot I(x, y)}{\sum_{x, y} I(x, y)}
$$

2. Subpixel Precision:
   - Refines the center using Gaussian fitting for higher accuracy.

---

### `flux_calc(sub_im, std_rad)`
Calculates the flux (total light) from the object using an aperture-based method:

1. Total Flux Calculation:
   - Defines a circular aperture of radius `std_rad` around the object's center.
   - The flux is the sum of pixel intensities within the aperture:
$$
\text{Flux} = \sum_{x, y \in \text{aperture}} I(x, y)
$$

4. Background Subtraction:
   - Corrects the flux by subtracting the background intensity:
 $$
 \text{Corrected Flux} = \text{Flux} - \text{Background}
 $$
   - The background intensity is estimated using the mean or median of pixels outside the aperture.

5. Flux Uncertainty:
   - Computes uncertainty using the standard deviation of pixel intensities:
 $$
 \sigma_{\text{flux}} = \sqrt{\sum_{x, y \in \text{aperture}} (\sigma_{\text{intensity}}^2)}
 $$
<img width="602" alt="image" src="https://github.com/user-attachments/assets/bbc29820-deb0-43a4-a10d-9dafba8da131">

6. Signal Clipping:
   - Removes extreme outliers based on statistical thresholds.

---

### `signal_to_noise(N_star, t, p, R, S_sky)`
Calculates the signal-to-noise ratio using the formula:
$$
\text{SNR} = \frac{N_{\text{star}} \cdot t}{\sqrt{N_{\text{star}} \cdot t + p \cdot (S_{\text{sky}} + R^2)}}
$$
Where:
- $N_{\text{star}}$: Count rate from the object (e-/s).
- $t$: Exposure time (s).
- $p$: Number of pixels in the aperture.
- $R$: Read noise (default = 100). This is based on the specifications of the large Ximea camera. Unable to find specifications for small Ximea camera. 
- $S_{\text{sky}}$: Background sky count (default = 0). This was chosen because the images were taken in a dark room

---

## Workflow
1. Load Image Data:
   - Load the (FITS) image file into a matrix format for processing.

2. Identify Coordinates of Interest:
   - Provide initial coordinates for objects in the image.

3. Extract Sub-Images:
   - Use `sub_image` to extract square regions around each object's center.

4. Analyze Sub-Images:
   - Find exact centers with `find_exact_center`.
   - Compute flux using `flux_calc`.
   - Calculate SNR using physical parameters and `signal_to_noise`.

5. Output Results:
   - Display and log flux, flux uncertainty, and SNR for each detected object.

---

## Inputs
- **`data`**: The image matrix containing raw pixel data.
- **`coords`**: List of tuples specifying initial positions of objects in the image.
- **SNR Parameters:**
  - `t`: Exposure time (default = 0.2 seconds).
  - `R`: Read noise (default = 100 e-).
  - `S_sky`: Sky background level (default = 10).

---

## Outputs
- A summary for each detected object including:
  - **Flux**: Total light collected from the object.
  - **Flux Uncertainty**: Statistical uncertainty in the flux calculation.
  - **SNR**: Ratio indicating the quality of the detected signal.

---

## Week 12 Report
### SNR with noise removal
For noise removal I used a combination of OpenCV histogram equalization and low pass filtering. This did not remove all noise, but only few pixels needed to be removed and a few cells needed to be added, so I did the remaining adjustments manually. 

<img width="457" alt="image" src="https://github.com/user-attachments/assets/f0a7eb44-f672-47ce-bb5d-6a60e3bb4a42">

Cell 1:
Flux = 29199.86666666667, Flux uncertainty = 156.67286811852827
Signal to noise ratio = 186.3747502507964

Cell 2:
Flux = 32518.43111111111, Flux uncertainty = 148.35183854569823
Signal to noise ratio = 219.1980323930677

Cell 3:
Flux = 7397.760000000006, Flux uncertainty = 112.34563814265096
Signal to noise ratio = 65.84821736120004

Cell 4:
Flux = 25424.87555555555, Flux uncertainty = 135.60028071258557
Signal to noise ratio = 187.49869411734758

Cell 5:
Flux = 33992.01333333333, Flux uncertainty = 154.40757125186306
Signal to noise ratio = 220.1447316199735

Cell 6:
Flux = 27787.675555555554, Flux uncertainty = 148.8716641253053
Signal to noise ratio = 186.65523569459575

<img width="472" alt="image" src="https://github.com/user-attachments/assets/5ae32661-8c9e-4749-9b48-200f1ebdc771">

Cell 1: 
Flux = 304905.3008220226, Flux uncertainty = 945.4518237735476, Signal to noise ratio = 322.4969196262852

Cell 2: 
Flux = 185420.08442689497, Flux uncertainty = 841.0128082173887, Signal to noise ratio = 220.47236690712416

Cell 3: 
Flux = 59509.27888799353, Flux uncertainty = 758.4066533267796, Signal to noise ratio = 78.46618779905716

Cell 4: 
Flux = 198479.94834810637, Flux uncertainty = 965.414499195535, Signal to noise ratio = 205.59039512405985

Cell 5: 
Flux = 174602.7942163934, Flux uncertainty = 855.273879414225, Signal to noise ratio = 204.14840020132315

Cell 6: 
Flux = 181020.04047276466, Flux uncertainty = 785.0437422092123, Signal to noise ratio = 230.58592883417095

**Attempted Autodetection**
  1. From last week, both cameras had noise possibly in the lens or in the camera itself. The noise is about 1-2px wide and carries an extremely high pixel value, which affected the detection of DAOStarfinder
    For small Ximea:

<img width="635" alt="image" src="https://github.com/user-attachments/assets/3b045056-3084-4665-b362-8b2d4f339c1e">

Same result for large Ximea:

<img width="364" alt="image" src="https://github.com/user-attachments/assets/4d76ecd2-e288-4d9c-9806-40c12b370412">
    
    I tried to adjust parameters such as threshold or reducing clustering (see the top left light source in the large Ximea has many overlap detection) but it still only detects noise rather than the real light sources from cells.
    I tired to add a maximum intensity value as well as parameters such as sharplo=0.2, sharphi=1.0 for DAOStarfinder which should filter out huge spikes in pixel values. This seemed to have some effect on the large Ximea camera, however, the noise appears to vary in pixel value, so it is hard to filter them out without filtering out the real light sources from cells.
Reducing max intensity:

<img width="298" alt="image" src="https://github.com/user-attachments/assets/4437ca97-3051-4dc4-961d-cd6862fcda9a">

Reducing it even more:

<img width="291" alt="image" src="https://github.com/user-attachments/assets/b0db9530-a122-42c6-b009-557c21382400">

As for the small Ximea camera, the DAOStarfinder cannot identify any real cell fluorescence. At the lowest possible max intensity cut off, it still only finds noise:

<img width="296" alt="image" src="https://github.com/user-attachments/assets/c441ee50-a170-427b-bde0-d8f91096368d">

If I set the threshold to be lower, then nothing is captured:

<img width="295" alt="image" src="https://github.com/user-attachments/assets/64cfdd60-5861-42cf-8669-094349c18b6f">
    
  2. Last week, with help from Adi, we managed to debug all syntax errors that prevented the code from running. 
  Although the code ran, it contained many logical errors. I outlined the brief fixes that may have resulted in slightly different calculations below:
  Recall calculations for ‘std’ value. I added edge cases that account for when the light source is very faint, resulting in a very flat Gaussian. It usually returns a very large std, but in the case it exceeds the boxsize/2, then std = boxsize/2. 
  Additionally, if the pixel values are too similar, which may happen if the aperture is very small ends up only capturing the white part, then autocalculating std with numpy may result in NaN. So if this happens, then I set std to 0.01.


