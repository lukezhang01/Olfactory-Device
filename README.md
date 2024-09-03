
Experiment 12 Flourescene Levels 
Compare captured light levels with the same Min/Max Scale:
Small Ximea 

<img width="351" alt="image" src="https://github.com/user-attachments/assets/c6b52fec-d618-4495-aa50-8ae680097ddd">

Large Ximea

<img width="355" alt="image" src="https://github.com/user-attachments/assets/43585229-512a-4615-b1e3-e1a710c90427">


Notably the large Ximea camera has a higher average pixel value.
Manually using DS9, I noticed a larger gap between the background and fluorescent pixel values. (*) That is, 
For the large Ximea, the fluorescence pixels are usually between 1000-1500 value while the background pixels are usually between 100-200 value 
For the small Ximea, the fluorescence pixels are usually between 40-60 value while the background pixels are usually between 30-50 value

SNR calculations
SNR Calculation Procedure:
For each pixel (x, y) that is assumed to be the center of a light source, create a subimage centred around (x, y). The subimage is 2*boxsize in width and height, where boxsize is a hyperparameter. See example below:

<img width="450" alt="image" src="https://github.com/user-attachments/assets/31a643e1-11aa-47e6-a490-946be800918e">

The std value is found by fitting all the intensity levels of each pixel in the subimage to a Gaussian distribution and taking the std of the distribution. 
Draw aperture
For each pixel (x, y), if the distance of (x, y) from the center is <= (1.5*std)**2, then it will be considered inside the aperture. Since (1.5*std)**2 is fixed for every image, then this will return a circle however, different images will have a different-sized aperture, depending on how much the light source deviates.
Draw annulus:
For each pixel (x, y), if the distance of (x, y) from the center is >= (2.5*std)**2 and <= (4*std)**2 then it will be considered inside the annulus. This will return a torus shape.
Flux calculation:
Sum of : (intensity levels in the aperture - average intensity level in annulus)
Uncertainty Calculation:

Finally SNR=(Flux/Flux uncertainty)
Fixes and updates to code:
Last week, with help from Adi, we managed to debug all syntax errors that prevented the code from running. 
Although the code ran, it contained many logical errors. I outlined the brief fixes that may have resulted in slightly different calculations below:
Recall calculations for ‘std’ value. I added edge cases that account for when the light source is very faint, resulting in a very flat Gaussian. It usually returns a very large std, but in the case it exceeds the boxsize/2, then std = boxsize/2. 
Additionally, if the pixel values are too similar, which may happen if the aperture is very small ends up only capturing the white part, then autocalculating std with numpy may result in NaN. So if this happens, then I set std to 0.01.

Attempted Autodetection
From last week, both cameras had noise possibly in the lens or in the camera itself. The noise is about 1-2px wide and carries an extremely high pixel value, which affected the detection of DAOStarfinder
For small Ximea

Same result for large Ximea

I tried to adjust parameters such as threshold or reducing clustering (see the top left light source in the large Ximea has many overlap detection) but it still only detects noise rather than the real light sources from cells.
I tired to add a maximum intensity value as well as parameters such as sharplo=0.2, sharphi=1.0 for DAOStarfinder which should filter out huge spikes in pixel values. This seemed to have some effect on the large Ximea camera, however, the noise appears to vary in pixel value, so it is hard to filter them out without filtering out the real light sources from cells.
Reducing max intensity:

Reducing it even more:

As for the small Ximea camera, the DAOStarfinder cannot identify any real cell fluorescence. At the lowest possible max intensity cut off, it still only finds noise:

If I set the threshold to be lower, then nothing is captured:

Manual Detection


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



Cell 1:
Flux = 3049053.0082202256, Flux uncertainty = 945.4518237735476
Signal to noise ratio = 3224.969196262852

Cell 2:
Flux = 1854200.8442689497, Flux uncertainty = 841.0128082173887
Signal to noise ratio = 2204.7236690712416

Cell 3:
Flux = 595092.7888799353, Flux uncertainty = 758.4066533267796
Signal to noise ratio = 784.6618779905716

Cell 4:
Flux = 1984799.4834810637, Flux uncertainty = 965.414499195535
Signal to noise ratio = 2055.9039512405984

Cell 5:
Flux = 1746027.942163934, Flux uncertainty = 855.273879414225
Signal to noise ratio = 2041.4840020132315

Cell 6:
Flux = 1810200.4047276466, Flux uncertainty = 785.0437422092123
Signal to noise ratio = 2305.8592883417095


Possible Next Steps
Improving autodetection so it is able to filter out noise more effectively (for the large Xiema only)
I plan to ask Adi for further guidance in using the photutils and DAOStarfinder libraries
Are there additional plots I should generate?
Rei told me that his contract ended at the end of July. Are there any tasks I should do with the actual optical device without Rei?

Confirm Flux, signal to noise ratio calculation
Get autodetection working

Number of pixels inside aperture 
Show mask of the aperture 
Flux of the pixels in aperture
mean/sd for the annulus - take all of the counts from the annulus
Mean - average intensity of the annulus

SNR 
S = Flux - number of pixels in aperture * background
N = sqrt(flux + number of pixels in aperture * sd_)

Intensity vs flux
Spatial distribution - show aperture mask
Redo SNR
