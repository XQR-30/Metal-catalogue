# README

This folder contains the data files accompanying the XQR-30 metal absorber catalog paper (Davies et al. 2023a). The files are arranged in four sub-folders.

## AbsorberCatalogs
This folder includes:

* The full metal absorber catalog (provided as a CSV file).
* Sub-folders for each individual quasar containing:
    * A catalog of metal absorbers in the spectrum of that quasar (provided as a CSV file).
    * Voigt profile fits for each individual absorption system (provided as a PDF file).

## AbsorptionPathTool
This folder contains the following files:

* ```absorption_path_tool.py```: Python3 tool to calculate the absorption path covered by any subset of the E-XQR-30 quasars for any of the primary ions (MgII, FeII, CII, CIV, SiIV, NV) or for any list of (one or more) transitions. The script uses the Planck 2018 cosmology by default but the Omega_m value can be changed (assumes flat cosmology).
* ```interface.py```: User interface for the ```absorption_path_tool.py``` script with customizable inputs. Open and read this file for more information.
* ```requirements.txt```: lists the package dependencies for the Python script (AstroPy 5.1 and NumPy 1.22.4). The dependencies can be automatically installed by running  
```pip install -r requirements.txt```
* ```transition_wavelengths.dat```: rest-frame wavelengths of common absorption lines. Used by ```absorption_path_tool.py```.

## ContinuumFits
This folder contains one PDF file per quasar showing the flux (beginning at the observed wavelength of the quasar Lya emission line), error spectrum, and the fitted continuum level.

## Tables
This folder contains machine-readable (CSV or JSON) versions of Tables B1, B2 and B4.

* ```redshifts_and_spectral_resolutions.csv```: Adopted quasar emission redshifts and measured spectral resolutions (presented in Table B1).
* ```masked_wavelength_regions.json```: Wavelength regions that were masked when performing the initial automatic search for absorption systems (Table B2).
* ```bal_wavelength_regions.json```: Wavelength regions in which Broad Absorption Line signatures were detected (originally published in Bischetti et al. 2022, Nature, 605, 244 and Bischetti et al. 2023, submitted) (presented in Table B4).

## XQR-30_data_movie.mp4
This is an animated version of Figure 5. The first frame shows only absorbers detected in the lowest redshift quasar, and each subsequent frame adds data for one additional quasar, sorted in ascending order of redshift.