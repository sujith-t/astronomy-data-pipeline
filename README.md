# <u>Pipeline to Estimate Metallicity and Element Abundance</u>

The purpose of the pipeline is to extract important spectral lines from galaxies to determine element abundance against hydrogen and also to estimate metallicity profile by identifying emission and absorption lines from galaxies

Upon execution it queries the SDSS DR19 and fetches plate, mjd, fiber id so that it can identify relevant FITS file, download locally and continue processing. Following are the important activities taking place:
1. Galaxy information is fetched from a local database in a batch with the user specified batch-size and offset in a commandline parameter
2. Using the RA, DEC stored in local DB; plate, mjd, fiber id are fetched from SDSS DR 19 query
3. Using the above 3, spectra FITS file is downloaded and the above 3 parameters are updated in the local DB as a point of reference if it required to download again
4. Downloaded FITS galaxy spectra is extracted and updated in the local database

### Spectral Lines Under Consideration
- H Alpha 6564.61, *Hα star formation*
- H Beta 4862.68
- O3 5007, *Metallicity*
- O3 4959
- O2 3727: 3727.09
- N2 6583: 6583, *Metallicity*
- N2 6548
- S2 6716, *Density*
- S2 6730, *Density*
- S3 9069
- S3 9532
- Ne 3868, *Ionization*
- He 4685, *Ho star/AGN*
- Fe 5200, *Stellar Population*
- O3 4363

### Reusable Library Functions
**SpectralProfiler** class has reusable functions for following purposes:

- *plot_line_wavelength_vs_flux()* Extracts data from a FITS file, shifts it into the galaxy's "rest-frame" to remove the expansion of the universe (redshift correction), cuts out a specific wavelength window, and plots a spectrum graph. This is useful to visually inspect specific emission lines to check data quality before doing automated calculations
- *__fit_line (Private)* Isolates a tiny region of data centered around an emission line and mathematically fits a Gaussian profile (a bell-shaped curve) to it using least-squares optimization. It then calculates the area under that curve. The total area under the curve equals the total raw flux (intensity) of that specific atomic emission line

S. Thillayampalam, "Morphological Characteristics of Galaxies in Studying Galaxy Evolution via Optimising Classification," 2025 IEEE International Conference on Future Machine Learning and Data Science (FMLDS), Los Angeles, CA, USA, 2025, pp. 460-466, doi: 10.1109/FMLDS67896.2025.00088.


