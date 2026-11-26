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
- *__find_rest_wavelength (Private)* Opens a FITS file, reads the observed wavelengths and fluxes, looks up the galaxy's redshift (z), and applies the formula λ<sub>rest</sub> = λ<sub>observed</sub> / (1 + z) to convert wavelengths to the galaxy's static frame of reference. It maps out exactly where atomic emission lines should be, eliminating the shifting effect caused by the galaxy moving away
- *__detect_emission_flux (Private)* Performs continuum subtraction. It finds areas of the spectrum with no emission lines (continuum windows), fits a 3rd-degree smooth polynomial curve through them to model the baseline starlight, and subtracts this baseline model from the overall spectrum. Finally, it triggers __fit_line on the remaining pure emission peaks. Gas emission line spikes sit on top of background starlight from the galaxy. To accurately measure the gas, need to peel away the starlight first
- *star_formation_rate()* Calculates how quickly the galaxy is birthing new stars (in solar masses per year, M⊙/yr. It utilizes the ratio of Hydrogen-Alpha (Hα) to Hydrogen-Beta (Hβ) lines to determine how much dust is blocking the light, corrects the Hα flux, calculates true luminosity based on distance, and applies standard Kennicutt (1998) scaling. It provides an index of active star birth explains how fast a galaxy is growing
- *__dust_bias_correction (Private)* Observed Hα / Hβ ratio (the Balmer Decrement) and compares it against the theoretical cosmic baseline of 2.86. If the observed ratio is higher, it infers the presence of dust and calculates an extinction correction multiplier for all lines. Dust absorbs blue light more heavily than red light. Removes the dust effect from the spectrum and provides true original light emitted by the gas cloud
- *__detect_iron_flux (Private)* Measures the average flux across two specific wavelength absorption windows for Iron (Fe<sub>5270</sub> and Fe<sub>5335</sub>). Unlike emission lines (which show hot gas), iron lines are absorption features identifies the older, solid stellar populations making up the bulk of the galaxy
- *__metallicity_empirical (Private)* Implements fallback equations (Pettini & Pagel 2004; Kobulnicky & Kewley 2004) to approximate oxygen abundance using strong lines like Oxygen (O<sub>2</sub), (O<sub>3</sub), Nitrogen (N<sub>2</sub), and Hydrogen (Hα, Hβ) when weak, direct lines can't be detected. Deep galaxies are often too faint for direct methods, making these statistical approximations important
- *element_abundance_profile()* The main engine in calculating the metallicity. It checks the signal-to-noise ratio of a rare temperature-sensitive line (O<sub>3</sub> λ<sub>4363</sub>). If it is clean, it runs the advanced, hyper-accurate "Direct Temperature Exact" physical method. If it's too noisy, it relies on empirical proxies (O3N2 and R23 calibrations). It translates these into relative abundances for Oxygen, Nitrogen, Carbon, Sulphur, Neon, and Iron relative to Hydrogen
- *__direct_metallicity_exact (Private)* Utilizes PyNeb to run an iterative numerical solver across a 5-level quantum atom model. It solves the physical feedback loop between gas electron density (*n<sub>e</sub>*) and electron temperature (*T<sub>e</sub>*) using physical atomic ratios (S<sub>2</sub>) and (O<sub>3</sub>) lines to yield an exact, non-empirical chemical profile. Gold standard of chemical calculation in nebular, providing accurate physical measurements instead of statistical guesses
- *corrected_emission_flux()* It opens a target FITS file, reads the pre-calculated emission line catalog (SPZLINE), extracts all target atomic elements, fills in missing features (like Iron or uncalculated weak lines) by scanning the raw spectrum data via private sub-methods, and returns a fully corrected, organized baseline dataset ready for profiling. Simply cleaning unstructured raw file inputs into pristine datasets

### References and Initial Dataset
**S. Thillayampalam**, *"Morphological Characteristics of Galaxies in Studying Galaxy Evolution via Optimising Classification,"* 2025 IEEE International Conference on Future Machine Learning and Data Science (FMLDS), Los Angeles, CA, USA, 2025, pp. 460-466, doi: 10.1109/FMLDS67896.2025.00088.


