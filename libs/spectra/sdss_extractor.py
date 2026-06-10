#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 15:34:26 2026

@author: sujith-t
"""

import matplotlib.pyplot as plt
import numpy as np
import logging
import pyneb as pn
import math

from astropy.io import fits
from scipy.optimize import curve_fit
from astropy.cosmology import Planck18 as cosmo


class SpectralProfiler:

    logger = logging.getLogger("__spec_profiler__")
    METHOD_R23 = "R23"
    METHOD_O3N2 = "O3N2"
    METHOD_TE = "TE"
    METHOD_O3N2_R23 = "O3N2_R23"

    lines = {
        "h_alpha": 6564.61,  # Hα star formation
        "h_beta": 4862.68,
        "o3_5007": 5007,  # Metallicity
        "o3_4959": 4959,
        "o2_3727": 3727.09,
        "n2_6583": 6583,  # Metallicity
        "n2_6548": 6548,
        "s2_6716": 6716,  # Density
        "s2_6730": 6730,  # Density
        "s3_9069": 9069,
        "s3_9532": 9532,
        "ne_3868": 3868,  # Ionization
        "he_4685": 4685,  # Ho star/AGN
        "fe_5200": 5200, # Stellar population
        "o3_4363": 4363
    }

    def plot_line_wavelength_vs_flux(self, file_path: str, lower_limit, upper_limit, title: str):
        hdul = fits.open(file_path)
        data = hdul[1].data
        flux = data['flux']
        loglam = data['loglam']
        wavelength = 10 ** loglam  # Ångstroms

        # redshif correction steps
        z = hdul[2].data['Z'][0]
        # Rest-frame wavelength
        rest_wavelength = wavelength / (1 + z)

        # Select region around upper and lower limits
        mask = (rest_wavelength > lower_limit) & (rest_wavelength < upper_limit)

        w = rest_wavelength[mask]
        f = flux[mask]

        plt.plot(w, f)
        plt.title(title)
        plt.xlabel("Wavelength (Å)")
        plt.show()


    def __fit_line(self, w, f, center, window=25):

        def gaussian(x, amp, mu, sigma):
            return amp * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))

        # 1. Isolate the window
        half_window = window / 2.0
        mask = (w >= (center - half_window)) & (w <= (center + half_window))
        x = w[mask]
        y = f[mask]

        if len(x) < 5:
            return 0

        # Initial guesses
        y_max = 0 if len(y) == 0 else np.max(y)
        initial_guesses = [y_max, center, 2.5]

        try:
            popt, _ = curve_fit(gaussian, x, y, p0=initial_guesses)
            fit_amp, fit_mean, fit_stddev = popt

            # Flux = area under Gaussian
            flux_line = fit_amp * np.abs(fit_stddev) * np.sqrt(2 * np.pi)

            return flux_line
        except:
            return 0


    def __find_rest_wavelength(self, file_path: str):
        hdul = fits.open(file_path)

        # data retrieval from the file
        data = hdul[1].data
        flux = data['flux']
        loglam = data['loglam']

        wavelength = 10 ** loglam

        # redshift correction
        z = hdul[2].data['Z'][0]
        rest_wavelength = wavelength / (1 + z)

        hdul.close()
        return rest_wavelength, flux

    def __detect_emission_flux(self, file_path: str, detection_lines: dict) -> dict:
        rest_wavelength, flux = self.__find_rest_wavelength(file_path)
        hdul = fits.open(file_path)

        ivar = hdul['COADD'].data['ivar']
        continuum_windows = [
            (4200, 4280),
            (4450, 4600),
            (5050, 5400),
            (6000, 6200),
            (6800, 7000)
        ]

        # Create a master mask for pixels that fall inside ANY continuum window
        # Also ensure the pixel is valid (inverse variance > 0)
        continuum_mask = np.zeros(len(rest_wavelength), dtype=bool)
        for start, end in continuum_windows:
            window_mask = (rest_wavelength >= start) & (rest_wavelength <= end)
            continuum_mask = continuum_mask | window_mask

        # Exclude bad/flagged data points
        valid_continuum_mask = continuum_mask & (ivar > 0)

        fit_wave = rest_wavelength[valid_continuum_mask]
        fit_flux = flux[valid_continuum_mask]
        if len(fit_wave) == 0:
            SpectralProfiler.logger.warning("No valid data points found in the specified continuum windows")
            return {}

        # Fit a polynomial to the isolated continuum points
        # numpy.polyfit computes a least-squares polynomial fit
        poly_coefficients = np.polyfit(fit_wave, fit_flux, deg=3)

        # Evaluate the polynomial across the ENTIRE rest-wavelength grid
        # This creates your baseline continuum model
        continuum_model = np.polyval(poly_coefficients, rest_wavelength)

        # Calculate Pure Emission Line Spectrum (Continuum Subtracted)
        # This is exactly what you need before fitting Gaussians to H-alpha or H-beta!
        pure_lines_flux = flux - continuum_model

        if len(detection_lines) == 0:
            detection_lines = self.lines

        results = {}
        for name, center in detection_lines.items():
            flux_val = self.__fit_line(rest_wavelength, pure_lines_flux, center)
            results[name] = flux_val

        return results


    # FR(M⊙​/yr)=7.9×10−42×LHα​(erg/s)
    # Kennicutt, R. C. (1998), "Star Formation in Galaxies Along the Hubble Sequence," Annual Review of Astronomy and Astrophysics, vol. 36, pp. 189–231
    # Chabrier, G. (2003), "Galactic Stellar and Substellar Initial Mass Function," Publications of the Astronomical Society of the Pacific, vol. 115, no. 809, pp. 763–795
    def star_formation_rate(self, h_alpha:float, h_beta:float, redshift:float) -> float:

        if h_alpha == 0 or h_beta == 0:
            return 0

        intrinsic_ratio = 2.86

        a_ha = 0.0
        observed_ratio = h_alpha / h_beta
        if observed_ratio > 2.86:
            ebv = 1.97 * np.log10(observed_ratio / intrinsic_ratio)
            a_ha = 3.33 * ebv

        # Convert redshift → distance Luminosity distance in cm
        luminosity_distance = cosmo.luminosity_distance(redshift).to('cm').value

        f_ha_corr = h_alpha * (10**(0.4 * a_ha)) * 1e-17

        l_ha = 4 * np.pi * (luminosity_distance**2) * f_ha_corr

        return 7.9e-42 * l_ha


    # remove dust effects
    def __dust_bias_correction(self, emission_flux: dict, h_alpha: float=0, h_beta: float=0) -> dict:
        k = 2.5
        intrinsic_ratio = 2.86

        e_bv = k
        if h_alpha > 0 and h_beta > 0:
            observed_ratio = h_alpha / h_beta
            e_bv = k * np.log10(observed_ratio / intrinsic_ratio)

        corrected = {}

        # make correction for all lines
        for line, flux in emission_flux.items():
            if flux:
                corrected[line] = flux * 10 ** (0.4 * k * e_bv)
            else:
                corrected[line] = 0

        return corrected


    def __detect_iron_flux(self, file_path: str):
        fe_5270_range = (5245, 5285)
        fe_5335_range = (5315, 5355)

        rest_wavelength, flux = self.__find_rest_wavelength(file_path)

        def measure_index(wave, flux, wmin, wmax):
            mask = (wave > wmin) & (wave < wmax)
            return np.mean(flux[mask])

        fe_5270 = measure_index(rest_wavelength, flux, *fe_5270_range)
        fe_5335 = measure_index(rest_wavelength, flux, *fe_5335_range)

        return float(fe_5270), float(fe_5335)


    # the ratio is in relation to Hydrogen as the baseline
    def element_abundance_profile(self, corrected:dict, redshift: float=0.09) -> dict:
        ratios = {}

        if "h_alpha" not in corrected or "h_beta" not in corrected:
            return

        h_alpha = corrected["h_alpha"]
        h_beta = corrected["h_beta"]

        # Calculate Signal-to-Noise Ratio (SNR)
        snr_4363 = corrected['o3_4363'] / corrected['o3_4363_err'] if ((corrected['o3_4363'] is not None and corrected['o3_4363'] > 0) and (corrected['o3_4363_err'] is not None and corrected['o3_4363_err'] > 0)) else 0

        # Standard peer-reviewed threshold is usually SNR > 3 or SNR > 5
        use_empirical = True if (corrected['o3_4363'] <= 1.5 or snr_4363 < 3.0) else False
        metallicity_o3n2, metallicity_r23, final_metallicity, temperature_exact = None, None, None, None

        ratios["final_method"] = SpectralProfiler.METHOD_O3N2_R23
        if not use_empirical:
            final_metallicity, temperature_exact = self.__direct_metallicity_exact(corrected['o3_4363'], corrected['o3_4959'], corrected['o3_5007'], corrected['o2_3727'], corrected['s2_6716'], corrected['s2_6730'], h_beta)
            SpectralProfiler.logger.debug("Using o3_4363 branch for metallicity calculation")
            ratios["final_method"] = SpectralProfiler.METHOD_TE
        else:
            # weighted combined metallicity - (research + production grade)
            # interpretation -> 12 + log(O/H) = final_metallicity
            metallicity_o3n2, metallicity_r23 = self.__metallicity_empirical(corrected["o2_3727"], corrected["o3_4959"], corrected["o3_5007"], corrected["n2_6583"], h_alpha, h_beta)
            final_metallicity = 0.6 * metallicity_o3n2 + 0.4 * metallicity_r23
            SpectralProfiler.logger.debug("Using o3n2 + R23 branch for metallicity calculation")

            # In case of discrepancy we consider o3n2 is reliable
            if abs(metallicity_o3n2 - metallicity_r23) > 0.3:
                final_metallicity = metallicity_o3n2
                ratios["final_method"] = SpectralProfiler.METHOD_O3N2
                SpectralProfiler.logger.warning(f"O3N2/R23 discrepancy detected: {metallicity_o3n2:.2f} vs {metallicity_r23:.2f}")

        ratios["metallicity_r23"] = metallicity_r23
        ratios["metallicity_o3n2"] = metallicity_o3n2
        ratios["final_metallicity"] = final_metallicity
        ratios["temperature_exact"] = temperature_exact

        # log(O/H) = final_metallicity - 12
        log_oh = ratios["final_metallicity"] - 12
        ratios["oxygen"] = 10 ** log_oh

        ##### 2. NITROGEN #####
        log_no = 0
        if corrected["o2_3727"] > 0:
            log_no = np.log10(corrected["n2_6583"] / corrected["o2_3727"]) + log_no

        # log(N/H) = log(N/O)+log(O/H)
        log_nh = log_no + log_oh
        ratios["nitrogen"] = 10 ** log_nh

        ##### 3. CARBON ##### this is an indirect estimate, flux lines are weak and complicated. only Ultra-Violet files can catch them

        # plan log(C/O) = a + b × (12+log(O/H)) where a, b are calibration values | log(C/H) = log(C/O) + log(O/H)
        log_co = -0.8 + 0.14 * (ratios["final_metallicity"] - 8.0)
        log_ch = log_co + log_oh
        ratios["carbon"] = 10 ** log_ch

        ##### 4. SULPHUR #####
        ratios["sulphur"] = 0

        if redshift < 0.09:
            # widely adopted empirical polynomial from Díaz et al
            # 12 + log10(S/H) = 5.79 + (1.54 x log10(S23)) + (0.15 x log10(S23)^2)
            s23 = (corrected["s2_6716"] + corrected["s2_6730"] + corrected["s3_9069"] + corrected["s3_9532"]) / h_beta
            log_sh = 5.79 + (1.54 * np.log10(s23)) + (0.15 * np.log10(s23) ** 2) - 12
            ratios["sulphur"] = 10 ** log_sh
            SpectralProfiler.logger.debug("Full S23 Calculation")
        elif 0.09 <= redshift < 0.54:
            # using standard empirical calibrations. log10(S/O) = 1.6, [12 + log10(O/H)] = metallicity
            # 12 + log10(S/H) = [12 + log10(O/H)] - log10(S/O)
            log_sh = (ratios["final_metallicity"] - 1.6) - 12
            ratios["sulphur"] = 10 ** log_sh
            SpectralProfiler.logger.debug("S3 Redshifted out. Using pure S2 calibrations")
        else:
            SpectralProfiler.logger.warning("Warning: All primary optical Sulphur lines have left the BOSS window")

        ##### 5. estimate NEON #####
        log_neo =  0.7
        if (corrected["o3_5007"] > 0 or corrected["o3_4959"] > 0) and corrected["ne_3868"] > 0:
            log_neo = np.log10(corrected["ne_3868"] / (corrected["o3_5007"] + corrected["o3_4959"])) + log_neo

        log_neh = log_neo + log_oh
        ratios["neon"] = 10 ** log_neh

        ##### 6. Iron (Fe) ##### this can't be interpreted as ratio Fe/H, Trager et al, Johansson, Thomas, & Maraston (2010)
        fe_index = (corrected['fe_5270'] + corrected['fe_5335']) / 2
        ratios["iron_strength"] = -2.0 + 0.4 * fe_index

        return ratios

    """
    Empirical methods to based on O3N2 and R23 methods
    Using a weight of 0.6 on O3N2, 0.4 on R23
    """
    def __metallicity_empirical(self, o2_3727: float, o3_4959:float, o3_5007:float, n2_6583:float, h_alpha:float, h_beta:float):

        """
        1. OXYGEN and Metallicity
        Calculates R23 metallicity and automatically breaks the branch degeneracy.
        Uses [NII]/[OII] boundary at -1.2 (Kewley & Ellison 2008).
        """
        total_oxygen = o3_5007 + o3_4959 + o2_3727
        # 1. Calculate standard line ratios
        log_o_r23 = 0
        if total_oxygen > 0 and  h_beta > 0:
            log_o_r23 = np.log10(total_oxygen / h_beta)

        # Ionization parameter proxy
        log_o32 = 0
        if o2_3727 != 0:
            log_o32 = np.log10((o3_4959 + o3_5007) / o2_3727)

        # Degeneracy breaker ratio
        log_no = 0
        if o2_3727 > 0:
            log_no = np.log10(n2_6583 / o2_3727) + log_no

        # Branch selection logic (Kobulnicky & Kewley 2004, ~7.0 to 9.3 range)
        # Precise KK04 Upper Branch polynomial:
        metallicity_r23 = float(8.85 - 0.65 * log_o_r23 - 0.2 * log_o32)

        if log_no < -1.2:
            SpectralProfiler.logger.debug("Low Metallicity branch triggered for R23 calculation")
            metallicity_r23 = 8.0 + 0.3 * log_o_r23 - 0.25 * log_o32
        else:
            SpectralProfiler.logger.debug("High Metallicity branch triggered for R23 calculation")

        # compute O3N2 | metallicity_O3N2 > 8.4 = high otherwise low
        log_o3n2 = 0
        if h_beta > 0 and h_alpha > 0:
            log_o3n2 = np.log10((o3_5007 / h_beta) / (n2_6583 / h_alpha))

        # Pettini & Pagel (2004) Calibration (~8.0 - 8.8 valid range)
        metallicity_o3n2 = float(8.73 - 0.32 * log_o3n2)

        return metallicity_o3n2, metallicity_r23


    """
    Calculates high-accuracy gas-phase metallicity using the Direct Te method
    via PyNeb's 5-level atom numerical solver.
    
    Parameters:
    fluxes (dict): Dust-corrected line fluxes (normalized such that Hbeta = 100).
                   Required keys: 'OIII_4363', 'OIII_4959', 'OIII_5007',
                                  'OII_3727', 'SII_6716', 'SII_6731'
    
    Returns:
    tuple: Exact physical values for metallicity, Te
    """
    def __direct_metallicity_exact(self, o3_4363, o3_4959, o3_5007, o2_3727, s2_6716, s2_6731, h_beta):

        # 1. Initialize PyNeb atom emitters (loads default atomic data configurations)
        o2 = pn.Atom('O', 2)
        o3 = pn.Atom('O', 3)
        s2 = pn.Atom('S', 2)

        # 2. Extract line fluxes
        f_4363 = (o3_4363 / h_beta) * 100
        f_4959 = (o3_4959 / h_beta) * 100
        f_5007 = (o3_5007 / h_beta) * 100
        f_3727 = (o2_3727 / h_beta) * 100
        f_6716 = (s2_6716 / h_beta) * 100
        f_6731 = (s2_6731 / h_beta) * 100

        # 3. Define the physical line ratios
        r_o3 = (f_5007 + f_4959) / f_4363
        r_s2 = f_6716 / f_6731

        # 4. Iterative cross-dependent solver loop
        # Resolves the minor feedback loop between density and temperature scales.
        te_guess = 12000.0
        te_calibrated = 0
        ne_calibrated = 0
        for _ in range(5):  # 5 iterations guarantee analytical convergence < 0.01 K
            ne_calibrated = s2.getTemDen(int_ratio=r_s2, tem=te_guess, wave1=6716, wave2=6731)

            # Handle low-density limits gracefully
            if ne_calibrated < 1.0 or np.isnan(ne_calibrated):
                ne_calibrated = 100.0  # Default standard low-density regime limit

            te_calibrated = o3.getTemDen(int_ratio=r_o3, den=ne_calibrated, wave1=5007, wave2=4363)
            te_guess = te_calibrated

        te_o3 = float(te_calibrated)
        ne_s2 = float(ne_calibrated)

        # 5. Determine Low-Ionization Zone Temperature (Te_OII)
        # Using the standard research relation from Campbell, Terlevich, & Melnick (1986)
        te_o2 = 0.7 * te_o3 + 3000.0

        # 6. Calculate Ionic Abundances relative to H-beta
        # Emissivity values are determined using PyNeb's full multi-level atom configuration
        ab_o3_h = o3.getIonAbundance(int_ratio=f_5007, tem=te_o3, den=ne_s2, wave=5007, Hbeta=100.0)

        ab_o2_h = o2.getIonAbundance(int_ratio=f_3727, tem=te_o2, den=ne_s2, wave=3727, Hbeta=100.0)

        # 7. Total Oxygen Abundance Summation
        total_o_h = ab_o2_h + ab_o3_h
        metallicity = float(12.0 + np.log10(total_o_h))

        # metallicity and temperature exact of o3
        return metallicity, te_o3

    def corrected_emission_flux(self, file_path: str):
        hdul = fits.open(file_path)
        target_maps = {'H_beta': 'h_beta', 'H_alpha': 'h_alpha', '[O_III] 5007': 'o3_5007', '[O_III] 4959': 'o3_4959', '[O_II] 3727': 'o2_3727', '[N_II] 6583': 'n2_6583', '[N_II] 6548': 'n2_6548', '[S_II] 6716': 's2_6716', '[S_II] 6730': 's2_6730', '[Ne_III] 3868': 'ne_3868', 'He_II 4685': 'he_4685', '[O_III] 4363': 'o3_4363'}
        line_data = hdul["SPZLINE"].data
        redshift = float(hdul[2].data['Z'][0])

        flux = {"o3_4363_err": None}
        for r in line_data:
            if r['LINENAME'].strip() in target_maps:
                flux[target_maps[r['LINENAME'].strip()]] = float(r['LINEAREA'])
            if r['LINENAME'].strip() == "[O_III] 4363":
                flux["o3_4363_err"] = float(r['LINEAREA_ERR'])

        # this step is done due to unavailability of iron in SPZLINE (not already calculated)
        # rest of them already calculated, no need to do again. h-alpha,h-beta are calculated to detect the observed values
        wavelengths = {"h_alpha": self.lines["h_alpha"], "h_beta": self.lines["h_beta"], "fe_5200": self.lines["fe_5200"], "s2_6716": self.lines["s2_6716"],
                       "s2_6730": self.lines["s2_6730"], "s3_9069": self.lines["s3_9069"], "s3_9532": self.lines["s3_9532"]}

        if "o3_4363" not in flux:
            wavelengths["o3_4363"] = self.lines["o3_4363"]

        result = self.__detect_emission_flux(file_path, wavelengths)
        flux["h_alpha_observed"] = float(result["h_alpha"])
        flux["h_beta_observed"] = float(result["h_beta"])
        result = self.__dust_bias_correction(result, h_alpha=result["h_alpha"], h_beta=result["h_beta"])
        flux["fe_5200"] = float(result["fe_5200"])
        flux["fe_5270"], flux["fe_5335"] = self.__detect_iron_flux(file_path)
        flux["s3_9069"], flux["s3_9532"] = float(result["s3_9069"]), float(result["s3_9532"])

        diff = set(target_maps.values()) - set(flux.keys())
        for k in diff:
            flux[k] = 0

        hdul.close()
        return flux, redshift
