#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 15:34:26 2026

@author: sujith-t
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.polynomial import Polynomial
from scipy.optimize import curve_fit
from astropy.cosmology import Planck18 as cosmo


class SpectralProfiler:

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


    def __fit_line(self, w, f, center, window=10):

        def __gaussian__(x, amp, mu, sigma):
            return amp * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))

        mask = (w > center - window) & (w < center + window)
        x = w[mask]
        y = f[mask]

        # local continuum = np.median(y)
        y = y - np.median(y)

        if len(x) < 5:
            return 0

        # Initial guesses
        amp_guess = max(np.max(y), 1e-3)
        mu_guess = center
        sigma_guess = 2

        try:
            popt, _ = curve_fit(__gaussian__, x, y, p0=[amp_guess, mu_guess, sigma_guess], bounds=([0, center - 3, 0.5], [np.inf, center + 3, 10]))
            amp, mu, sigma = popt

            # Flux = area under Gaussian
            flux_line = amp * sigma * np.sqrt(2 * np.pi)

            return flux_line
        except:
            return 0


    def __find_rest_wavelength(self, file_path: str):
        hdul = fits.open(file_path)

        # data retreval from the file
        data = hdul[1].data
        flux = data['flux']
        loglam = data['loglam']
        ivar = data['ivar']  # inverse variance

        wavelength = 10 ** loglam

        # Convert ivar to error
        error = np.zeros_like(ivar)
        mask = ivar > 0
        error[mask] = 1 / np.sqrt(ivar[mask])

        # redshift correction
        z = hdul[2].data['Z'][0]
        rest_wavelength = wavelength / (1 + z)

        return rest_wavelength, flux

    def __detect_emission_flux(self, file_path: str, detection_lines: dict) -> dict:
        rest_wavelength, flux = self.__find_rest_wavelength(file_path)

        # Remove stellar background -> in a real spectra steller background is also embedded
        # Fit continuum (exclude emission regions roughly)
        p = Polynomial.fit(rest_wavelength, flux, deg=3)
        continuum = p(rest_wavelength)
        flux_cont_sub = flux - continuum

        if len(detection_lines) == 0:
            detection_lines = {
                "h_alpha": 6563,  # Hα star formation
                "h_beta": 4861,
                "o3_5007": 5007,  # Metallicity
                "o3_4959": 4959,
                "o2_3727": 3727,
                "n2_6584": 6584,  # Metallicity
                "n2_6548": 6548,
                "s2_6716": 6716,  # Density
                "s2_6731": 6731,  # Density
                "ne_3869": 3869,  # Ionization
                "he_4686": 4686,  # Ho star/AGN
                "fe_5200": 5200 # Stellar population
            }

        results = {}
        for name, center in detection_lines.items():
            flux_val = self.__fit_line(rest_wavelength, flux_cont_sub, center)
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
            ebv = 1.97 * np.log10(observed_ratio / 2.86)
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

        return fe_5270, fe_5335


    # the ratio is in reltion to Hydrogen as the baseline
    def element_abundance_profile(self, corrected:dict):
        ratios = {}

        if "h_alpha" not in corrected or "h_beta" not in corrected:
            return

        ##### 1. OXYGEN and Matalicity #####
        # compute R23
        total_oxygen = corrected["o3_5007"] + corrected["o3_4959"] + corrected["o2_3727"]
        oxygen_r23 = 0
        if total_oxygen != 0 and  corrected["h_beta"] != 0:
            oxygen_r23 = total_oxygen / corrected["h_beta"]

        # compute O3N2 | metallicity_O3N2 > 8.4 = high otherwise low
        o3n2 = 0
        if corrected["h_beta"] > 0 and corrected["h_alpha"] > 0:
            o3n2 = np.log10((corrected["o3_5007"] / corrected["h_beta"]) / (corrected["n2_6584"] / corrected["h_alpha"]))
        ratios["metallicity_o3n2"] = 8.73 - 0.32 * o3n2

        ratios["metallicity_r23"] = 7.5 + 0.8 * np.log10(oxygen_r23)
        if ratios["metallicity_o3n2"] > 8.4:
            ratios["metallicity_r23"] = 9.2 - 0.3 * np.log10(oxygen_r23)

        # weighted combined metallicity - (research + production grade)
        # interpretation -> 12 + log(O/H) = final_metallicity
        ratios["final_metallicity"] = 0.6 * ratios["metallicity_o3n2"] + 0.4 * ratios["metallicity_r23"]

        # log(O/H) = final_metallicity - 12
        log_oh = ratios["final_metallicity"] - 12
        ratios["oxygen"] = 10 ** log_oh

        ##### 2. NITROGEN #####
        # calculate log(N/O) ratio and add a calibration offset of 0.05
        log_no = 0.05
        if corrected["o2_3727"] > 0:
            log_no = np.log10(corrected["n2_6584"] / corrected["o2_3727"]) + log_no

        # log(N/H) = log(N/O)+log(O/H)
        log_nh = log_no + log_oh
        ratios["nitrogen"] = 10 ** log_nh

        ##### 3. CARBON ##### this is an indirect estimate, flux lines are weak and complicated. only Ultra-Violet files can catch them

        # plan log(C/O) = a + b × (12+log(O/H)) where a, b are calibration values | log(C/H) = log(C/O) + log(O/H)
        log_co = -0.8 + 0.14 * (ratios["final_metallicity"] - 8.0)
        log_ch = log_co + log_oh
        ratios["carbon"] = 10 ** log_ch

        ##### 4. SULPHUR #####
        sulphur_oxygen_ratio = 0.025
        ratios["sulphur"] = ratios["oxygen"] * sulphur_oxygen_ratio

        ##### 5. estimate NEON #####
        log_neo =  0.7
        if corrected["o3_5007"] != 0 or corrected["o3_4959"]:
            log_neo = np.log10(corrected["ne_3869"] / (corrected["o3_5007"] + corrected["o3_4959"])) + log_neo

        log_neh = log_neo + log_oh
        ratios["neon"] = 10 ** log_neh

        ##### 6. Iron (Fe) ##### this can't be interpreted as ratio Fe/H
        fe_index = (corrected['fe_5270'] + corrected['fe_5335']) / 2
        ratios["iron_strength"] = -2.0 + 0.4 * fe_index

        return ratios

    def corrected_emission_flux(self, file_path: str):
        hdul = fits.open(file_path)
        target_maps = {'H_beta': 'h_beta', 'H_alpha': 'h_alpha', '[O_III] 5007': 'o3_5007', '[O_III] 4959': 'o3_4959', '[O_II] 3727': 'o2_3727', '[N_II] 6583': 'n2_6584', '[N_II] 6548': 'n2_6548', '[S_II] 6716': 's2_6716', '[S_II] 6730': 's2_6731', '[Ne_III] 3868': 'ne_3869', 'He_II 4685': 'he_4686'}
        line_data = hdul["SPZLINE"].data
        redshift = hdul[2].data['Z'][0]

        flux = {}
        for r in line_data:
            if r['LINENAME'].strip() in target_maps:
                flux[target_maps[r['LINENAME'].strip()]] = r['LINEAREA']

        # this step is done due to unavailability of iron in SPZLINE (not already calculated)
        # rest of them already calculated, no need to do again
        wavelengths = {"h_alpha": 6563, "h_beta": 4861, "fe_5200": 5200}
        result = self.__detect_emission_flux(file_path, wavelengths)
        flux["h_alpha_obs"] = result["h_alpha"]
        flux["h_beta_obs"] = result["h_beta"]
        result = self.__dust_bias_correction(result, h_alpha=result["h_alpha"], h_beta=result["h_beta"])
        flux["fe_5200"] = result["fe_5200"]
        flux["fe_5270"], flux["fe_5335"] = self.__detect_iron_flux(file_path)

        diff = set(target_maps.values()) - set(flux.keys())
        for k in diff:
            flux[k] = 0

        hdul.close()
        return flux, redshift
