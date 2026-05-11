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

file_location = "/var/project/astronomy-data-pipeline/scripts/galaxy_spectrum.fits"


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


    def fit_line(self, w, f, center, window=10):

        def __gaussian__(x, amp, mu, sigma):
            return amp * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))

        mask = (w > center - window) & (w < center + window)
        x = w[mask]
        y = f[mask]

        if len(x) < 5:
            return None

        # Initial guesses
        amp_guess = np.max(y)
        mu_guess = center
        sigma_guess = 2

        try:
            popt, _ = curve_fit(__gaussian__, x, y, p0=[amp_guess, mu_guess, sigma_guess])
            amp, mu, sigma = popt

            # Flux = area under Gaussian
            flux_line = amp * sigma * np.sqrt(2 * np.pi)

            return flux_line
        except:
            return None


    def find_rest_wavelength(self, file_path: str):
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


    def detect_emission_flux(self, file_path: str) -> dict:
        rest_wavelength, flux = self.find_rest_wavelength(file_path)

        # Remove stellar background -> in a real spectra steller background is also embedded
        # Fit continuum (exclude emission regions roughly)
        p = Polynomial.fit(rest_wavelength, flux, deg=3)
        continuum = p(rest_wavelength)
        flux_cont_sub = flux - continuum

        detection_lines = {
            "H_Alpha": 6563,  # Hα str formation
            "H_Beta": 4861,
            "O3_5007": 5007,  # Metallicity
            "O3_4959": 4959,
            "O2_3727": 3727,
            "N2_6584": 6584,  # Metallicity
            "N2_6548": 6548,
            "S2_6716": 6716,  # Density
            "S2_6731": 6731,  # Density
            "Ne_3869": 3869,  # Ionization
            "He_4686": 4686,  # Ho star/AGN
            "Fe_5200": 5200 # Stellar population
        }

        results = {}
        for name, center in detection_lines.items():
            flux_val = self.fit_line(rest_wavelength, flux_cont_sub, center)
            results[name] = flux_val

        return results


    # using Kennicutt relation to find star formation rate FR(M⊙​/yr)=7.9×10−42×LHα​(erg/s)
    def star_formation_rate(self, file_path: str):
        # read the file
        hdul = fits.open(file_path)

        result = self.detect_emission_flux(file_path)

        observed_ratio = result["H_Alpha"] / result["H_Beta"]
        intrinsic_ratio = 2.86

        e_bv = 2.5 * np.log10(observed_ratio / intrinsic_ratio)

        # simple correction factor
        correction = 10 ** (0.4 * e_bv)

        h_alpha_flux = result["H_Alpha"] * correction

        # redshift
        z = hdul[2].data['Z'][0]

        # Convert redshift → distance Luminosity distance in cm
        luminocity_distance = cosmo.luminosity_distance(z).to('cm').value

        # flu -> luminocity L=4πdL2​×F
        luminocity_h_alpha = 4 * np.pi * luminocity_distance ** 2 * h_alpha_flux

        return 7.9e-42 * luminocity_h_alpha


    # remove dust effects
    def dust_bias_correction(self, emission_flux: dict) -> dict:
        k = 2.5
        intrinsic_ratio = 2.86

        if "H_Alpha" not in emission_flux or "H_Beta" not in emission_flux:
            print("H_Alpha or H_Beta not found in the emission_flux.")
            return

        observed_ratio = emission_flux["H_Alpha"] / emission_flux["H_Beta"]
        e_bv = k * np.log10(observed_ratio / intrinsic_ratio)

        corrected = {}

        # make correction for all lines
        for line, flux in emission_flux.items():
            if flux:
                corrected[line] = flux * 10 ** (0.4 * k * e_bv)
            else:
                corrected[line] = 0

        return corrected


    def detect_approximate_iron(self, file_path: str):
        Fe5270_range = (5245, 5285)
        Fe5335_range = (5315, 5355)

        rest_wavelength, flux = self.find_rest_wavelength(file_path)

        def measure_index(wave, flux, wmin, wmax):
            mask = (wave > wmin) & (wave < wmax)
            return np.mean(flux[mask])

        Fe5270 = measure_index(rest_wavelength, flux, *Fe5270_range)
        Fe5335 = measure_index(rest_wavelength, flux, *Fe5335_range)

        Fe_index = (Fe5270 + Fe5335) / 2
        log_FeH = -2.0 + 0.4 * Fe_index

        return log_FeH


    # the ratio is in reltion to Hydroden as the baseline
    def element_abundance_profile(self, file_path: str):
        ratios = {}

        result = self.detect_emission_flux(file_path)
        corrected = self.dust_bias_correction(result)

        if "H_Alpha" not in corrected or "H_Beta" not in corrected:
            return

        ##### 1. OXYGEN and Matalicity #####
        # compute R23
        total_oxygen = corrected["O3_5007"] + corrected["O3_4959"] + corrected["O2_3727"]
        oxygen_R23 = total_oxygen / corrected["H_Beta"]

        # compute O3N2 | metallicity_O3N2 > 8.4 = high otherwise low
        o3n2 = np.log10((corrected["O3_5007"] / corrected["H_Beta"]) / (corrected["N2_6584"] / corrected["H_Alpha"]))
        ratios["metallicity_O3N2"] = 8.73 - 0.32 * o3n2

        ratios["metallicity_R23"] = 7.5 + 0.8 * np.log10(oxygen_R23)
        if ratios["metallicity_O3N2"] > 8.4:
            ratios["metallicity_R23"] = 9.2 - 0.3 * np.log10(oxygen_R23)

        # weighted combined metalicity - (research + production grade)
        # intepretation -> 12 + log(O/H) = final_metallicity
        ratios["final_metallicity"] = 0.6 * ratios["metallicity_O3N2"] + 0.4 * ratios["metallicity_R23"]

        # log(O/H) = final_metallicity - 12
        log_OH = ratios["final_metallicity"] - 12
        ratios["oxygen"] = 10 ** log_OH

        ##### 2. NITROGEN #####
        # calculate log(N/O) ratio and add a calibration offset of 0.05
        log_NO = np.log10(corrected["N2_6584"] / corrected["O2_3727"]) + 0.05

        # log(N/H) = log(N/O)+log(O/H)
        log_NH = log_NO + log_OH
        ratios["nitrogen"] = 10 ** log_NH

        ##### 3. CARBON ##### this is an indirect estimate, flux lines are weak and complicated. only Ultra-Violet files can catch them

        # plan log(C/O) = a + b × (12+log(O/H)) where a, b are calibration values | log(C/H) = log(C/O) + log(O/H)
        log_CO = -0.8 + 0.14 * (ratios["final_metallicity"] - 8.0)
        log_CH = log_CO + log_OH
        ratios["carbon"] = 10 ** log_CH

        ##### 4. SULFER #####
        log_sulpher_hbeta = np.log10((corrected["S2_6716"] + corrected["S2_6731"]) / corrected["H_Beta"])

        # oxygen ionization rato estimating: ICF(S) ≈ O/O+
        oxygen_ratio = total_oxygen / corrected["O2_3727"]
        log_SH = log_sulpher_hbeta + np.log10(oxygen_ratio)
        ratios["sulpher"] = 10 ** log_SH

        ##### 5. stimate NEON #####
        log_NeO = np.log10(corrected["Ne_3869"] / (corrected["O3_5007"] + corrected["O3_4959"])) + 0.7
        log_NeH = log_NeO + log_OH
        ratios["neon"] = 10 ** log_NeH

        ##### 6. Iron (Fe) ##### this can't be intepretted as ratio Fe/H
        ratios["iron_strength"] = self.detect_approximate_iron(file_path)

        return ratios
