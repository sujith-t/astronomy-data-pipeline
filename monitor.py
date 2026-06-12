import os
import numpy as np

from astroquery.sdss import SDSS
from astropy.io import fits

def file_inspect_errors(mjd, plate_id, fiber_id):
    file_name =  f"{mjd}-{plate_id}-{fiber_id}.fits"

    if not os.path.exists(file_name):
        try:
            # Download FITS file
            sp = SDSS.get_spectra(plate=plate_id, mjd=mjd,
                                  fiberID=fiber_id, data_release=19)
            sp[0].writeto(file_name, overwrite=True)
        except Exception as e:
            print(f"Failed to download FITS file for galaxy : {e}")
            return

    ZWARNING_FLAGS = {
        1 << 0: "SKY_FIBER - Fiber is a sky fiber, not a galaxy",
        1 << 1: "TOO_MANY_OUTLIERS - Too many points rejected in fit",
        1 << 2: "NO_LINES - No emission/absorption lines found",
        1 << 3: "SMALL_DELTA_CHI2 - Chi-squared minimum is too shallow",
        1 << 4: "NEGATIVE_MODEL - Best-fit stellar template went negative",
        1 << 5: "UNPLUGGED - Physical fiber was not plugged in",
        1 << 6: "BAD_TARGET - Target was flagged as bad during selection",
        1 << 7: "NODATA - Missing or corrupted pixel data",
        1 << 8: "MULTIPLE_AMBIGUOUS - Multiple distinct redshift solutions fit well"
    }

    SPPIXMASK_FLAGS = {
        1 << 16: "BADTRACE - Bad profile fit / trace",
        1 << 17: "BADCHIP - Bad pixel on CCD detector",
        1 << 18: "BADANDMASK - Masked out in all individual exposures",
        1 << 19: "BADORMASK - Masked out in at least one individual exposure",
        1 << 20: "SCATTEREDLIGHT - Scattered light significant at this pixel",
        1 << 21: "CROSSTALK - Electronic cross-talk from another fiber",
        1 << 22: "NOSKY - Sky level unknown at this wavelength (ivar=0)",
        1 << 23: "BRIGHTSKY - Sky level is extremely bright relative to object"
    }

    hdul = fits.open(file_name)
    header = hdul[0].header
    coadd_data = hdul[1].data
    flux = coadd_data['flux']
    loglam = coadd_data['loglam']
    ivar = coadd_data['ivar']
    and_mask = coadd_data['and_mask']
    or_mask = coadd_data['or_mask']

    zwarning = header.get("ZWARNING", None)
    if zwarning is None and len(hdul) > 2:
        summary_data = hdul[2].data
    if "ZWARNING" in summary_data.names:
        zwarning = summary_data["ZWARNING"][0]

    if zwarning is None:
            print("Could not find ZWARNING in this FITS file structure.")
    elif zwarning == 0:
        print("⚡ ZWARNING = 0: No problems detected! The pipeline thinks this spectrum is clean.")
    else:
        print(f"⚠️ ZWARNING = {zwarning}: Issues found! Decoding bitmask...\n")
        print("Triggered Flags:")

        # Use bitwise AND (&) to check which flags are flipped to '1'
        for bit_value, description in ZWARNING_FLAGS.items():
            if zwarning & bit_value:
                print(f"  - [Bit Value {bit_value}]: {description}")

    wavelengths = 10 ** loglam

    negative_indices = np.where(flux < 0)[0]

    print(f"Total pixels in spectrum: {len(flux)}")
    print(f"Number of pixels with negative flux: {len(negative_indices)}\n")

    if len(negative_indices) > 0:
        print(f"{'Observed Wave (Å)':<20}{'Flux Value':<15}{'Noise (1/sqrt(ivar))':<22}{'Triggered Bitmask Flags'}")
    print("-" * 90)

    # Check the first 10 negative occurrences as a sample
    errors = []
    for idx in negative_indices[:10]:
        wave = wavelengths[idx]
        flx = flux[idx]

        # Calculate error bar. If ivar is 0, noise is infinite (completely untrustworthy)
        err = 1.0 / np.sqrt(ivar[idx]) if ivar[idx] > 0 else float('inf')

        # We test 'and_mask' first as it is the most robust pixel mask
        pixel_mask_value = and_mask[idx]

        # Decode what flags are raised for this specific negative pixel
        triggered_flags = []
        for bit_value, description in SPPIXMASK_FLAGS.items():
            if pixel_mask_value & bit_value:
                # Extract just the short name from our description dictionary
                triggered_flags.append(description.split(" - ")[0])

        # Format flags string
        flags_str = ", ".join(triggered_flags) if triggered_flags else "None (Routine statistical noise)"

        print(f"{wave:<20.2f}{flx:<15.4f}{err:<22.4f}{flags_str}")
        issue = {"wave": wave, "flux": flx, "error": err, "flags": flags_str}
        errors.append(issue)

    if len(negative_indices) > 10:
        print(f"\n... and {len(negative_indices) - 10} more negative pixels.")

    return errors, file_name


def detect_instrument_errors(target: list[dict], file_name:str):
    wl = [float(r["wave"]) for r in target]

    hdul = fits.open(file_name)
    #header = hdul[0].header
    coadd_data = hdul[1].data
    flux = coadd_data['flux']
    loglam = coadd_data['loglam']
    ivar = coadd_data['ivar']
    and_mask = coadd_data['and_mask']
    #or_mask = coadd_data['or_mask']

    wavelengths = 10 ** loglam

    for tw in wl:
        # Find the pixel closest to your target wavelength
        closest_idx = np.argmin(np.abs(wavelengths - tw))

        found_wave = wavelengths[closest_idx]
        flx = flux[closest_idx]
        pixel_ivar = ivar[closest_idx]
        mask_val = and_mask[closest_idx]

        # Test for SDSS's specific blue-end calibration flags
        triggered_problems = []

        # Check if the pipeline already flagged it as a bad CCD chip zone
        if mask_val & (1 << 17):
            triggered_problems.append("BADCHIP (Hardware defect)")
        # Check if it was masked out in individual exposures
        if mask_val & (1 << 18):
            triggered_problems.append("BADANDMASK (Unreliable data)")
        # Check if the pipeline thinks sky subtraction broke
        if mask_val & (1 << 22):
            triggered_problems.append("NOSKY (Sky calibration failed)")

        # If no flags are tripped but ivar is tiny, it's a S/N limitation
        if not triggered_problems and pixel_ivar < 0.01:
            triggered_problems.append("Pure Statistical Noise (Low S/N at Blue Edge)")

        problems_str = ", ".join(triggered_problems) if triggered_problems else "Clean Pixel (Pure Random Variance)"

        print(f"{tw:<12} {found_wave:<20.2f} {flx:<10.3f} {pixel_ivar:<10.3f} {problems_str}")


# test file and values
errors, file_name = file_inspect_errors(52000, 288, 137)
detect_instrument_errors(errors, file_name)
