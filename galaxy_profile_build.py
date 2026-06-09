#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on May  13 18:30:21 2026

@author: sujith-t
"""

import os
import astropy.units as u
import logging
import time

from astroquery.sdss import SDSS
from astropy import coordinates as coords

from dotenv import load_dotenv
from libs.util.db import MySQLUtil
from libs.spectra.sdss_extractor import SpectralProfiler
from monitor import calculate_direct_metallicity_exact

load_dotenv()

logging.basicConfig(level=os.getenv("LOG_LEVEL"), format="%(asctime)s \t %(name)s \t %(levelname)s \t %(message)s")

# initializing the logging in the main script
logger = logging.getLogger(__name__)

if os.getenv("EXEC_ENABLED") != "1":
    logger.warning("Execution halted. To enable set 1 in .env")
    exit()

def populate_galaxy_spectra_flux(start_position=0, no_records=50000):
    start_time = time.time()
    db_util = MySQLUtil(os)
    profiler = SpectralProfiler()

    q = "SELECT obj_id, ra, declination FROM galaxy_catalog WHERE plate_id IS NULL LIMIT %s, %s"
    rows = db_util.fetch_all(q, [start_position, no_records])

    # download
    def __download_spec_file__(id, ra, dec):
        pos = coords.SkyCoord(ra, dec, unit="deg")
        logger.info(f"Processing galaxy {id} at position ({ra}, {dec})")

        # Query spectroscopy (DR19)
        spec = SDSS.query_region(pos, radius=2*u.arcsec, spectro=True, data_release=19)
        if spec is None:
            q = "UPDATE galaxy_catalog SET plate_id = %s WHERE obj_id = %s"
            p = [-1, id]
            db_util.execute(q, p, commit=True)
            logger.debug(f"Couldn't find the matching plate_id, mjd and fiber_id for galaxy {id}")
            return

        file_name = id + ".fits"
        try:
            # Download FITS file
            sp = SDSS.get_spectra(matches=spec)
            sp[0].writeto(file_name, overwrite=True)
        except Exception as e:
            logger.error(f"Failed to download FITS file for galaxy {id}: {e}")
            return

        return spec

    # extract and save flux
    def __extract_save_flux__(id, spec):
        file_name = id + ".fits"
        q = "UPDATE galaxy_catalog SET fiber_id = %s, plate_id = %s, mjd = %s WHERE obj_id = %s"
        p = [int(spec["fiberID"].value[0]), int(spec["plate"].value[0]), int(spec["mjd"].value[0]), id]
        db_util.execute(q, p)

        cef, z = profiler.corrected_emission_flux(file_name)
        q = ("INSERT INTO galaxy_spectra_flux (obj_id, h_alpha, h_alpha_observed, h_beta, h_beta_observed, o3_5007, "
             "o3_4959, o3_4363, o2_3727, n2_6583, n2_6548, s2_6716, s2_6730, "
             "s3_9069, s3_9532, ne_3868, he_4685, fe_5200, fe_5270, fe_5335, o3_4363_err) "
             "VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)")
        p = [id, cef["h_alpha"], cef["h_alpha_observed"], cef["h_beta"], cef["h_beta_observed"], cef["o3_5007"],
             cef["o3_4959"], cef["o3_4363"], cef["o2_3727"], cef["n2_6583"], cef["n2_6548"], cef["s2_6716"], cef["s2_6730"],
             cef["s3_9069"], cef["s3_9532"], cef["ne_3868"], cef["he_4685"], cef["fe_5200"], cef["fe_5270"], cef["fe_5335"], cef["o3_4363_err"]]
        db_util.execute(q, p)

        q = "UPDATE sdss_meta SET redshift = %s WHERE obj_id = %s"
        p = [z, id]
        db_util.execute(q, p, commit=True)

        # Delete FITS file after processing
        if os.path.exists(file_name):
            os.remove(file_name)


    # main execution controlled here
    for obj_id, obj_ra, obj_dec in rows:
        spectra = __download_spec_file__(obj_id, obj_ra, obj_dec)

        if spectra is None:
            continue

        __extract_save_flux__(obj_id, spectra)

    db_util.close()
    end_time = time.time()
    logger.info(f"Processing completed in {end_time - start_time:.2f} seconds for {no_records} records")


def proximate_metallicity_profile(start_position=0, no_records=50000):
    start_time = time.time()
    db_util = MySQLUtil(os)
    profiler = SpectralProfiler()

    q = "SELECT * FROM galaxy_spectra_flux WHERE obj_id NOT IN (SELECT obj_id FROM metallicity_profile) LIMIT %s, %s"
    rows = db_util.fetch_all(q, [start_position, no_records], col_names=True)

    for r in rows:
        q = "SELECT redshift FROM sdss_meta WHERE obj_id = %s"
        z, = db_util.fetch_one(q, [r['obj_id']])
        metallicity = profiler.element_abundance_profile(r, z)
        #metallicity = calculate_direct_metallicity_exact(r['o3_4363'], r['o3_4959'], r['o3_5007'], r['o2_3727'], r['s2_6716'], r['s2_6730'], r['h_beta'])
        print(f"Metallicity for obj_id {r['obj_id']}: {metallicity}")


# now invoke
option = os.getenv("EXEC_LEVEL").lower()
if "f" in option:
    logger.info("Starting to populate spectroscopic data")
    populate_galaxy_spectra_flux(no_records=1)

if "m" in option:
    logger.info("Metallicity calculation is commencing")
    proximate_metallicity_profile(no_records=1)
#id = "587722981741363294"
#file_path = id + ".fits"
#profiler = SpectralProfiler()
#flux, z = profiler.corrected_emission_flux(file_path)

#x = profiler.element_abundance_profile(flux)
#x = profiler.star_formation_rate(flux['h_alpha_observed'], flux['h_beta_observed'], z)
#print(flux[["o2_3727", "ne_3868"]])

# https://www.google.com/search?q=how+do+i+estimate+sulphur+in+a+galaxy+with+SDSS+DR19+dataset&gs_lcrp=EgZjaHJvbWUyBggAEEUYOTIHCAEQIRiPAjIHCAIQIRiPAtIBCTI5NjkzajBqN6gCALACAA&sourceid=chrome&ie=UTF-8&udm=50&fbs=ADc_l-aN0CWEZBOHjofHoaMMDiKpaEWjvZ2Py1XXV8d8KvlI3sbM0Xv-BZKE_VrZb6-djVhOtEikQxVPaqfTMUjZRiWtTDeRmGZV61HbyZRtNy50cjiZ3QFCBrqF00rks-mkLZ5cyenQSeK2Ozzsqht8OSu_KYQHCBA3WlMuTfBn5oe0ikwz2bViWzFTkECg2EsGQaXbXKBDD1XwRgZ7p-ssYEfbcSQUVw&aep=10&ntc=1&mstk=AUtExfCV-aUzvKtsgXFNGl4LslJyWgTFc1FH2Q0skIEP39DgdhQ2LieMUbxLUSWvEWRXfjqDNbRmhcvT5FqsK7XqNQ-SKurq7qv05KgVsn9nLUHOGb-ZIowniexI_JOuufUwsyccse8aBVzxgAQ5wBcjENWvJJxy2BqFidMwPIeY1OX6bTH3mrvFkmSPGYYO_J0zbdH7eoMuQX4OEN7cXTfFSSNXBalgsJP5ZOMUr_7QBkPpefnxolvjjHDyOQ&aioh=3&csuir=1&cs=0&mtid=UGoVaueVE9z0juMPr5-S0AE
