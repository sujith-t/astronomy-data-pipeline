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

load_dotenv()

logging.basicConfig(level=os.getenv("LOG_LEVEL"), format="%(asctime)s \t %(name)s \t %(levelname)s \t %(message)s")

# initializing the logging in the main script
logger = logging.getLogger(__name__)

def populate_galaxy_spectra_flux(start_position=0, no_records=50000):
    start_time = time.time()
    db_util = MySQLUtil(os)
    profiler = SpectralProfiler()

    q = "SELECT obj_id, ra, declination FROM galaxy_catalog WHERE plate_id IS NULL LIMIT %s, %s"
    rows = db_util.fetch_all(q, [start_position, no_records])

    for id, ra, dec in rows:
        pos = coords.SkyCoord(ra, dec, unit="deg")
        logger.info(f"Processing galaxy {id} at position ({ra}, {dec})")

        # Query spectroscopy (DR19)
        spec = SDSS.query_region(pos, radius=2*u.arcsec, spectro=True, data_release=19)
        if spec is None:
            q = "UPDATE galaxy_catalog SET plate_id = %s WHERE obj_id = %s"
            p = [-1, id]
            db_util.execute(q, p, commit=True)
            logger.debug(f"Couldn't find the matching plate_id, mjd and fiber_id for galaxy {id}")
            continue

        file_name = id + ".fits"
        try:
            # Download FITS file
            sp = SDSS.get_spectra(matches=spec)
            sp[0].writeto(file_name, overwrite=True)
        except Exception as e:
            logger.error(f"Failed to download FITS file for galaxy {id}: {e}")
            continue

        q = "UPDATE galaxy_catalog SET fiber_id = %s, plate_id = %s, mjd = %s WHERE obj_id = %s"
        p = [int(spec["fiberID"].value[0]), int(spec["plate"].value[0]), int(spec["mjd"].value[0]), id]
        db_util.execute(q, p)

        cef, z = profiler.corrected_emission_flux(file_name)
        q = ("INSERT INTO galaxy_spectra_flux (obj_id, h_alpha, h_alpha_observed, h_beta, h_beta_observed, o3_5007, "
             "o3_4959, o2_3727, n2_6583, n2_6548, s2_6716, s2_6730, "
             "s3_9069, s3_9532, ne_3868, he_4685, fe_5200, fe_5270, fe_5335) "
             "VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)")
        p = [id, cef["h_alpha"], cef["h_alpha_observed"], cef["h_beta"], cef["h_beta_observed"], cef["o3_5007"],
             cef["o3_4959"], cef["o2_3727"], cef["n2_6583"], cef["n2_6548"], cef["s2_6716"], cef["s2_6730"],
             cef["s3_9069"], cef["s3_9532"], cef["ne_3868"], cef["he_4685"], cef["fe_5200"], cef["fe_5270"], cef["fe_5335"]]
        db_util.execute(q, p)

        q = "UPDATE sdss_meta SET redshift = %s WHERE obj_id = %s"
        p = [z, id]
        db_util.execute(q, p, commit=True)

        # Delete FITS file after processing
        if os.path.exists(file_name):
            os.remove(file_name)

    db_util.close()
    end_time = time.time()
    logger.info(f"Processing completed in {end_time - start_time:.2f} seconds for {no_records} records")

# now invoke
populate_galaxy_spectra_flux(no_records=50)
