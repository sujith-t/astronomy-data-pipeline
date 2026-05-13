#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on May  13 18:30:21 2026

@author: sujith-t
"""

import os
import astropy.units as u

from astroquery.sdss import SDSS
from astropy import coordinates as coords

from dotenv import load_dotenv
from libs.util.db import MySQLUtil
from libs.spectra.extractor import SpectralProfiler

load_dotenv()

start_position = 0
no_records = 1


# Save first spectrum


db_util = MySQLUtil(os)
profiler = SpectralProfiler()

q = "SELECT obj_id, ra, declination FROM galaxy_catalog LIMIT %s, %s"
rows = db_util.fetch_all(q, [start_position, no_records])

for id, ra, dec in rows:
    pos = coords.SkyCoord(ra, dec, unit="deg")

    # Query spectroscopy (DR17)
    spec = SDSS.query_region(pos, radius=2*u.arcsec, spectro=True)

    # Download FITS file
    file_name = id + ".fits"
    #sp = SDSS.get_spectra(matches=spec)
    #sp[0].writeto(file_name, overwrite=True)

    q = "UPDATE galaxy_catalog SET fiber_id = %s, plate_id = %s, mjd = %s WHERE obj_id = %s"
    p = [int(spec["fiberID"].value[0]), int(spec["plate"].value[0]), int(spec["mjd"].value[0]), id]
    db_util.execute(q, p, commit=True)

    result = profiler.element_abundance_profile(file_name)
    #result = profiler.dust_bias_correction(result)
    print(result)
    # Delete FITS file after processing
    #if os.path.exists(file_name):
        #os.remove(file_name)

db_util.__close__()
