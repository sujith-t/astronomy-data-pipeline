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
    db_util.execute(q, p)

    corrected = profiler.dust_bias_correction(profiler.detect_emission_flux(file_name))
    q = "INSERT INTO galaxy_spectra_flux VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)"
    p = [id, corrected["H_Alpha"], corrected["H_Beta"], corrected["O3_5007"], corrected["O3_4959"], corrected["O2_3727"],
         corrected["N2_6584"], corrected["N2_6548"], corrected["S2_6716"], corrected["S2_6731"], corrected["Ne_3869"], corrected["He_4686"], corrected["Fe_5200"]]
    db_util.execute(q, p, commit=True)

    result = profiler.element_abundance_profile(file_name)
    # the factor would mean the 10^6 atoms of Hydrogen ex. for 309.43 oxygen atoms per 10^6 hydrogen atoms
    factor = 10 ** 6
    oxygen = factor * result["oxygen"],
    nitrogen = factor * result["nitrogen"]
    carbon = factor * result["carbon"]
    sulphur = factor * result["sulphur"]
    neon = factor * result["neon"]
    iron_strength = result["iron_strength"]

    q = "INSERT INTO metallicity_profile VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)"
    p = [id, result["metallicity_O3N2"], result["metallicity_R23"], result["final_metallicity"], oxygen, nitrogen,
         carbon, sulphur, neon, iron_strength]
    #db_util.execute(q, p, commit=True)



    #result = profiler.dust_bias_correction(result)
    print(id, " done")
    # Delete FITS file after processing
    #if os.path.exists(file_name):
        #os.remove(file_name)

db_util.__close__()
