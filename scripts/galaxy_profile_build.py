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

def build_galaxy_profile(start_position=0, no_records=50000):
    db_util = MySQLUtil(os)
    profiler = SpectralProfiler()

    q = "SELECT obj_id, ra, declination FROM galaxy_catalog WHERE plate_id IS NULL LIMIT %s, %s"
    rows = db_util.fetch_all(q, [start_position, no_records])

    for id, ra, dec in rows:
        pos = coords.SkyCoord(ra, dec, unit="deg")

        # Query spectroscopy (DR17)
        spec = SDSS.query_region(pos, radius=2*u.arcsec, spectro=True)

        file_name = id + ".fits"
        try:
            # Download FITS file
            sp = SDSS.get_spectra(matches=spec)
            sp[0].writeto(file_name, overwrite=True)
        except Exception as e:
            print(e)
            continue

        q = "UPDATE galaxy_catalog SET fiber_id = %s, plate_id = %s, mjd = %s WHERE obj_id = %s"
        p = [int(spec["fiberID"].value[0]), int(spec["plate"].value[0]), int(spec["mjd"].value[0]), id]
        db_util.execute(q, p)



        corrected = profiler.dust_bias_correction(profiler.detect_emission_flux(file_name))
        q = "INSERT INTO galaxy_spectra_flux VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)"
        p = [id, corrected["H_Alpha"], corrected["H_Beta"], corrected["O3_5007"], corrected["O3_4959"], corrected["O2_3727"],
             corrected["N2_6584"], corrected["N2_6548"], corrected["S2_6716"], corrected["S2_6731"], corrected["Ne_3869"], corrected["He_4686"], corrected["Fe_5200"]]
        db_util.execute(q, p)

        result = profiler.element_abundance_profile(file_name)
        # the factor would mean the 10^6 atoms of Hydrogen ex. for 309.43 oxygen atoms per 10^6 hydrogen atoms
        factor = 10 ** 6
        oxygen = float(factor * result["oxygen"])
        nitrogen = float(factor * result["nitrogen"])
        carbon = float(factor * result["carbon"])
        sulphur = float(factor * result["sulphur"])
        neon = float(factor * result["neon"])
        iron_strength = float(result["iron_strength"])
        star_formation_rate = profiler.star_formation_rate(file_name) #77118407393508496

        q = "INSERT INTO metallicity_profile VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)"
        p = [id, result["metallicity_O3N2"], result["metallicity_R23"], result["final_metallicity"], oxygen, nitrogen,
             carbon, sulphur, neon, iron_strength, star_formation_rate]
        db_util.execute(q, p, commit=True)

        # Delete FITS file after processing
        if os.path.exists(file_name):
            os.remove(file_name)

        print(id, " done")

    db_util.close()


# now invoke
build_galaxy_profile(no_records=10000)