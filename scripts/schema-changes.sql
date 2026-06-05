/* alter table galaxy_catalog */
ALTER TABLE galaxy_catalog
ADD plate_id VARCHAR(10),
ADD mjd VARCHAR(10),
ADD fiber_id VARCHAR(10);

/* creating table to save dust bias corrected spectra flux values */
CREATE TABLE galaxy_spectra_flux (
    obj_id VARCHAR(20) NOT NULL,
    h_alpha DOUBLE DEFAULT NULL,
    h_alpha_observed DOUBLE DEFAULT NULL,
    h_beta DOUBLE DEFAULT NULL,
    h_beta_observed DOUBLE DEFAULT NULL,
    o3_5007 DOUBLE DEFAULT NULL,
    o3_4959 DOUBLE DEFAULT NULL,
    o3_4363 DOUBLE DEFAULT NULL,
    o2_3727 DOUBLE DEFAULT NULL,
    n2_6583 DOUBLE DEFAULT NULL,
    n2_6548 DOUBLE DEFAULT NULL,
    s2_6716 DOUBLE DEFAULT NULL,
    s2_6730 DOUBLE DEFAULT NULL,
    s3_9069 DOUBLE DEFAULT NULL,
    s3_9532 DOUBLE DEFAULT NULL,
    ne_3868 DOUBLE DEFAULT NULL,
    he_4685 DOUBLE DEFAULT NULL,
    fe_5200 DOUBLE DEFAULT NULL,
    fe_5270 DOUBLE DEFAULT NULL,
    fe_5335 DOUBLE DEFAULT NULL,
    PRIMARY KEY (`obj_id`)
);

CREATE TABLE metallicity_profile (
    obj_id VARCHAR(20) NOT NULL,
    o3n2_metallicity DOUBLE DEFAULT NULL,
    r23_metallicity DOUBLE DEFAULT NULL,
    final_metallicity DOUBLE DEFAULT NULL,
    oxygen_hydrogen_ratio DOUBLE DEFAULT NULL,
    nitrogen_hydrogen_ratio DOUBLE DEFAULT NULL,
    carbon_hydrogen_ratio DOUBLE DEFAULT NULL,
    sulphur_hydrogen_ratio DOUBLE DEFAULT NULL,
    neon_hydrogen_ratio DOUBLE DEFAULT NULL,
    iron_strength_index DOUBLE DEFAULT NULL,
    star_form_rate DOUBLE DEFAULT NULL,
    o3_temp_exact DOUBLE DEFAULT NULL,
    PRIMARY KEY (`obj_id`)
);
