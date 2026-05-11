/* alter table galaxy_catalog */
ALTER TABLE galaxy_catalog
ADD plate_id VARCHAR(10),
ADD mjd VARCHAR(10),
ADD fiber_id VARCHAR(10);

/* creating table to save dust bias corrected spectra flux values */
CREATE TABLE galaxy_spectra_flux (
    obj_id VARCHAR(20) NOT NULL,
    h_alpha DOUBLE DEFAULT NULL,
    h_beta DOUBLE DEFAULT NULL,
    o3_5007 DOUBLE DEFAULT NULL,
    o3_4959 DOUBLE DEFAULT NULL,
    o2_3727 DOUBLE DEFAULT NULL,
    n2_6584 DOUBLE DEFAULT NULL,
    n2_6548 DOUBLE DEFAULT NULL,
    s2_6716 DOUBLE DEFAULT NULL,
    s2_6731 DOUBLE DEFAULT NULL,
    he_3869 DOUBLE DEFAULT NULL,
    he_4686 DOUBLE DEFAULT NULL,
    fe_5200 DOUBLE DEFAULT NULL,
    PRIMARY KEY (`obj_id`)
);
