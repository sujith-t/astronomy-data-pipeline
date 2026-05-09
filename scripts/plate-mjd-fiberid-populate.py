from astroquery.sdss import SDSS
from astropy import coordinates as coords
import astropy.units as u

# Your coordinates
ra = 182.92526245117188
dec = -1.092357039

pos = coords.SkyCoord(ra, dec, unit="deg")

# Query spectroscopy (DR17)
spec = SDSS.query_region(
    pos,
    radius=2*u.arcsec,
    spectro=True
)
#plate, mjd, fiberID
# Download FITS file
sp = SDSS.get_spectra(matches=spec)

# Save first spectrum
sp[0].writeto("galaxy_spectrum.fits", overwrite=True)