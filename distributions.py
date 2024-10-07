from astropy.io import fits
import astropy.coordinates as coord
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

coord.galactocentric_frame_defaults.set("v4.0")

d = fits.open('mwsall-pix-fuji.fits')
#d = fits.open('H06_030deg_mock.fits')
#d.info()

s = d['SPTAB'].data
m = d['FIBERMAP'].data
g = d['GAIA'].data

w = (g['PARALLAX']/g['PARALLAX_ERROR'] > 10.) & (g['PARALLAX'] > 0.)
lenw = len(np.where(w)[0])

icrs = coord.SkyCoord(
    ra = m['TARGET_RA'][w] * u.deg, 
    dec = m['TARGET_DEC'][w] * u.deg, 
    distance = 1./g['PARALLAX'][w] * u.kpc,
    pm_ra_cosdec = g['PMRA'][w] * u.mas / u.yr,
    pm_dec = g['PMDEC'][w] * u.mas / u.yr,
    radial_velocity = s['RV_ADOP'][w] * u.km / u.s 
    )

galcen_frame = coord.Galactocentric()
galcen = icrs.transform_to(galcen_frame)

plt.plot(galcen.v_x,galcen.v_y,',')
plt.xlabel('Vx (km/s)')
plt.ylabel('Vy (km/s)')
plt.xlim((-500,500))
plt.ylim((-500,500))
plt.title('Gaia sources with positive parallax good to 10%')
plt.show()

w = (s['TEFF'] > 4000.)  & (s['FEH'] > -4.)
plt.clf()
plt.hist(s['feh'][w],bins=100)
plt.xlim((-4.9,0.5))
plt.xlabel('[Fe/H]')
plt.title('Sources with Teff>4000 K and [Fe/H]> -4.')
plt.show()

