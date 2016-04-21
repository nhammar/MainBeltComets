from astropy.io import fits
from ossos import wcs
from ossos import storage
from pyraf import iraf
import os

def make_iraf_psf(fname):

    with fits.open(fname) as stp:   # Read in image info placed in header previously
        header = stp[0].header
        ra = header['RA']
        dec = header['DEC']
        for ext in range(len(stp)): # Look in each header section
            try:
                header = stp[ext].header
                expnum = header['EXPNUM']
            except:
                continue

    stamp = fits.open(fname)
    for ext in range(len(stamp)):   # Establish size of image and location of source in image (source should be centered in cutout)
        try:
            world_coords = wcs.WCS(stamp[ext].header)
            image_coords = world_coords.sky2xy(float(ra),float(dec))
            image_size = stamp[ext].data.shape
        except:
            continue

    # Some cutouts are split over multiple CCDs, meaning the source may not actually be in some of these small image 'pieces'
    if 0<image_coords[0]<image_size[0] and 0<image_coords[1]<image_size[1]:
       
        ccdnum = stamp[ext].header.get('EXTNO', None) 

        psf_file = storage.get_file(expnum, ccd=ccdnum, ext='psf.fits') # Get psf file for image from storage

        iraf.noao()
        iraf.digiphot()
        iraf.daophot(_doprint=0)
        psf_image = psf_file.replace('psf','seepsf')

        # Seepsf creates psf_image
        iraf.seepsf(psf_file, psf_image, xpsf=image_coords[0], ypsf=image_coords[1])
        os.unlink(psf_file)

        return psf_image

    else:
        return None
