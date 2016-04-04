from astropy.io import fits
from ossos import wcs
from ossos import storage
from pyraf import iraf
import os

def make_iraf_psf(fname):

    hdulist = fits.open(fname)

    # Pixel value amplitude of brightest source, rough guess of brightness of asteroid for passing to Trippy
    amplitude = hdulist[1].data.max() - hdulist[1].data.mean()

    with fits.open(fname) as stp:   # Read in image infor placed in header previously
        header = stp[0].header
        objnum = header['OBJNUM']
        ra = header['RA']
        dec = header['DEC']
        ra_rate = header['RARATE']
        dec_rate = header['DECRATE']

        for ext in range(len(stp)): # Look in each header section for these values
            try:
                header = stp[ext].header
                expnum = header['EXPNUM']
                exptime = header['EXPTIME']
                zpt = header['PHOTZP']
                pixscale = header['CD1_1']*3600
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

        with fits.open(psf_file) as p:
            for ext in range(len(p)): # Look in each header section for radius of PSF
                try:
                    header = p[ext].header
                    psfrad = header['PSFRAD']
                except:
                    continue

        iraf.noao()
        iraf.digiphot()
        iraf.daophot(_doprint=0)
        psf_image = psf_file.replace('psf','seepsf')

        # Tell seepsf to create a psf image of size 2*psfrad, +5 pixels for padding since Trippy assumes larger images
        # +5 also ensures an image with dimensions of an odd number, which is what Trippy wants
        iraf.seepsf(psf_file, psf_image, xpsf=image_coords[0], ypsf=image_coords[1], dimension = 2*psfrad+5)
        os.unlink(psf_file)

        return objnum, expnum, exptime, zpt, pixscale, ra_rate, dec_rate, amplitude, psfrad, psf_image

    else:
        return objnum, expnum, exptime, zpt, pixscale, ra_rate, dec_rate, amplitude, psfrad, None
