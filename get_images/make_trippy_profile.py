import os
import math
from astropy.io import fits
import numpy as np
from trippy import psf, pill, MCMCfit, bgFinder
from sep_phot import sep_phot

import pylab as pyl
from stsci import numdisplay
from astropy.visualization import interval

def make_trippy_profile(psf_image, stamp):

    with fits.open(psf_image) as seepsf:        # Get pixel values and Gaussian fit value par1 from seepsf image
        for ext in range(len(seepsf)):
            try:
                header = seepsf[ext].header
                psfdata = seepsf[ext].data      # PSF image data
                par1 = header['PAR1']
            except:
                continue

    with fits.open(stamp) as stp:       # Postage stamp image
        header = stp[0].header
        ra_rate = header['RARATE']      # Values placed in header ext 0 previously
        dec_rate = header['DECRATE']
        for ext in range(len(stp)):     # Look in each header section for these values
            try:
                header = stp[ext].header
                stampdata = stp[ext].data       # Stamp image data
                exptime = header['EXPTIME']
                zpt = header['PHOTZP']
                pixscale = header['PIXSCAL1']
            except:
                continue

    imgdim = 501        # Should be odd number so that we have a central pixel

    halfdim = (imgdim-1)/2

    emptyimg = (np.random.poisson(1000, (imgdim, imgdim)) - 1000).astype('float64')    # Create empty array with "noise", place psf in center of larger "image" for giving to trippy

    leftsize = (psfdata.shape[0]-1)/2     # Pixels left of center
    rightsize = (psfdata.shape[0]+1)/2    # Pixels right of center

    emptyimg[halfdim-leftsize : halfdim+rightsize, halfdim-leftsize : halfdim+rightsize] += psfdata     # Plant PSF from seepsf in center of empty image

    fwhm = 2*math.sqrt(2*math.log(2))*par1      # Get full width/half max using par1 of Gaussian PSF

    rate = math.sqrt(ra_rate**2+dec_rate**2)
    angle = math.degrees(math.atan2(dec_rate, ra_rate))      # Angle w.r.t. horizontal, between +/-90
    
    if angle>90:        # Since angle can only be in range +/-90
        angle -= 180

    if angle <-90:
        angle += 180

    if dec_rate<0:
        angle=-1*angle

    # the center of the source is the center of the image
    seecen = np.array([emptyimg.shape[1]/2.+0.5])

    # Generate model of size 61x61, alpha and beta are moffat fit values, 1.75 determined experimentally when matching moffat to Gaussian shape
    goodPSF=psf.modelPSF(np.arange(61), np.arange(61), alpha=1.75*par1, beta=2)
    goodPSF.genLookupTable(emptyimg, seecen, seecen)    # Generate lookup table for model PSF, centered at the center of the psf image
    
    '''z2 = goodPSF.lookupTable.max()
    z1 = goodPSF.lookupTable.min()
    normer=interval.ManualInterval(z1,z2)
    pyl.imshow(normer(goodPSF.lookupTable))
    pyl.show()'''

    goodPSF.line(rate, angle, exptime/3600., pixScale=pixscale, useLookupTable=True)    # Generate trailed/line PSF

    # Use if you want to save image of the trippy psf model
    #modelimg = psf_image.replace('seepsf','model')
    #goodPSF.psfStore(modelimg)

    stampdata = stampdata.astype('float64')     # Stamp data is originally ints
    bgf=bgFinder.bgFinder(stampdata)            # Calculate and subtract background
    bg=bgf.smartBackground()
    stampdata-=bg

    centroid = sep_phot(stampdata, 0.002, 3)        # Use Kristi's code to get centroid info using SEP

    r_close = 1000                                        # Arbitrary number that will intentionally be too large 

    for i in range(len(centroid['x'])):             # For all the objects identified in the stamp by SEP, find the one closest to the center of the image
        dist = (centroid['x'][i]-stampdata.shape[0]/2)**2 + (centroid['y'][i]-stampdata.shape[1]/2)**2
        if dist < r_close:
            r_close = dist
            closest = i
        
    xcen, ycen = centroid['x'][closest], centroid['y'][closest]     # Use the centroid that was closest to the center, since it will correspond to the asteroid

    print 'Using centroid {} {} (in numpy coordinates)'.format(xcen,ycen)

    phot=pill.pillPhot(stampdata)       # Perform trailed source photometry, calculate SNR, total flux
    phot(xcen, ycen, radius=1.1*fwhm, l=(exptime/3600.)*rate/pixscale, a=angle, skyRadius=4*fwhm, width=6*fwhm, zpt=zpt, exptime=exptime)
    phot.SNR(verbose=True)

    fitter = MCMCfit.MCMCfitter(goodPSF, stampdata)     # Fit trailed model to background subtracted stamp image
    fitter.fitWithModelPSF(xcen, ycen, m_in=-1, bg=bg, useLinePSF=True,fitWidth=15,nWalkers=20,nStep=20,nBurn=20,useErrorMap=True)
    (fitPars,fitRange)=fitter.fitResults(0.67)      # 0.67 = fit to 1 sigma
    print '1-sigma range on Best Point', fitRange

    modelImage=goodPSF.plant(fitPars[0],fitPars[1],fitPars[2],stampdata,addNoise=False,useLinePSF=True,returnModel=True)        # Model the trailed source
    
    pyl.imshow(modelImage)
    pyl.show()

    removed=goodPSF.remove(fitPars[0],fitPars[1],fitPars[2],stampdata,useLinePSF=True)      # Subtract the model from the stamp

    HDU=fits.PrimaryHDU(removed)
    List=fits.HDUList([HDU])
    psf_removed = psf_image.replace('seepsf','removed')
    os.unlink(psf_image)
    List.writeto(psf_removed,clobber=True)

    print 'Asteroid PSF source flux: {}'.format(np.sum(modelImage))

    z1 = stampdata[stampdata.shape[0]/2-30:stampdata.shape[0]/2 + 30,stampdata.shape[1]/2-30:stampdata.shape[1]/2 + 30].min()
    z2 = stampdata[stampdata.shape[0]/2-30:stampdata.shape[0]/2 + 30,stampdata.shape[1]/2-30:stampdata.shape[1]/2 + 30].max()
    normer=interval.ManualInterval(z1,z2)
    pyl.imshow(normer(stampdata))
    pyl.show()
    pyl.imshow(normer(removed))
    pyl.show()
