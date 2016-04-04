import math
from trippy import psf, pill, MCMCfit
from astropy.io import fits
import numpy as np

def make_trippy_profile(objnum, expnum, exptime, zpt, pixscale, ra_rate, dec_rate, amplitude, psfrad, psf_image):

    with fits.open(psf_image) as seepsf:
        for ext in range(len(seepsf)):
            try:
                data = seepsf[ext].data
                header = seepsf[ext].header
                par1 = header['PAR1']
            except:
                continue

    fwhm = 2*math.sqrt(2*math.log(2))*par1

    rate = math.sqrt(ra_rate**2+dec_rate**2)
    angle = math.degrees(math.atan2(dec_rate, math.fabs(ra_rate)))

    # Check that the array sizes and X/Y location of the PSF are odd numbers
    # X and Y are reversed in Numpy

    print data.shape[1], data.shape[0]

    if data.shape[1]%2==0 :
        xdim = data.shape[1] + 1 
    else:
        xdim = data.shape[1] 

    if data.shape[0]%2==0 :
        ydim  = data.shape[0] + 1 
    else:
        ydim = data.shape[0] 

    print xdim, ydim, data.shape

    # the center of the PSF is the center of the seepsf output.
    xcen = np.array([data.shape[1]/2.])
    ycen = np.array([data.shape[0]/2.])
    print xcen, ycen

    goodPSF=psf.modelPSF(np.arange(xdim), np.arange(ydim), alpha=1.75*par1, beta=2, repFact=10)
    goodPSF.genLookupTable(data, xcen, ycen)
    goodPSF.line(rate, angle, exptime/3600., pixScale=pixscale, useLookupTable=True)

    goodPSF.computeRoundAperCorrFromPSF(np.linspace(1.1*fwhm,3.9*fwhm,10))
    roundAperCorr=goodPSF.roundAperCorr(1.1*fwhm)
    goodPSF.computeLineAperCorrFromTSF(np.linspace(1.1*fwhm,3.9*fwhm,10), l=(exptime/3600.)*rate/pixscale, a=angle)
    lineAperCorr=goodPSF.lineAperCorr(1.1*fwhm)

    goodPSF.psfStore('psf.fits')
    phot=pill.pillPhot(data,repfact=10)

    phot(xcen, ycen, radius=1.1*fwhm, l=(exptime/3600.)*rate/pixscale, a=angle, skyRadius=4*fwhm, width=6*fwhm, zpt=zpt, exptime=exptime)
    phot.SNR()
    phot.computeRoundAperCorrFromSource(xcen, ycen, np.linspace(1.1*fwhm,3.9*fwhm,10), 4*fwhm, 6*fwhm)

    Data = stamp - phot.bg

    fitter = MCMCfit.MCMCfitter(goodPSF, Data)
    fitter.fitWithModelPSF(xcen, ycen, m_in=amplitude, bg=phot.bg, useLinePSF=True)
    (fitPars,fitRange)=fitter.fitResults(0.67)
    modelImage=goodPSF.plant(fitPars[0],fitPars[1],fitPars[2],Data,useLinePSF=True)
    removed=goodPSF.remove(fisPars[0],fitPars[1],fitPars[2],Data,useLinePSF=True)
