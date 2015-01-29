import os
import sep
import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.io import ascii
import argparse
from drizzlepac import pixtosky


# Identify objects in each postage stamp


def main(): 
    
    parser = argparse.ArgumentParser(
        description='For an input .fits image, aperture size, threshold, and output file: preforms photometry')
    parser.add_argument("--ossin",
                        action="store",
                        default="3330_stamps_gt8/",
                        help="The directory in getImages/3330/ with input .fits files for astrometry/photometry measurements.")
# CREATE A DIRECTORY FOR OUTPUT?
#    parser.add_argument("--output", "-o",
#                        action="store",
#                        default="/Users/admin/Desktop/MainBeltComets/getImages/sep_phot.txt",   
#                        help='Location and name of output file containing image photometry values.')
    parser.add_argument("--radius", '-r',
                        action='store',
                        default=6.0,
                        help='aperture (degree) of circle for photometry.')
    parser.add_argument("--thresh", '-t',
                            action='store',
                            default=1.5,
                            help='threshold value.')
                            
# PARSE IN MEASURED X AND Y COORDINATES
    
    
    args = parser.parse_args()
    global ap
    ap = float(args.radius)
    global th
    th = float(args.thresh)
    
    global x
    x = 387
    global y
    y = 422
    
    for image in os.listdir('/Users/admin/Desktop/MainBeltComets/getImages/{}'.format(args.ossin)):
        if image.endswith('.fits') == True:
            with fits.open('{}/{}'.format(args.ossin, image)) as hdulist:
                # doesnt work for 304757_1692837p_242.335075_-14.024202.fits specifically for default values
                print "Doing photometry on image %s " % image
                #print hdulist.info()
                if (hdulist[0].data == None):
                    table1 = dosep(hdulist[1].data)
                    table2 = dosep(hdulist[2].data)
                    table = vstack([table1, table2])
                    ascii.write(table, '{}_info.txt'.format(image))
                    compare(table, image)
                    # what if more than two ccds???
                else:
                    table0 = dosep(hdulist[0].data)
                    ascii.write(table0, '{}_info.txt'.format(image))
                    compare(table, image)
                #print data
        
def dosep(data):
        
    # Measure a spatially variable background of some image data (numpy array)
    bkg = sep.Background(data) #, mask=mask, bw=64, bh=64, fw=3, fh=3) # optional parameters
        
    # Evaluate the spatially variable background and RMS:
    back = bkg.back() # creates an array, same shape and type as data
    rms = bkg.rms()   # creates an array, same shape and type as data
        
    # Directly subtract the background from the data in place
    bkg.subfrom(data)
    bkg.globalback    # Global "average" background level
    bkg.globalrms     # Global "average" RMS of background
        
    # for the background subtracted data, detect objects in data given some threshold
    # ****CHOOSE APPROPRIATE THRESHOLD****
    thresh = th * bkg.globalrms    # ensure the threshold is high enough wrt background        
    objs = sep.extract(data, thresh)
    #print len(objs)
    #print objs['x'][0] # print flux-wieghted 

    # calculate the Kron radius for each object, then we perform elliptical aperture photometry within that radius
    kronrad, krflag = sep.kron_radius(data, objs['x'], objs['x'], objs['y'], objs['b'], objs['theta'], ap)
    flux, fluxerr, flag = sep.sum_ellipse(data, objs['x'], objs['x'], objs['y'], objs['b'], objs['theta'], 2.5*kronrad, subpix=1)
    flag |= krflag  # combine flags into 'flag'
    
    # mask = np.zeros(data.shape, dtype=np.bool)
    # sep.mask_ellipse(mask, objs['x'], objs['y'], obs['a'], objs['b'], objs['theta'], r=3.)

    # Specify a per-pixel "background" error and a gain. This is suitable when the data have been background subtracted.
    # ***WHAT IS THE IMAGE GAIN? 1.67 ?***
    # flux, fluxerr, flag = sep.sum_circle(data, objs['x'], objs['y'], 3.0, err=bkg.globalrms, gain=1.0)

    # write to ascii table
    table = Table([objs['x'], objs['y'], flux], names=('x', 'y', 'flux'))
    return table

def compare(table, image):
    # compare predicted RA and DEC to that measured by sep photometry
    # get predicted RA and DEC from text output from getImages
    with open('test_images.txt') as infile:
        for line in infile.readlines()[1:]:
            assert len(line.split()) > 0
            objectname = line.split()[0]
            image = line.split()[1]
            pRA = float(line.split()[3])
            pDEC = float(line.split()[4])
            print " Predicted RA and DEC for object {} in image {}: {}  {}".format(objectname, image, pRA, pDEC)
    
    # parse through table and get RA and DEC closest to measured-by-eye coordinates
    # compare to predicted
    for row in table:
        x_max = x + 2
        x_min = x - 2
        y_max = y + 2
        y_min = y - 2
        if (float(row['x']) < x_max) & (float(row['x']) > x_min) & (float(row['y']) < y_max) & (float(row['y']) > y_min):
            mRA_pix = float(row['x'])
            mDEC_pix = float(roy['y'])
    mRA, mDEC = pixtosky.xy2rd('{}'.format(image), mRA_pix, mDEC_pix)
    print " Measured RA and DEC for object {} in image {}: {}  {}".format(objectname, image, mRA, mDEC)
        

if __name__ == '__main__':
    main()
