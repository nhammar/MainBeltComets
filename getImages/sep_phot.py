import os
import sep
import urllib2 as url
import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.io import ascii
from astropy.time import Time
import argparse
from scipy.spatial import cKDTree
import math

import sys
sys.path.append('/Users/admin/Desktop/MainBeltComets/getImages/ossos_scripts/')

import ossos_scripts.wcs as wcs
from ossos_scripts.storage import get_astheader
from getImages import get_image_info

def main(): 
    ''' 
    Preforms photometry on .fits files given an input of family name and object name
    Identifies object in image from predicted coordinates, magnitude (and eventually shape)
    Assumes files organised as:
    dir_path_base/familyname/familyname_objectname/*.fits       - images to do photometry on
      or dir_path_base/familyname/familyname_stamps/*.fits
    dir_path_base/familyname/*_images.txt                       - list of image exposures, predicted RA and DEC, dates etc.
    '''
    
    parser = argparse.ArgumentParser(
                        description='For an object in an asteroid family, parses AstDys for a list of members, \
                        parses MPC for images of those objects from CFHT/MegaCam in specific filter and exposure time,\
                        cuts out postage stamps images of given radius (should eventually be uncertainty ellipse), \
                         preforms photometry on a specified object given an aperture size and threshold, \
                        and then selects the object in the image from the predicted coordinates, magnitude, and eventually shape')
    parser.add_argument("--family", '-f',
                        action="store",
                        default="3330",
                        help="Asteroid family name. Usually the asteroid number of the largest member.")
    parser.add_argument("--aperture", '-a',
                        action='store',
                        default=10.0,
                        help='aperture (degree) of circle for photometry.')
    parser.add_argument("--thresh", '-t',
                        action='store',
                        default=5.0,
                        help='threshold value for photometry (sigma above background).')
    parser.add_argument("--object", '-o',
                        action='store',
                        default=None,
                        help='The object to preform photometry on.')
    parser.add_argument("--filter",
                        action="store",
                        default='r',
                        dest="filter",
                        choices=['r', 'u'],
                        help="passband: default is r'")
    parser.add_argument('--type',
                        default='p',
                        choices=['o', 'p', 's'], 
                        help="restrict type of image (unprocessed, reduced, calibrated)")
                            
    args = parser.parse_args()
    
    find_objects_by_phot(args.family, args.object, float(args.aperture), float(args.thresh), args.filter, args.type)

def find_objects_by_phot(familyname, objectname, ap=10.0, th=5.0, filtertype='r', imagetype='p'):
    
    print objectname
    if objectname == None:
        image_list = get_image_info(familyname, filtertype, imagetype)
        objectname = image_list[0]
        print "WARNING: Object name not specified, searching for {}".format(objectname)
    print objectname
    
    # initiate directories
    family_dir, object_dir, output_dir = init_dirs(familyname, objectname)

    with open('{}/{}_r{}_t{}_output.txt'.format(object_dir, objectname, ap, th), 'w') as outfile:
        outfile.write("{:>3s} {:>5s} {:>17s} {:>10s} {:>17s} {:>10s} {:>17s}\n".format(
            "Image", "RA", "RA_pix", "DEC", "DEC_pix", "flux", "meas_mag"))        

    # get range of magnitudes of object over time span that the images were taken
    print "----- Querying JPL Horizon's ephemeris for apparent magnitudes -----"
    mag_list_jpl = mag_query_jpl(objectname, step=1)
    
    # preform the photometry and identify object
    print "----- Preforming photometry on all images of {} in family {} -----".format(objectname, familyname)
    for file in os.listdir('{}'.format(object_dir)):
        if file.endswith('.fits') == True:
            expnum_p = file.split('_')[1]

            with fits.open('{}/{}'.format(object_dir, file)) as hdulist:
                print " Preforming photometry on image {} ".format(expnum_p)
                
                if hdulist[0].data is None: # what if more than 2ccd mosaic? could just be aperture values?
                    try:
                        zeropt = fits.getval('{}/{}'.format(object_dir, file), 'PHOTZP', 1)
                        table1 = sep_phot(hdulist[1].data, ap, th)
                        table2 = sep_phot(hdulist[2].data, ap, th)
                        table = vstack([table1, table2])
                        #ascii.write(table, os.path.join(output_dir, '{}_phot.txt'.format(expnum_p))) # write all phot data to file in directory familyname/famlyname_objectname/sep_phot_output
                        astheader = hdulist[0].header
                    except LookupError: # maybe not correct error type?
                        print "ERROR: no PHOTZP in header, skipping image "
                    
                else:
                    try:
                        zeropt = fits.getval('{}/{}'.format(object_dir, file), 'PHOTZP', 0)
                        table = sep_phot(hdulist[0].data, ap, th)
                        astheader = hdulist[0].header
                        #ascii.write(table, os.path.join(output_dir, '{}_phot.txt'.format(expnum_p))) # write all phot data to file in directory familyname/famlyname_objectname/sep_phot_output
                    except LookupError:
                        print "ERROR: no PHOTZP in header, skipping image "
                
                object_data = comp_coords(table, expnum_p, astheader, zeropt, mag_list_jpl)
                assert len(object_data) > 0
                with open('{}/{}_r{}_t{}_output.txt'.format(object_dir, objectname, ap, th), 'a') as outfile:
                    try:
                        outfile.write('{} {} {} {} {} {} {}\n'.format(
                                object_data[0], object_data[1], object_data[2], object_data[3], 
                                object_data[4], object_data[5], object_data[6]))
                    except:
                        print "ERROR: cannot write to outfile"    

def mag_query_jpl(objectname, step=1, su='d'):
    '''
    Constructs a URL to query JPL Horizon's for apparent magnitude in a date range
    '''
    # from familyname_images.txt get date range of images for objectname
    date_range = []
    with open('{}/{}'.format(family_dir, imageinfo)) as infile:
        for line in infile.readlines()[1:]:
            assert len(line.split()) > 0
            if objectname == line.split()[0]:
                date_range.append(float(line.split()[5]))
    date_range_t = Time(date_range, format='mjd', scale='utc')
    assert  len(date_range_t.iso) > 0
    time_start = ((date_range_t.iso[0]).split())[0] + ' 00:00:00.0'
    time_end = ((date_range_t.iso[-1]).split())[0] + ' 00:00:00.0'
    
    if time_start == time_end:
        print "WARNING: only searching for one day"
        time_end_date = (((date_range_t.iso[-1]).split())[0]).split('-')
        day_add_one = int(time_end_date[2])+1
        if day_add_one < 10:
            day = '0{}'.format(day_add_one)
        else:
            day = day_add_one
        time_end = '{}-{}-{} 00:00:00.0'.format(time_end_date[0], time_end_date[1], day)
    
    print " Date range in query: {} -- {}".format(time_start, time_end)
    
    # change date format from 01-01-2001 00:00 to 01-Jan-2001 00:00
    date_start = change_date(time_start)
    date_end = change_date(time_end)
    
    s = '9' # select parameter for apparent magnitude
    
    # form URL pieces that Horizon needs for its processing instructions
    urlArr = ["http://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&COMMAND=",
              '',
              "&MAKE_EPHEM='YES'&TABLE_TYPE='OBSERVER'&START_TIME=",
              '',
              "&STOP_TIME=",
              '',
              "&STEP_SIZE=",
              '',
              "&QUANTITIES=" + s,
              "&CSV_FORMAT='YES'"]
              
    # change the object name, start and end times, and time step into proper url-formatting
    url_style_output = []
    for obj in [objectname, time_start, time_end]:
        os = obj.split()
        if len(os) > 1:
            ob = "'" + os[0] + '%20' + os[1] + "'"
        else:
            ob =  "'" + objectname + "'"
        url_style_output.append(ob)
    step = "'" + str(step) + "%20" + su + "'"
     
    # URL components
    urlArr[1] = url_style_output[0]  # formatted object name
    urlArr[3] = url_style_output[1]  # start time
    urlArr[5] = url_style_output[2]  # end time
    urlArr[7] = step  # timestep   
    urlStr = "".join(urlArr)  # create the url to pass to Horizons
       
    # Query Horizons; if it's busy, wait and try again in a minute
    done = 0
    while not done:
        urlHan = url.urlopen(urlStr)
        urlData = urlHan.readlines()
        urlHan.close()
        if len(urlData[0].split()) > 1:
            if "BUSY:" <> urlData[0].split()[1]:
                done = 1
            else:
                print urlData[0],
                print "Sleeping 60 s and trying again"
                time.sleep(60)
        else:
            done = 1   
    
    # parse through urlData for indexes of start and end dates    
    try: 
        for idx, line in enumerate(urlData):
            date_jpl = line.split()[0]+' '+(line.split()[1]).strip(',')
            if date_start == date_jpl:
                index_start = idx
    except:
        print "index start could not be obtained"
        index_start = 69
    try:
        for idx, line in enumerate(urlData):  #testing
            date_jpl = line.split()[0]+' '+(line.split()[1]).strip(',')
            if date_end == date_jpl:
                index_end = idx
    except:
        print "index end could not be obtained"
        index_end = 70
        
    # for indexes from start to end dates, get apparent magnitude values
    mag_list = []
    for line in urlData[index_start:index_end+1]:
        assert len(line.split()) > 0    
        try:
            mag_list.append(float((line.split()[4]).strip(',')))
        except:
            print "WARNING: could not get magnitude from JPL line"
    
    assert len(mag_list) > 0
    return mag_list
       
def change_date(date):
    '''
    Convert time format 01-01-2001 00:00:00.00 to 01-Jan-2001 00:00
    '''
    date_split = date.split('-')
    date_strip = (date_split[2]).split()
    month = int(date_split[1])
    month_name = ['nan', 'Jan', "Feb", 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    for i in range(0,13):
        if i == month:
            month_jpl = month_name[i]
    date_new = date_split[0]+'-'+month_jpl+'-'+date_strip[0]+' 00:00'
    
    return date_new    
    
def sep_phot(data, ap, th):
    ''' 
    Preforms photometry by SEP, similar to source extractor 
    '''

    # Measure a spatially variable background of some image data (numpy array)
    try:    
        bkg = sep.Background(data) #, mask=mask, bw=64, bh=64, fw=3, fh=3) # optional parameters
    except:
        data = data.byteswap(True).newbyteorder()
        bkg = sep.Background(data) #, mask=mask, bw=64, bh=64, fw=3, fh=3) # optional parameters
        
    # Evaluate the spatially variable background and RMS:
    back = bkg.back() # creates an array, same shape and type as data
    rms = bkg.rms()   # creates an array, same shape and type as data
        
    # Directly subtract the background from the data in place
    bkg.subfrom(data)
    bkg.globalback    # Global "average" background level
    bkg.globalrms     # Global "average" RMS of background
        
    # for the background subtracted data, detect objects in data given some threshold
    thresh = th * bkg.globalrms    # ensure the threshold is high enough wrt background        
    objs = sep.extract(data, thresh)

    # calculate the Kron radius for each object, then we perform elliptical aperture photometry within that radius
    kronrad, krflag = sep.kron_radius(data, objs['x'], objs['y'], objs['a'], objs['b'], objs['theta'], ap)
    flux, fluxerr, flag = sep.sum_ellipse(data, objs['x'], objs['y'], objs['a'], objs['b'], objs['theta'], 2.5*kronrad, subpix=1)
    flag |= krflag  # combine flags into 'flag'
   
    # write to ascii table
    table = Table([objs['x'], objs['y'], flux], names=('x', 'y', 'flux'))
    return table

def comp_coords(septable, expnum_p, astheader, zeropt, mag_list_jpl):
    '''
    Compares predicted RA and DEC to that measured by sep photometry
    Selects nearest neighbour object from predicted coordinates as object of interest
    Compares measured apparent magnitude to predicted
    '''

    tree = cKDTree(zip((np.array(septable['x'])).ravel(), (np.array(septable['y'])).ravel()))
    
    with open('{}/{}'.format(family_dir, imageinfo)) as infile:
        for line in infile.readlines()[1:]:
            assert len(line.split()) > 0
            expnum_p_fromfile = line.split()[1]
            
            pvwcs = wcs.WCS(astheader)
            
            # for entries in *_images.txt that correspond to images of the object
            if expnum_p_fromfile == expnum_p:
                objectname = line.split()[0]
                pRA = float(line.split()[3])
                pDEC = float(line.split()[4])
                expnum = (line.split()[1]).rstrip('p')
                
                pRA_pix, pDEC_pix = pvwcs.sky2xy(pRA, pDEC) # convert from WCS to pixels
                #print " Predicted RA and DEC: {}  {}".format(pRA, pDEC)
                #print "  in pixels: {} {}".format(pRA_pix, pDEC_pix)
                
                # parse through table and get RA and DEC closest to predicted coordinates (in pixels)
                
                coords = np.array([pRA_pix, pDEC_pix])
                d_list, i_list = tree.query(coords, k=10)

                mRA_pix = None
                for i in i_list:
                    flux = septable[i][2]
                    try:
                        mag_sep = -2.5*math.log10(flux)+zeropt
                        mean = np.mean(mag_list_jpl)
                        maxmag = np.amax(mag_list_jpl)
                        minmag = np.amin(mag_list_jpl)
                        
                        if ( 1 > maxmag - minmag):
                            if (abs(mag_sep - mean) < 1):
                                mRA_pix = septable[i][0]
                                mDEC_pix = septable[i][1]
                                break
                        else:
                            if (abs(mag_sep - mean) < maxmag - minmag):
                                mRA_pix = septable[i][0]
                                mDEC_pix = septable[i][1]  
                                break 
                    except:
                        None
                
                if mRA_pix == None:
                    print "WARNING: Magnitude condition could not be satisfied"
                    break
                
                #print "   Flux, mag: {}, {}".format(flux, mag_sep)     
                mRA, mDEC = pvwcs.xy2sky(mRA_pix, mDEC_pix) # convert from pixels to WCS
                #print " Measured RA and DEC: {}  {}".format(mRA, mDEC)
                #print "  in pixels: {} {}".format(mRA_pix, mDEC_pix)
                
                diffRA = mRA - pRA
                diffDEC = mDEC - pDEC
                #print " Difference: {} {}".format(diffRA, diffDEC)
                
                

                        
                return expnum_p, mRA, mRA_pix, mDEC, mDEC_pix, flux, mag_sep   

def init_dirs(familyname, objectname):
    global imageinfo
    imageinfo = familyname+'_images.txt'
    
    # would like dir_path_base to be specified
    dir_path_base = '/Users/admin/Desktop/MainBeltComets/getImages/'
    global family_dir
    family_dir = os.path.join(dir_path_base, familyname)
    if os.path.isdir(family_dir) == False:
        print "Invalid family name"
    stamps_dir = os.path.join(family_dir, familyname+'_stamps')
    object_dir = os.path.join(family_dir, familyname+'_'+objectname)
    if os.path.isdir(object_dir) == False:
        print "Invalid object name or object directory does not exist, searching through {}/{}_stamps instead".format(familyname, familyname)
        object_dir = stamps_dir
    output_dir = os.path.join(object_dir, 'sep_phot_output')
    if os.path.isdir(output_dir) == False:
        os.makedirs(output_dir)
    
    return family_dir, object_dir, output_dir
    
if __name__ == '__main__':
    main()
    