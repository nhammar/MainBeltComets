# From a list of identified MBCs, query and parse the SSOIS to get images from CFHT/MegaCam in r and u filters

import datetime
import os
import warnings
from astropy.io import ascii
from astropy.time import Time
import requests
import sys
import argparse
import requests
from requested_band_images import parse_ssois_return


# *** input is a line separated list of object names or reference numbers (dep. on ephemeris)

# IDENTIFY PARAMETERS FOR QUERY OF SSOIS FROM INPUT

# From the given input, make list of MBCs to query
#   ephemeris, name, date range, resolve image extension, resolve to x,y, positional uncertainty?
    ''' http:// www3.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/cadcbin/ssos/ssosclf.pl?lang=en; object=elst-pizarro; search=bynameMPC;
    epoch1=2013+01+01; epoch2=2015+1+16; eellipse=; eunits=arcseconds; extres=yes; xyres=yes; format=tsv '''



parser = argparse.ArgumentParser(description='Run SSOIS on a given set of MBCs and return the available images in a particular filter.')

parser.add_argument("--filter", "-f",
                    action="store",
                    default='r',
                    dest="filter",
                    choices=['r', 'u'],
                    help="passband: default is r'")
parser.add_argument("--listin",
                        action="store",
                        default="mbc.txt",
                        help='vospace dbaseclone containerNode')
parser.add_argument("--dbimages",
                        action="store",
                        default="vos:OSSOS/dbimages",
                        help='vospace dbimages containerNode')
parser.add_argument('--type',
                        default='p',
                        choices=['o', 'p', 's'], 
                        help="restrict type of image (unprocessed, reduced, calibrated)")
parser.add_argument("--output", "-o",
                        action="store",
                        default="/Users/admin/Desktop/band.txt",   
                        help='Location and name of output file containing image IDs.')

args = parser.parse_args()

with open(mbc_file) as infile: 
    filestr = infile.read()
input_mbc_lines = filestr.split('\n') # array of MBCs to query

    ''' little johnny 133P doesnt like the fame that comes with being an identified comet, 
    and wants to be left alone. Accordingly, you won't find and images of him from CFHT in 
    the MPC ephemeris. You can sometimes find him under to pen-name elst-pizarro though '''

    ''' test the code for one object at first, need for loop to cycle through input_mbc_lines later '''
    
object_name = 'elst-pizarro'

# Confirm that the input is the proper format to search for the appropriate ephemeris

# print input_mbc_lines # more elegant check to be done later

# From the given input, identify the desired filter and rename appropriately

if args.filter.lower().__contains__('r'):
    args.filter = 'r.MP9601'  # this is the old (standard) r filter for MegaCam
if args.filter.lower().__contains__('u'):
    args.filter = 'u.MP9301'

# Define time period of image search, basically while MegaCam in operation

search_start_date=Time('2013-01-01', scale='utc')   # epoch1=2013+01+01
search_end_date=Time('2017-01-01', scale='utc')     # epoch2=2017+1+1

# Define URL to search the SSOIS from

SSOS_URL = "http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/cadcbin/ssos/ssosclf.pl"    # no clf?

# Query a specific ephemeris depending on object naming convention. CADC : search = bynameCADC , MPC : search = bynameMPC
# Specify search parameters

resolve_extension = yes
resolve_position = yes
RESPONSE_FORMAT = 'tsv'

    ''' refer to ssos ParamDictBuilder and Query '''

params = dict(format = RESPONSE_FORMAT, epoch1 = str(search_start_date), epoch2 = str(search_end_date), search = bynameMPC, 
                    eunits = arcseconds, extres = resolve_extension, xyres = resolve_position, object = object_name )
logger.debug("{}\n".format(params))

# QUERY THE SSOIS DATABASE WITH ^ PARAMETERS

    ''' differnce between request.get() and request.post() 
        GET - Requests data from a specified resource
        POST - Submits data to be processed to a specified resource '''
        
headers = {'User-Agent': 'OSSOS Target Images'} # still thinking about what this is for....

response = requests.post(SSOS_URL, data=params, headers=headers)

# Confirm that the results are of the appropriate format, else return error statements

''' copied straight from ssos Query.get() '''

logger.info(response.url)
assert isinstance(self.response, requests.Response)

assert (self.response.status_code == requests.codes.ok )

lines = self.response.content
# note: spelling 'occured' is in SSOIS
if len(lines) < 2 or "An error occured getting the ephemeris" in lines:
    raise IOError(os.errno.EACCES, "call to SSOIS failed on format error")

if os.access("backdoor.tsv", os.R_OK):
    lines += open("backdoor.tsv").read()


# GET/REVIEW THE DATA RETURNED FROM THE SEARCH
    
obs_in_filter = requested_band_images.parse_ssois_return(lines, camera_filter=args.filter)


# Parse the data for Image, Instrument, Filter and create table for each object
# Download files with appropriate names? ie object+SSOISfilename ??

if len(obs_in_filter) > 0:
    with open(args.output, 'a') as outfile:
        for line in obs_in_filter:
            outfile.write("{:>10s} {:>10s} {:>2s} {:6.1f} {:6.1f} {:8.16f} {:8.16f} {:>10s} {:>10s}\n".format(orbit.name,
                            line['Image'], str(line['Ext']), line['X'], line['Y'], line['Object_RA'], line['Object_Dec'],
                            Time(line['MJD'], format='mjd', scale='utc'), line['Filter']))
        outfile.write('\n')

# Confirm that there is coordinates for each image
# Add information to table? or file corresponding to the image name?



