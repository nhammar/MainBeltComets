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

# *** input is a line separated list of object names or reference numbers (dep. on ephemeris)

# IDENTIFY PARAMETERS FOR QUERY OF SSOIS FROM INPUT

# From the given input, make list of MBCs to query

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
    the MPC ephemeris '''

    ''' test the code for one object at first, need for loop to cycle through input_mbc_lines later '''
payload = {'name': '176P'}

# Confirm that the input is the proper format to search for the appropriate ephemeris

print input_mbc_lines # more elegant check to be done later

# From the given input, identify the desired filter and rename appropriately

if args.filter.lower().__contains__('r'):
    args.filter = 'r.MP9601'  # this is the old (standard) r filter for MegaCam
if args.filter.lower().__contains__('u'):
    args.filter = 'u.MP9301'

# Define time period of image search, basically while MegaCam in operation

search_start_date=Time('2013-01-01', scale='utc')
search_end_date=Time('2017-01-01', scale='utc')

# Define URL to search the SSOIS from

SSOS_URL = "http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/cadcbin/ssos/ssosclf.pl"

# QUERY THE SSOIS DATABASE WITH ^ PARAMETERS

    ''' differnce between request.get() and request.post() 
        GET - Requests data from a specified resource
        POST - Submits data to be processed to a specified resource '''

query = request.post(SSOS_URL, data=payload, search_start_date=search_start_date, search_end_date=search_end_date)

# Query a specific ephemeris depending on object naming convention
    ''' not sure how to select a specific ephermeris a.t.m. '''
# Confirm that the results are of the appropriate format, else return error statements

# GET/REVIEW THE DATA RETURNED FROM THE SEARCH

    ''' how to I keep the search results? '''
    

query_content = StringIO(query.content)

# Parse the data for Image, Instrument, Filter and create table for each object
# Download files with appropriate names? ie object+SSOISfilename ??

# Confirm that there is coordinates for each image
# Add information to table? or file corresponding to the image name?



