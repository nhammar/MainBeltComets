import sys
import argparse
from sys import argv

import datetime
from astropy.time import Time
from astropy.io import ascii

parser = argparse.ArgumentParser(
        description='Run SSOIS on a given set of minor planets and return the available images in a particular filter.')

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

infile = open(args.listin)

if args.filter.lower().__contains__('r'):
    args.filter = 'r.MP9601'  # this is the old (standard) r filter for MegaCam
if args.filter.lower().__contains__('u'):
    args.filter = 'u.MP9301'
    
search_start_date=Time('2013-01-01', scale='utc')
search_end_date=Time('2017-01-01', scale='utc'))


