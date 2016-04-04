import argparse
from cStringIO import StringIO
import getpass
from astropy.io import fits
import requests
import os
import vos
from astropy.time import Time
import urllib2 as url
import numpy as np
from astropy.table import Table, Column
import pandas as pd
import sys
from ossos import storage
from ossos import coding
from ossos import mpc
from ossos import util
from astropy.coordinates import SkyCoord
from astropy import units
import config

_TARGET = "TARGET"
_DIR_PATH_BASE = config._DIR_PATH_BASE
_IMAGE_LISTS = config._IMAGE_LISTS
_STAMPS_DIR = config._STAMPS_DIR
_VOS_PATH = 'vos:OSSOS/ActiveAsteroids/postage_stamps'

BASEURL = "http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/vospace/auth/synctrans"

_INPUT = 'object_images.txt'

"""
Retrieval of cutouts of the FITS images associated with the CFHT/MegaCam detections.
Takes a table (get_images.py output) as input
An Example URL for cutouts from OSSOS (not CFHT/MegaCam)
http://www.canfar.phys.uvic.ca/vospace/auth/synctrans?
TARGET=vos://cadc.nrc.ca~vospace/OSSOS/dbimages/1625356/1625356p.fits&
DIRECTION=pullFromVoSpace&
PROTOCOL=ivo://ivoa.net/vospace/core%23httpget&
view=cutout&
cutout=CIRCLE+ICRS+242.1318+-12.4747+0.05
"""


def main():
    parser = argparse.ArgumentParser(
        description='Parse an familyname.images.txt input file (get_images.py output) and create links in the postage stamp directory '
                    'that allow retrieval of cutouts of the FITS images associated with the CHT/MegaCam detections. '
                    'Cutouts are defined on the WCS ra/dec of the object position.')
    parser.add_argument("--member", '-m',
                        action="store",
                        default=None,
                        help='Single object number for input .txt files.')
    parser.add_argument("--radius", '-r',
                        action='store',
                        default=0.005,
                        help='Radius (degrees) of cutout postage stamp.')

    args = parser.parse_args()

    os.path.exists(_STAMPS_DIR) or os.makedirs(_STAMPS_DIR)

    if args.member != None:
        # CADC PERMISSIONS
        username = raw_input("CADC username: ")
        password = getpass.getpass("CADC password: ")
        get_stamps(str(args.member), username, password, args.radius)
    else:
        print "Please supply an object name"

def get_stamps(member, username, password, radius):
    """
    Cutout postage stamps of objects from VOSpace OSSOS mosaic, upload to VOSpace again
    """

    print "----- Cutting postage stamps of object {} -----".format(member)

    image_list_path = '{}/{}_{}'.format(_IMAGE_LISTS, member, _INPUT)

    try:
        table = pd.read_table(image_list_path, usecols=[0, 1, 3, 4, 7, 8], header=0, names=['object', 'expnum', 'ra', 'dec', 'rarate', 'decrate'],
                              sep='\t', dtype={'Object': object, 'Image': object})
    except pd.parser.CParserError:
        table = pd.read_table(image_list_path, usecols=[0, 1, 3, 4, 7, 8], header=0, names=['object', 'expnum', 'ra', 'dec', 'rarate', 'decrate'],
                              sep=' ', dtype={'Object': object, 'Image': object})

    for row in range(len(table)):
        postage_stamp_filename = "{}_{}_{:8f}_{:8f}.fits".format(table['object'][row],
                                                                 table['expnum'][row],
                                                                 table['ra'][row],
                                                                 table['dec'][row])
        cutout(member, table['object'][row], table['expnum'][row], table['ra'][row],
                   table['dec'][row], radius, table['rarate'][row], table['decrate'][row])


def cutout(member, object_name, image, ra, dec, radius, ra_rate, dec_rate):

    global vos_dir

    if member == 'none':
        vos_dir = '{}/nofamily'.format(_VOS_PATH)   #Planned to sort stamps on VOS into groups, but this didn't happen and so currently everything will be saved to /hasfamily
    else:
        vos_dir = '{}/hasfamily'.format(_VOS_PATH)

    radius = radius * units.deg

    expnum = image.strip('p')  # only use processed images
    target = storage.get_uri(expnum)
    coord = SkyCoord(unit="deg", ra=ra, dec=dec)

    hdulist = storage.ra_dec_cutout(target, coord, radius)

    hdulist[0].header['OBJNUM'] = (object_name,'object')    #Put image info into header for later use
    hdulist[0].header['RA'] = (ra,'degrees')
    hdulist[0].header['DEC'] = (dec,'degrees')
    hdulist[0].header['RARATE'] = (ra_rate,'arcsec/hr')
    hdulist[0].header['DECRATE'] = (dec_rate,'arcsec/hr')

    postage_stamp_filename = "{}_{}_{:8f}_{:8f}.fits".format(object_name, image, float(ra), float(dec))

    print postage_stamp_filename

    hdulist.writeto("{}/{}".format(_STAMPS_DIR, postage_stamp_filename), output_verify='warn', clobber=True)
    del hdulist

    try:
        storage.copy('{}/{}'.format(_STAMPS_DIR, postage_stamp_filename),
                     '{}/{}'.format(vos_dir, postage_stamp_filename))
    except Exception, e:
        print e


def write_to_file(object_name, expnum):

    with open('{}/cutout_error.txt'.format(_IMAGE_LISTS), 'a') as outfile:
        outfile.write('{} {}\n'.format(object_name, expnum))

if __name__ == '__main__':
    main()
