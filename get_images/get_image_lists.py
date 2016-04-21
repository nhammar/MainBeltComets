import os
from astropy.io import ascii
from astropy.time import Time
from astropy import units
import requests
import argparse
import sys
import time
import find_family
from ossos_scripts.ssois_query import Query
import config
from ossos import storage, horizons

_DIR_PATH_BASE = config._DIR_PATH_BASE
_FAMILY_LISTS = config._FAMILY_LISTS
_OUTPUT_DIR = config._OUTPUT_DIR 

def main():
    """
    Input asteroid family, filter type, and image type to query SSOIS
    """

    parser = argparse.ArgumentParser(description='Run SSOIS and return the available images in a particular filter.')

    parser.add_argument("--filter",
                        action="store",
                        default='r',
                        dest="filter",
                        choices=['r', 'u'],
                        help="Passband: default is r.")
    parser.add_argument("--family", '-f',
                        action="store",
                        default=None,
                        help='List of objects to query.')
    parser.add_argument("--member", '-m',
                        action="store",
                        default=None,
                        help='Member object of family to query.')

    args = parser.parse_args()

    if args.family != None and args.member == None:
        get_family_info(str(args.family), args.filter)
    elif args.family == None and args.member != None:
        get_member_info(str(args.member), args.filter)
    else:
        print "Please input either a family or single member name"


def get_family_info(familyname, filtertype='r', imagetype='p'):
    """
    Query the ssois ephemeris for images of objects in a given family. Then parse through for desired image type, 
    filter, exposure time, and telescope instrument
    """

    # establish input
    family_list = '{}/{}_family.txt'.format(_FAMILY_LISTS, familyname)

    if os.path.exists(family_list):
        with open(family_list) as infile:
            filestr = infile.read()
        object_list = filestr.split('\n')  # array of objects to query
    elif familyname == 'all':
        object_list = find_family.get_all_families_list()
    else:
        object_list = find_family.find_family_members(familyname)

    for obj in object_list[0:len(object_list) - 1]:  # skip header lines
        get_member_info(obj, filtertype)


def get_member_info(object_name, filtertype='r', imagetype='p'):
    """
    Query the ssois ephemeris for images of a given object. Then parse through for desired image type, 
    filter, exposure time, and telescope instrument
    """

    # From the given input, identify the desired filter and rename appropriately                    Replace this?
    if filtertype.lower().__contains__('r'):
        filtertype = 'r.MP9601'  # this is the old (standard) r filter for MegaCam
    if filtertype.lower().__contains__('u'):
        filtertype = 'u.MP9301'

    # Define time period of image search, basically while MegaCam in operation
    search_start_date = Time('2013-01-01', scale='utc')  # epoch1=2013+01+01
    search_end_date = Time('2017-01-01', scale='utc')  # epoch2=2017+1+1

    print "----- Searching for images of object {}".format(object_name)

    query = Query(object_name, search_start_date=search_start_date, search_end_date=search_end_date)
    try:
        objects = parse_ssois_return(query.get(), object_name, imagetype, camera_filter=filtertype)
    except IOError:
        print "Sleeping 30 seconds"
        time.sleep(30)
        objects = parse_ssois_return(query.get(), object_name, imagetype, camera_filter=filtertype)

    # Setup output, label columns
    if len(objects)>0:
        output = '{}/{}_object_images.txt'.format(_OUTPUT_DIR, object_name)
        with open(output, 'w') as outfile:
            outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                "Object", "Image", "Exp_time", "RA (deg)", "Dec (deg)", "time", "filter", "RA rate (\"/hr)", "Dec rate (\"/hr)"))

        mjds = []
        for line in objects:
            mjds.append(float(line['MJD'])) #Have to convert elements to floats
        start_time = Time(min(mjds), format='mjd') - 1.0*units.minute
        stop_time = Time(max(mjds), format='mjd') + 1.0*units.minute

        #Query Horizons once to establish position values over given time period, then give it a current time which it interpolates withl
        body = horizons.Body(object_name, start_time=start_time, stop_time=stop_time, step_size=10 * units.minute)

        for line in objects:
            with open(output, 'a') as outfile:
                time = Time(line['MJD'], format='mjd', scale='utc')
                time.format = 'iso'
                body.current_time = time
                p_ra = body.coordinate.ra.degree  
                p_dec = body.coordinate.dec.degree 
                ra_dot = body.ra_rate.to(units.arcsecond/units.hour).value
                dec_dot = body.dec_rate.to(units.arcsecond/units.hour).value
                try:
                    outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        object_name, line['Image'], line['Exptime'], p_ra, p_dec,
                        Time(line['MJD'], format='mjd', scale='utc'), line['Filter'], ra_dot, dec_dot))
                except Exception, e:
                    print "Error writing to outfile: {}".format(e)


def parse_ssois_return(ssois_return, object_name, imagetype, camera_filter='r.MP9601',
                       telescope_instrument='CFHT/MegaCam'):
    """
    Parse through objects in ssois query and filter out images of desired filter, type, exposure time, and instrument
    """

    assert camera_filter in ['r.MP9601', 'u.MP9301']

    ret_table = []
    good_table = 0

    table_reader = ascii.get_reader(Reader=ascii.Basic)
    table_reader.inconsistent_handler = _skip_missing_data
    table_reader.header.splitter.delimiter = '\t'
    table_reader.data.splitter.delimiter = '\t'
    table = table_reader.read(ssois_return)

    for row in table:
        # Excludes the OSSOS wallpaper.
        # note: 'Telescope_Insturment' is a typo in SSOIS's return format
        if not 'MegaCam' in row['Telescope_Insturment']:
            continue
        if not storage.exists(storage.get_uri(row['Image'][:-1])):  #Check if image of object exists in OSSOS observations
            continue
        if not str(row['Image_target']).startswith('WP'):
           good_table += 1
           ret_table.append(row)

    if good_table > 0:
        print " %d images found" % good_table

    return ret_table


def _skip_missing_data(str_vals, ncols):
    """
    add a extra column if one is missing, else return None.
    """
    if len(str_vals) == ncols - 1:
        str_vals.append('None')
        return str_vals
    # else:
        # raise ValueError("not enough columns in table")
        # print '  Not enough columns in data table'


if __name__ == '__main__':
    main()
