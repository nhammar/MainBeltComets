from astropy.coordinates import SkyCoord
from ossos import horizons
from astropy import units

def get_coords(object_name, mid_time):
    """
    Queries the JPL Horizon's ephemeris for rate of change of RA and DEC for a specific day
    """

    start_time = mid_time - 1.0*units.minute
    stop_time = mid_time + 1.0*units.minute

    step_size = 1.0*units.minute

    body = horizons.Body(object_name, start_time, stop_time, step_size)

    body.current_time = mid_time

    ra = body.coordinate.ra.degree
    dec = body.coordinate.dec.degree
    ra_dot_cos_dec = body.ra_rate
    dec_dot = body.dec_rate


    return ra, dec, ra_dot_cos_dec, dec_dot
