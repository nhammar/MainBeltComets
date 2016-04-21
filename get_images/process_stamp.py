import argparse
import make_iraf_psf as mkpsf
import make_trippy_profile as mktrp
import glob
import config

_STAMPS_DIR = config._STAMPS_DIR

def main():

    parser = argparse.ArgumentParser(description='Produce IRAF PSF, give it to Trippy.')

    parser.add_argument("--stamp", '-s',
                        action="store",
                        default="all",
                        help='Path to stamp (usually postage_stamps/{name}.fits) to process')

    args = parser.parse_args()

    if args.stamp == 'all':
        for fname in glob.glob("{}/*.fits".format(_STAMPS_DIR)):
            process(fname)
    else:
        process(args.stamp)


def process(stamp):

    psf_image = mkpsf.make_iraf_psf(stamp)
    if psf_image != None:        
        mktrp.make_trippy_profile(psf_image, stamp)
    else:
        pass


if __name__ == '__main__':
    main()
