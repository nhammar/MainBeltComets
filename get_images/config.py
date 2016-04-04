import os
import shutil
import sys

_DIR_PATH_BASE = os.path.dirname(os.path.abspath(__file__))	# Path base

_FAMILY_LISTS = '{}/family_lists'.format(_DIR_PATH_BASE)	# Directory to write family lists to
_OUTPUT_DIR = '{}/image_lists'.format(_DIR_PATH_BASE)       # Directory to write image lists to

_IMAGE_LISTS = _OUTPUT_DIR	                                # Directory where image lists are written to, must match output_dir above
_STAMPS_DIR = '{}/postage_stamps'.format(_DIR_PATH_BASE)    # Directory to write postage stamps to

_CONFIG_DIR = '.AA_config'                                  # Name for config directory

user_config_path = os.path.join(os.path.expanduser('~'), _CONFIG_DIR)
user_config = os.path.join(user_config_path, 'config.py')   # Create config in user's home directory

if os.path.exists(user_config):                             # If a config exists in the user's home directory already, import it
    if os.path.abspath(__file__) != user_config_path:
        sys.path.append(user_config)
        import config
else:
    os.makedirs(user_config)                                # otherwise use this one but create a copy in the user's home directory
    shutil.copy(os.path.basename(__file__), user_config_path)
