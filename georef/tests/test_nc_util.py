"""
test_nc_util.py

Routines to test georef's nc_info.py code
"""

import os
import georef
import georef.nc_util


def test_write_sample_netCDF_file():
    """Test that can write a simple geo-referenced netCDF file"""
    test_filename = './test_georef.nc'
    try:
        assert not os.path.isfile(test_filename)
    except AssertionError:
        print('test file exists, skipping test: {}'.format(test_filename))
        return
    georef.nc_util.write_simple_netCDF_file(test_filename, overwrite=True)
    os.remove(test_filename)
