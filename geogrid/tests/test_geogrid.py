"""
test_geogrid.py

Tests for geogrid's GeoGrid class

Note: some of these tests require a local netCDF file
"""

import geogrid
import numpy as np


epsg_3411_as_wkt = """
	PROJCS["NSIDC Sea Ice Polar Stereographic North",
    GEOGCS["Unspecified datum based upon the Hughes 1980 ellipsoid",
        DATUM["Not_specified_based_on_Hughes_1980_ellipsoid",
            SPHEROID["Hughes 1980",6378273,298.279411123064,
                AUTHORITY["EPSG","7058"]],
            AUTHORITY["EPSG","6054"]],
        PRIMEM["Greenwich",0,
            AUTHORITY["EPSG","8901"]],
        UNIT["degree",0.0174532925199433,
            AUTHORITY["EPSG","9122"]],
        AUTHORITY["EPSG","4054"]],
    PROJECTION["Polar_Stereographic"],
    PARAMETER["latitude_of_origin",70],
    PARAMETER["central_meridian",-45],
    PARAMETER["false_easting",0],
    PARAMETER["false_northing",0],
    UNIT["metre",1,
        AUTHORITY["EPSG","9001"]],
    AXIS["Easting",SOUTH],
    AXIS["Northing",SOUTH],
    AUTHORITY["EPSG","3411"]]
    """

epsg_3411_as_proj4 = '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +a=6378273 +b=6356889.449 +units=m +no_defs'


def test_can_declare_GeoGrid():
    """Test can initialize a GeoGrid from scratch"""
    xdim = 3
    ydim = 4
    data_values = np.arange(12.).reshape((ydim, xdim))
    x_values = np.arange(xdim)
    y_values = np.arange(ydim)

    gg = geogrid.GeoGrid(data_values, x_values, y_values, 3411)
    gg = geogrid.GeoGrid(data_values, x_values, y_values, epsg_3411_as_wkt)
    gg = geogrid.GeoGrid(data_values, x_values, y_values, epsg_3411_as_proj4)


'''
def test_can_declare_GeoGrid_by_fname_tidx():
    """Test can initialize a GeoGrid from scratch"""
    nc_varstring = 'NETCDF:"../../utils/nh_0630.nc":TB'
    nc_tindex = 0
    gg = geogrid.GeoGrid_by_ncref(nc_varstring, nc_tindex)
'''
