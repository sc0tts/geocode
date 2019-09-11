"""
test_geogrid.py

Tests for geogrid's GeoGrid class

Note: some of these tests require a local netCDF file
"""

import os
import georef
import numpy as np
import netCDF4
import datetime


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

default_georef_testnc_fname = \
    '/home/scotts/bdi_v1/georef/test_ncfile.nc'
default_georef_testnc_fname_2 = \
    '/home/scotts/bdi_v1/georef/utils/nh_0630.nc'
default_georef_testnc_fname_3 = \
    '/home/scotts/bdi_v1/georef/utils/nh_0051.nc'

default_georef_testnc_varname = 'tb'
default_georef_testnc_varname_2 = 'TB'
default_georef_testnc_varname_3 = 'sea_ice_concentration'


def test_can_declare_GeoGrid():
    """Test can initialize a GeoGrid from scratch"""
    xdim = 3
    ydim = 4
    data_values = np.arange(12.).reshape((ydim, xdim))
    x_values = np.arange(xdim)
    y_values = np.arange(ydim)

    gg = georef.GeoGrid(data_values, x_values, y_values, 3411)
    gg = georef.GeoGrid(data_values, x_values, y_values, epsg_3411_as_wkt)
    gg = georef.GeoGrid(data_values, x_values, y_values, epsg_3411_as_proj4)


def test_can_declare_GeoGrid_by_fname_tidx():
    """Test can initialize a GeoGrid from scratch"""
    #nc_varstring = 'dumb'
    #nc_varstring = 'NETCDF:"./utils/nh_0630.nc":TB'
    nc_varstring = 'NETCDF:\"{}\":{}'.format(default_georef_testnc_fname,
                                             default_georef_testnc_varname)
    nc_tindex = 0
    gg = georef.geogrid.GeoGrid_by_ncref(nc_varstring, nc_tindex)


def test_read_from_netCDF_file():
    """Test that we can read from a specified netCDF file"""
    raw_nc_fname = default_georef_testnc_fname
    quoted_nc_fname = '\"{}\"'.format(raw_nc_fname)
    if not os.path.isfile(raw_nc_fname):
        print('test nc file does not exist, skipping: {}'.format(
            raw_nc_fname))
    nc_varname = default_georef_testnc_varname
    nc_timeindex = 9
    nc_filestring = 'NETCDF:{}:{}'.format(quoted_nc_fname, nc_varname)
    nc_geogrid = georef.geogrid.GeoGrid_by_ncref(nc_filestring, nc_timeindex)
    #print('type(nc_geogrid): {}'.format(type(nc_geogrid)))
    assert str(type(nc_geogrid)) == "<class 'georef.geogrid.GeoGrid'>"


def test_read_from_netCDF_file_alt():
    """Test that we can read from a specified netCDF file"""
    raw_nc_fname = default_georef_testnc_fname_2
    quoted_nc_fname = '\"{}\"'.format(raw_nc_fname)
    if not os.path.isfile(raw_nc_fname):
        print('test nc file does not exist, skipping: {}'.format(
            raw_nc_fname))
    nc_varname = default_georef_testnc_varname_2
    nc_timeindex = 0
    nc_filestring = 'NETCDF:{}:{}'.format(quoted_nc_fname, nc_varname)
    nc_geogrid = georef.geogrid.GeoGrid_by_ncref(nc_filestring, nc_timeindex)
    #print('type(nc_geogrid): {}'.format(type(nc_geogrid)))
    assert str(type(nc_geogrid)) == "<class 'georef.geogrid.GeoGrid'>"


def test_get_ncvarname_from_stdname():
    """Test determination of nc variable name given fname and standard name"""
    # Note: hard-coding to 0051-equivalent file
    raw_nc_fname = '/home/scotts/bdi_v1/georef/utils/nh_0051.nc'
    if not os.path.isfile(raw_nc_fname):
        print('test nc file does not exist, skipping: {}'.format(
            raw_nc_fname))
        return
    stdname = 'sea_ice_area_fraction'
    expected_varname = 'conc'
    nc_varname = georef.geogrid.get_nc_var_stdname(raw_nc_fname, stdname)
    assert nc_varname == expected_varname


def test_can_get_time_from_nc_file():
    """Test retrieval of time values from netCDF file"""
    # Homemade .nc file
    nc_fname = default_georef_testnc_fname
    if not os.path.isfile(nc_fname):
        print('Testing nc file does not exist, skipping: {}'.format(nc_fname))
        return

    nc_time_values = georef.geogrid.get_nc_timevalues(nc_fname)

    assert (10,) == nc_time_values.shape

    if len(nc_time_values) > 1:
        assert type(nc_time_values[0]) == \
               type(datetime.datetime(2019, 1, 1, 12, 0))

    # NSIDC 0630 .nc file
    nc_fname = default_georef_testnc_fname_2
    if not os.path.isfile(nc_fname):
        print('Testing nc file does not exist, skipping: {}'.format(nc_fname))
        return

    nc_time_values = georef.geogrid.get_nc_timevalues(nc_fname)

    assert (1,) == nc_time_values.shape

    # Verify that we get a datetime.datetime object
    assert type(nc_time_values[0]) == \
           type(datetime.datetime(2019, 1, 1, 12, 0))

    # NSIDC 0051 .nc file
    nc_fname = default_georef_testnc_fname_3
    if not os.path.isfile(nc_fname):
        print('Testing nc file does not exist, skipping: {}'.format(nc_fname))
        return

    nc_time_values = georef.geogrid.get_nc_timevalues(nc_fname)

    assert (1,) == nc_time_values.shape

    # Verify that we get a datetime.datetime object
    assert type(nc_time_values[0]) == \
           type(datetime.datetime(2019, 1, 1, 12, 0))


def test_find_specific_nc_time_index():
    """Test that we can return the time value with a specified index"""
    # Homemade .nc file
    nc_fname = default_georef_testnc_fname
    if not os.path.isfile(nc_fname):
        print('Testing nc file does not exist, skipping: {}'.format(nc_fname))
        return

    """
    specified_date = datetime.date(2018, 1, 8)
    nc_specific_index = georef.geogrid.get_specific_nc_timeindex(
        nc_fname, specified_date)
    assert nc_specific_index == 1

    """
    specified_date = datetime.datetime(2018, 1, 8, 12, 0)
    nc_specific_index = georef.geogrid.get_specific_nc_timeindex(
        nc_fname, specified_date)
    assert nc_specific_index == 1
