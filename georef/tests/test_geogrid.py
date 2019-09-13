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
import osr
import gdal


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


def simple_GeoGrid_values():
    """Returns a simple set of values that can initialize a GeoGrid"""
    xdim = 5
    ydim = 6

    x0 = -163
    dx = 1

    y0 = 64
    dy = -0.5

    geotransform = (x0, dx, 0., y0, 0., dy)

    data_values = np.arange(xdim * ydim).reshape((ydim, xdim))
    x_values = np.linspace(
        x0 + 0.5 * dx,
        x0 + (xdim - 0.5) * dx,
        xdim)
    y_values = np.linspace(
        y0 + 0.5 * dy,
        y0 + (ydim - 0.5) * dy,
        ydim)

    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    wkt = srs.ExportToWkt()

    # Note: GDAL doesn't have an Int32 type...
    if data_values.dtype == np.int64:
        data_values = data_values.astype(np.int32)

    """
    # use local create_geotiff() function
    test_fn = 'intest.tif'
    create_geotiff(test_fn,
        x_values, y_values, data_values, wkt, geotransform)
    raise SystemExit('Stopping after creating geotiff')
    """


    return data_values, x_values, y_values, wkt


def create_geotiff(out_fn, x_values, y_values, data_values, wkt, geotransform):
    """Quick create_geotiff"""
    ds = gdal.GetDriverByName('GTiff').Create(
        out_fn, len(x_values), len(y_values), 1, 1)
    ds.SetGeoTransform(geotransform)
    ds.SetProjection(wkt)
    ds.GetRasterBand(1).WriteArray(data_values)
    ds.FlushCache()


def WKT_are_Proj4_equivalent(wkt1, wkt2):
    """Test that these WKT strings have equivalent proj4 strings"""
    srs1 = osr.SpatialReference()
    srs1.ImportFromWkt(wkt1)
    proj1 = srs1.ExportToProj4()

    srs2 = osr.SpatialReference()
    srs2.ImportFromWkt(wkt2)
    proj2 = srs1.ExportToProj4()

    return proj1 == proj2


def test_can_declare_GeoGrid():
    """Test can initialize a GeoGrid from scratch"""
    xdim = 3
    ydim = 4
    data_values = np.arange(12.).reshape((ydim, xdim))
    x_values = np.arange(xdim)
    y_values = np.arange(ydim)

    test_srs = osr.SpatialReference()
    test_srs.ImportFromEPSG(3411)

    gg = georef.GeoGrid(data_values, x_values, y_values, 3411)
    assert test_srs.ExportToWkt() == gg.wkt

    gg = georef.GeoGrid(data_values, x_values, y_values, epsg_3411_as_wkt)
    assert test_srs.ExportToWkt() == gg.wkt

    gg = georef.GeoGrid(data_values, x_values, y_values, epsg_3411_as_proj4)
    assert WKT_are_Proj4_equivalent(test_srs.ExportToWkt(), gg.wkt)


def test_can_declare_GeoGrid_by_fname_tidx():
    """Test can initialize a GeoGrid from scratch"""
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

    specified_date = datetime.datetime(2018, 1, 8, 12, 0)
    nc_specific_index = georef.geogrid.get_specific_nc_timeindex(
        nc_fname, specified_date)
    assert nc_specific_index == 1


def test_bad_projdef_yields_Warningo():
    """Test that a bad projdef raises a Warning"""
    try:
        dvals, xvals, yvals, wkt = simple_GeoGrid_values()
        gg = georef.GeoGrid(dvals, xvals, yvals, 'non-proj4_string',
                           warn=True)
    except Warning:
        return

    raise ValueError('Warning not raised for bad proj description')


def test_geogrid_to_VRT():
    """Test that we can create a valid in-memory GDAL model from geogrid"""
    dvals, xvals, yvals, wkt = simple_GeoGrid_values()
    print('wkt: {}'.format(wkt))
    gg = georef.GeoGrid(dvals.astype(np.float32), xvals, yvals, wkt, warn=True)
    assert gg.isComplete()

    mem_gg = georef.geogrid.geogrid_as_gdalInMem(gg)
    gg_from_mem = georef.geogrid.geogrid_from_gdalInMem(mem_gg)

    assert gg_from_mem.isComplete()


def test_reproject_geogrid():
    """Test that we can reproject a geogrid"""
    dvals, xvals, yvals, wkt = simple_GeoGrid_values()
    gg = georef.GeoGrid(dvals.astype(np.float32), xvals, yvals, wkt)

    print('in_geotransform:\n  {}'.format(gg.geotransform))
    gg.saveAsGeotiff('test_gg_preproj.tif')

    # Reprojection parameters
    out_xdim = 100
    out_ydim = 100

    out_epsg = 4326
    out_geotransform = (-163., 0.05, 0.0, 65.0, 0.0, -0.05)

    #out_epsg = 3411
    #out_geotransform = (-2737500.0, 2000.0, 0.0, 1412500.0, 0.0, -2000.0)

    gg_reproj = georef.geogrid.reproject_GeoGrid(
        gg, out_epsg, out_xdim, out_ydim, out_geotransform)

    assert gg_reproj.isComplete()
    gg_reproj.saveAsGeotiff('test_reproj.tif')

    """
    gg2 = gg_reproj
    #print('gg2.x:\n  {}'.format(gg2.x))
    #print('gg2.y:\n  {}'.format(gg2.y))
    print('gg2.data:\n  {}'.format(gg2.data_array))
    print('gg2.min:\n  {}'.format(gg2.data_array.min()))
    print('gg2.max:\n  {}'.format(gg2.data_array.max()))
    print('gg2.geotransform:  {}'.format(gg2.geotransform))
    print('gg2.xdim:  {}'.format(gg2.xdim))
    print('gg2.ydim:  {}'.format(gg2.ydim))
    print('right edge x: {}'.format(
        gg2.x[0] + (gg2.xdim - 0.5) * gg2.geotransform[1]))
    print('bottom edge y: {}'.format(
        gg2.y[0] + (gg2.ydim - 0.5) * gg2.geotransform[5]))
    gg_reproj.saveAsGeotiff('test_gg_reproj.tif')
    """

    """
    print('out_xdim: {}'.format(out_xdim))
    print('out_ydim: {}'.format(out_ydim))
    print('gg_reproj xdim: {}'.format(len(gg_reproj.x)))
    print('gg_reproj ydim: {}'.format(len(gg_reproj.y)))

    print('reprojed array:\n  {}'.format(gg_reproj.data_array))
    """


'''
def test_saveasGeotiff():
    """Test that we can save a GeoGrid as a geotiff"""
    overwrite=True
    test_fn = 'test_geotiff.tif'
    if overwrite or not os.path.isfile(test_fn):
        dvals, xvals, yvals, wkt = simple_GeoGrid_values()
        gg = georef.GeoGrid(dvals, xvals, yvals, wkt)
        gg.saveAsGeotiff(test_fn)
    else:
        print('Skipping saveAsGeotiff test because testfile already exists')


def test_convert_grid_corners():
    """Test can convert grid corners to GeoPoints"""
    dvals, xvals, yvals, wkt = simple_GeoGrid_values()
    gg = georef.GeoGrid(dvals, xvals, yvals, wkt)
    i0=1
    j0=1
    print('x[{}]: {}'.format(i0, gg.x[i0]))
    print('y[{}]: {}'.format(j0, gg.y[j0]))
    print('d[{}, {}]: {}'.format(i0, j0, gg.data_array[j0, i0]))
    gp = georef.GeoPoint(gg.srs, gg.x[i0], gg.y[j0], gg.data_array[j0, i0])

'''
