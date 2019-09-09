"""
geogrid.py

Provicdes GeoGrid class
"""
from osgeo import osr
import gdal
import gdalconst
import numpy as np
from pprint import pprint
import netCDF4


class GeoGrid(object):
    """geocode's located-array class"""
    def __init__(self, array=None, x=None, y=None, proj_definition=None):
        """GeoGrid is an array with projection"""
        self._srs = osr.SpatialReference()

        if proj_definition is not None:
            try:
                if self._srs.ExportToWkt() == "" and \
                   isinstance(proj_definition, int):
                    self._srs.ImportFromEPSG(proj_definition)
            except TypeError:
                print('Ignoring TypeError when importing from epsg')

            try:
                if self._srs.ExportToWkt() == "" and \
                   isinstance(proj_definition, str) and \
                   'PROJ' in proj_definition:
                    self._srs.ImportFromWkt(proj_definition)
            except TypeError:
                print('Ignoring TypeError when importing from WKT')

            try:
                if self._srs.ExportToWkt() == "" and \
                   isinstance(proj_definition, str) and \
                   '+proj' in proj_definition:
                    self._srs.ImportFromProj4(proj_definition)
            except TypeError:
                print('Ignoring TypeError when importing from Proj4')

        if self._srs.ExportToWkt() == "":
            #print('Warning: SRS not determined in GeoGrid initialization')
            raise Warning('SRS not determined in GeoGrid initialization')

        self._x = x
        self._y = y

        self._data_array = array
        try:
            assert self._data_array.shape[1] == len(self._x)
        except AssertionError:
            raise ValueError(
                'length of x ({}) does not correspond to array dim ({})'
                .format(self._data_array.shape[1], len(self._x)))
        try:
            assert self._data_array.shape[0] == len(self._y)
        except AssertionError:
            raise ValueError(
                'length of y ({}) does not correspond to array dim ({})'
                .format(self._data_array.shape[0], len(self._y)))


    @property
    def data_array(self):
        """Return the array of values"""
        return self._data_array


    @property
    def srs(self):
        """Return the srs"""
        return self._srs


    @property
    def x(self):
        """Return the x values"""
        return self._x


    @property
    def y(self):
        """Return the y values"""
        return self._y


def remove_quotes(fname_string):
    """Remove quote marks from beginning and end of string"""
    return fname_string.replace('"', '').replace("'", '')


def GeoGrid_by_ncref(nc_varstring, nc_timeindex=0):
    """Return a GeoGrid by reading a netCDF file"""
    try:
        ds = gdal.Open(nc_varstring, gdalconst.GA_ReadOnly)
        assert ds is not None
    except AssertionError:
        raise ValueError('could not open netCDF file: {}'.format(nc_varstring))

    nc_string, nc_fname, nc_fieldname = nc_varstring.split(':')
    print('nc_string: {}'.format(nc_string))
    print('nc_fname: {}'.format(nc_fname))
    print('nc_fieldname: {}'.format(nc_fieldname))

    # Field name is specifiec as the nc_fname
    # Each slice of the 2D array is a single Raster here
    # ...
    dataset_array = ds.ReadAsArray()
    print('dataset_array.shape: {}'.format(dataset_array.shape))
    return
    try:
        data_array = ds.ReadAsArray()[nc_timeindex]
    except ValueError:
        raise ValueError('value error!')

    #pprint(dir(gdal))

    nc_fname_raw = remove_quotes(nc_fname)

    print(data_array.shape)
    if len(data_array.shape) == 2:
        try:
            assert nc_timeindex == 0
        except AssertionError:
            #raise ValueError('array is 2D, but timeindex is not 0')
            print('array is 2D, but timeindex is not 0')
    else:
        assert nc_timeindex < ds.RasterCount

    print('data_array shape: {}'.format(data_array.shape))
    print('         (9, 10): {}'.format((9, 10)))
    (ydim, xdim) = (9, 10)
    ydim, xdim = (9, 10)
    das = data_array.shape
    print('das: {}'.format(das))
    (ydim, xdim) = das
    print('nc_varstring: {}'.format(nc_varstring))
    #print('xdim: {}'.format(xdim))
    #print('ydim: {}'.format(ydim))
    return

    # E.g.: (-3850000.0, 25000.0, 0.0, 5850000.0, 0.0, -25000.0)
    ds_geotransform = ds.GetGeoTransform()
    try:
        assert ds_geotransform[2] == 0
        assert ds_geotransform[4] == 0
    except AssertionError:
        raise ValueError('only 1D x,y values are supported')

    x = np.linspace(
        ds_geotransform[0] + 0.5 * ds_geotransform[1],
        ds_geotransform[0] + (xdim - 0.5) * ds_geotransform[1],
        xdim)
    y = np.linspace(
        ds_geotransform[3] + 0.5 * ds_geotransform[5],
        ds_geotransform[3] + (ydim - 0.5) * ds_geotransform[5],
        ydim)

    return GeoGrid(data_array, x, y, ds.GetProjection())


def get_nc_var_stdname(nc_fname, nc_stdvarname):
    """Return the name of ncfile var with given std name"""
    # Note: this may yield the wrong results if more than one var
    #       has the same standard name in a file?
    # Note: Have seen two different decodings of netcdf files
    #  - SUBDATASET...
    #  - <varname>#<varmetadata>

    try:
        ds = netCDF4.Dataset(nc_fname, 'r')
    except ValueError:
        print('ValueError trying to open: {}'.format(nc_fname))

    ds_info = gdal.Info(nc_fname)
    ds_info_lines = ds_info.split('\n')

    nc_varline = None
    nc_alt_var_name = None

    for i, line in enumerate(ds_info_lines):
        if ' {} '.format(nc_stdvarname) in line:
            if nc_varline is None:
                nc_varline = i - 1
            else:
                print('Warning: additional stdname also found')

        if '#standard_name={}'.format(nc_stdvarname) in line:
            nc_alt_var_name = line.split('#')[0]

    if nc_alt_var_name is not None:
        return nc_alt_var_name.lstrip()
    elif nc_varline is not None:
        colon_splits = ds_info_lines[nc_varline].split(':')
        nc_varname = colon_splits[-1]
        return nc_varname
    else:
        # Could not find standard name
        raise ValueError('Could not find stdname {} in {}'.format(
            nc_stdvarname, nc_fname))
