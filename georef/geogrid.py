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
import datetime


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

    # Field name is specifiec as the nc_fname
    # Each slice of the 2D array is a single Raster here
    # ...
    dataset_array = ds.ReadAsArray()
    if len(dataset_array.shape) == 3:
        assert nc_timeindex < ds.RasterCount
        dataset_array = dataset_array[nc_timeindex]

    assert len(dataset_array.shape) == 2

    ydim, xdim = dataset_array.shape

    # Sample geotransform values:
    #     (-3850000.0, 25000.0, 0.0, 5850000.0, 0.0, -25000.0)
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

    return GeoGrid(dataset_array, x, y, ds.GetProjection())


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


def get_nc_timevalues(fname, time_varname='time'):
    """Returns the time values as datetime objects"""
    nc_fid = netCDF4.Dataset(fname, 'r')
    try:
        assert nc_fid.variables[time_varname] is not None
    except AssertionError:
        print('nc file does not have variable/dimension: {}'.format(
                  time_varname))
        raise ValueError('time variable not found')

    nc_time_values = nc_fid.variables[time_varname]
    time_values = netCDF4.num2date(nc_time_values[:],
                                   nc_time_values.units,
                                   nc_time_values.calendar)

    return time_values


def time_difference(tval1, tval2):
    """Compute the difference between two time values
    If either input is a datetime.date, then the difference is in days,
    otherwise the difference is in seconds
    """
    assert isinstance(tval1, datetime.date)
    assert isinstance(tval2, datetime.date)

    if not isinstance(tval1, datetime.datetime) or \
       not isinstance(tval2, datetime.datetime):
        # Comparing dates
        compare_value_1 = datetime.date(tval1.year, tval1.month, tval1.day)
        compare_value_2 = datetime.date(tval2.year, tval2.month, tval2.day)

        days_diff = abs((compare_value_1 - compare_value_2).days)
        return days_diff


    # Otherwise, we are comparing datetimes...
    if not isinstance(tval1, datetime.datetime):
        compare_value_1 = \
            datetime.datetime(tval1.year, tval1.month, tval1.day, 12, 0)
    else:
        compare_value_1 = tval1

    if not isinstance(tval2, datetime.datetime):
        compare_value_2 = \
            datetime.datetime(tval2.year, tval2.month, tval2.day, 12, 0)
    else:
        compare_value_2 = tval2

    return abs((compare_value_1 - compare_value_2).total_seconds())


def get_specific_nc_timeindex(fname,
                              time_value,
                              time_varname='time'):
    """Return the index of the nearest value
    if time_value is a datetime.date, find first nearest-day
    if time_value is a datetime.datetime, find first nearest-seconds
    """

    assert isinstance(time_value, datetime.date)

    nc_fid = netCDF4.Dataset(fname, 'r')
    time_values = netCDF4.num2date(nc_fid.variables[time_varname][:],
                                   nc_fid.variables[time_varname].units,
                                   nc_fid.variables[time_varname].calendar)
    nearest_index = 0
    nearest_value = time_values[nearest_index]
    nearest_diff = time_difference(nearest_value, time_value)
    for idx, tvalue in enumerate(time_values):
        this_diff = time_difference(tvalue, time_value)
        if this_diff < nearest_diff:
            nearest_index = idx
            nearest_value = time_values[idx]
            nearest_diff = this_diff

    return nearest_index


"""
Reference section

Typeical Python object dir() of a dimension in a georeferenced nc dataset:


    time_values = nc_fid.variables[time_varname]
    dir(time_values)

    ['__array__', '__class__', '__delattr__', '__delitem__', '__dir__',
    '__doc__', '__eq__', '__format__', '__ge__', '__getattr__',
    '__getattribute__', '__getitem__', '__gt__', '__hash__', '__init__',
    '__init_subclass__', '__le__', '__len__', '__lt__', '__ne__', '__new__',
    '__orthogonal_indexing__', '__reduce__', '__reduce_ex__', '__repr__',
    '__setattr__', '__setitem__', '__sizeof__', '__str__', '__subclasshook__',
    '__unicode__', '_assign_vlen', '_check_safecast', '_cmptype', '_enumtype',
    '_get', '_getdims', '_getname', '_grp', '_grpid', '_has_lsd',
    '_iscompound', '_isenum', '_isprimitive', '_isvlen', '_name',
    '_ncstring_attrs__', '_nunlimdim', '_put', '_toma', '_use_get_vars',
    '_varid', '_vltype', 'always_mask', 'assignValue', 'axis', 'calendar',
    'chartostring', 'chunking', 'datatype', 'delncattr', 'dimensions', 'dtype',
    'endian', 'filters', 'getValue', 'get_dims', 'get_var_chunk_cache',
    'getncattr', 'group', 'mask', 'name', 'ncattrs', 'ndim', 'renameAttribute',
    'scale', 'set_always_mask', 'set_auto_chartostring', 'set_auto_mask',
    'set_auto_maskandscale', 'set_auto_scale', 'set_collective',
    'set_ncstring_attrs', 'set_var_chunk_cache', 'setncattr',
    'setncattr_string', 'setncatts', 'shape', 'size', 'standard_name', 'units',
    'use_nc_get_vars']

"""
