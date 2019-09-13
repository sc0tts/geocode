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


gdal.UseExceptions()

def assign_projection_to_srs(srs, projdef, verbose=False, warn=True):
    """Define the srs by interpreting 'projdef' as epsg code, WKT, or PROJ4"""
    try:
        if srs.ExportToWkt() == "" and \
           isinstance(projdef, int):
            srs.ImportFromEPSG(projdef)
    except TypeError:
        if verbose:
            print('Ignoring TypeError when importing from epsg')

    try:
        if srs.ExportToWkt() == "" and \
           isinstance(projdef, str) and \
           'PROJ' in projdef:
            srs.ImportFromWkt(projdef)
    except TypeError:
        if verbose:
            print('Ignoring TypeError when importing from WKT')

    try:
        if srs.ExportToWkt() == "" and \
           isinstance(projdef, str) and \
           '+proj' in projdef:
            srs.ImportFromProj4(projdef)
    except TypeError:
        if verbose:
            print('Ignoring TypeError when importing from Proj4')

    if warn:
        if srs.ExportToWkt() == "":
            if verbose:
                print('projdef was: {}'.format(projdef))
            raise Warning('SRS not determined in GeoGrid initialization')


def get_gdal_datatype(in_datatype):
    """Returns the GDAL data type for this data type"""
    if in_datatype == 'float64':
        return gdalconst.GDT_Float64
    elif in_datatype == 'float32':
        return gdalconst.GDT_Float32
    elif in_datatype == 'int32':
        return gdalconst.GDT_Int32
    else:
        raise ValueError(
            'Unrecognized data type in get_gdal_datatype():\n  {}'.format(
                in_datatype))


class GeoGrid(object):
    """geocode's located-array class"""
    def __init__(self, array=None, x=None, y=None, proj_definition=None,
                verbose=False, warn=False, nodataval=None):
        """GeoGrid is an array with projection"""
        self._srs = osr.SpatialReference()

        if proj_definition is not None:
            assign_projection_to_srs(self._srs, proj_definition,
                                     verbose=verbose, warn=warn)

        if warn:
            if self._srs.ExportToWkt() == "":
                raise Warning('SRS not determined in GeoGrid initialization')

        if nodataval is None and isinstance(array[0, 0], float):
            self._NoDataValue = np.nan
        elif nodataval is not None:
            self._NoDataValue = nodataval
        else:
            self._NoDataValue = None

        self._x = x
        self._y = y
        self._set_geotransform(x, y)

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

    def _set_geotransform(self, xvals, yvals):
        """Compute the geotransform string from the x and y values"""
        xdim = len(xvals)
        ydim = len(yvals)

        x0 = xvals[0]
        y0 = yvals[0]

        dx = xvals[1] - xvals[0]
        dy = yvals[1] - yvals[0]

        x_leftedge = x0 - dx / 2
        y_topedge = y0 + dx / 2

        xlast = x0 + (xdim -1) * dx
        ylast = y0 + (ydim -1) * dy

        assert abs(xlast - xvals[xdim - 1]) < \
               abs(max(xlast, xvals[xdim - 1])) / 1000.

        self._geotransform = (x_leftedge, dx, 0., y_topedge, 0., dy)


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

    @property
    def xdim(self):
        """Return the number of x values"""
        return len(self._x)

    @property
    def ydim(self):
        """Return the number of y values"""
        return len(self._y)

    @property
    def geotransform(self):
        """Return the 'geotransform' equivalent of the coordinates"""
        return self._geotransform

    def saveAsGeotiff(self, geotiff_fname, overwrite=True, verbose=False):
        """Save the GeoGrid as a geotiff with the given fname"""
        if not overwrite and os.path.exists(geotiff_fname):
            raise ValueError('output geotiff already exists:\n  {}'.format(
                geotiff_fname))

        # Note: Geotiffs can't be 64-bit floats
        src_datatype = get_gdal_datatype(self._data_array.dtype)
        if src_datatype == gdalconst.GDT_Float64:
            if verbose:
                print('  converting data from Float64 to Float32')
            src_datatype = gdalconst.GDT_Float32
            local_data_array = self._data_array.astype(np.float32)
        else:
            local_data_array = self._data_array

        # Encode the source as a Geotiff
        src = gdal.GetDriverByName('Gtiff').Create(
            geotiff_fname,
            len(self._x),
            len(self._y),
            1,
            src_datatype)

        src.SetGeoTransform(self._geotransform)
        src.SetProjection(self.srs.ExportToWkt())

        src_raster = src.GetRasterBand(1)
        if self._NoDataValue is not None:
            src_raster.SetNoDataValue(self._NoDataValue)

        src_raster.WriteArray(local_data_array)
        src_raster.FlushCache()

        src_raster = None


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


def warp_GeoGrid(geogrid_in, srs_string,
                     out_xdim=None, out_ydim=None, out_geotransform=None,
                     out_nodata_value=None, interp_method=None):
    """warp geogrid_in onto srs_out"""
    in_srs = geogrid_in.srs
    out_srs = osr.SpatialReference()
    assign_projection_to_srs(out_srs, srs_string)

    if interp_method is None:
        """
        I think the options are:
            GRA_Average
            GRA_Bilinear
            GRA_Cubic
            GRA_CubicSpline
            GRA_NearestNeighbour
        """
        interp_method = gdalconst.GRA_NearestNeighbour
    elif isinstance(interp_method, str) and \
            interp_method.lower() in ('nearest', 'nearestneighbor',):
        interp_method = gdalconst.GRA_NearestNeighbour
    elif isinstance(interp_method, str) and \
            interp_method.lower() in ('linear', 'bilinear',):
        interp_method = gdalconst.GRA_Bilinear
    elif isinstance(interp_method, str) and \
            interp_method.lower() in ('spline', 'cubicspline',):
        interp_method = gdalconst.GRA_CubicSpline
    elif isinstance(interp_method, str) and \
            interp_method.lower() in ('avg', 'average',):
        interp_method = gdalconst.GRA_Average
    else:
        raise ValueError('interp_method not recognized: {}'.format(
            interp_method))

    # Autodetection of transform coords is not supported yet:
    if out_xdim is None or \
       out_ydim is None or \
       out_geotransform is None:
        raise ValueError(
            'Currently need to specify out_xdim, out_ydim, out_geotransform')

    if out_nodata_value is None:
        try:
            out_nodata_value = geogrid_in.NoDataValue()
        except AttributeError:
            out_nodata_value = np.nan

    src_array_dtype = geogrid_in.data_array.dtype

    if src_array_dtype == 'float64':
        src_datatype = gdalconst.GDT_Float64
    elif src_array_dtype == 'int32':
        src_datatype = gdalconst.GDT_Int32
    else:
        raise ValueError('Could not determine src_datatype:\n{}'.format(
            src_array_dtype))


    src_srs = geogrid_in.srs
    #dst_srs = out_srs
    dst_srs = src_srs

    src_geotransform = geogrid_in.geotransform
    dst_geotransform = geogrid_in.geotransform

    # I think life will just be easier if we only work in Float32....
    dst_datatype = gdalconst.GDT_Float32

    src = gdal.GetDriverByName('MEM').Create(
        '', geogrid_in.xdim, geogrid_in.ydim, 1, src_datatype)

    dst = gdal.GetDriverByName('MEM').Create(
        '', out_xdim, out_ydim, 1, dst_datatype)

    print('dst init RasterXSize: {}'.format(dst.RasterXSize))
    print('dst init RasterYSize: {}'.format(dst.RasterYSize))

    src.SetGeoTransform(src_geotransform)
    src.SetProjection(geogrid_in.srs.ExportToWkt())

    dst.SetGeoTransform(out_geotransform)
    dst.SetProjection(out_srs.ExportToWkt())
    #dst.SetGeoTransform(src_geotransform)
    #dst.SetProjection(geogrid_in.srs.ExportToWkt())

    src_raster = src.GetRasterBand(1)
    src_raster.SetNoDataValue(out_nodata_value)
    src_raster.WriteArray(geogrid_in.data_array)

    print('src_raster:')
    print(src_raster.ReadAsArray()[:])

    dst_raster = dst.GetRasterBand(1)
    dst_raster.SetNoDataValue(out_nodata_value)

    print('src geotransform: {}'.format(src.GetGeoTransform()))
    print('dst geotransform: {}'.format(dst.GetGeoTransform()))

    print('out_srs: {}'.format(out_srs.ExportToWkt()))
    print('about to warp')
    gdal.Warp(dst, src)

    #gdal.warpImage(src,
    #                    dst,
    #                    geogrid_in.srs.ExportToWkt(),
    #                    out_srs.ExportToWkt(),
    #                    interp_method)

    print('done with warp')
    print('dst_raster:')
    print(dst_raster.ReadAsArray()[:])

    # Now, create a GeoGrid object from the remainder...
    print('out_xdim: {}'.format(out_xdim))
    print('out_ydim: {}'.format(out_ydim))
    print('out_geot: {}'.format(out_geotransform))
    print('RasterXSize: {}'.format(dst.RasterXSize))
    print('RasterYSize: {}'.format(dst.RasterYSize))

    xvals = np.linspace(
        out_geotransform[0] + 0.5 * out_geotransform[1],
        out_geotransform[0] + (out_xdim - 0.5) * out_geotransform[1],
        out_xdim)
    yvals = np.linspace(
        out_geotransform[3] + 0.5 * out_geotransform[5],
        out_geotransform[3] + (out_ydim - 0.5) * out_geotransform[5],
        out_ydim)
    data_array = dst_raster.ReadAsArray()
    print('dst array shape: {}'.format(data_array.shape))
    print('dst geotransform: {}'.format(dst.GetGeoTransform()))
    print('src array min: {}'.format(geogrid_in.data_array.min()))
    print('src array max: {}'.format(geogrid_in.data_array.max()))
    print('dst array min: {}'.format(data_array.min()))
    print('dst array max: {}'.format(data_array.max()))

    print(xvals)

    return_gg = GeoGrid(data_array, xvals, yvals, out_srs.ExportToWkt())
    return return_gg


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
