"""
geogrid.py

Provicdes GeoGrid class
"""
from osgeo import osr
import gdal
import gdalconst
import numpy as np


class GeoGrid(object):
    """geocode's located-array class"""
    pass

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


def GeoGrid_by_ncref(nc_varstring, nc_timeindex):
    """Return a GeoGrid by reading a netCDF file"""
    try:
        ds = gdal.Open(nc_varstring, gdalconst.GA_ReadOnly)
        assert ds is not None
    except AssertionError:
        raise ValueError('could not open netCDF file: {}'.format(nc_varstring))

    data_array = ds.GetRasterBand(1).ReadAsArray()
    array_shape = data_array.shape

    if len(array_shape) == 2:
        try:
            assert nc_timeindex == 0
        except AssertionError:
            raise ValueError('array is 2D, but timeindex is not 0')

    (ydim, xdim) = array_shape

    print('array shape: {}'.format(array_shape))
    print('len array shape: {}'.format(len(array_shape)))

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
        ds_geotransform[3] + (xdim - 0.5) * ds_geotransform[5],
        ydim)

    return GeoGrid(data_array, x, y, ds.GetProjection())
