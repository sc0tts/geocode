"""
geogrid.py

Provicdes GeoGrid class
"""
from osgeo import osr
import netCDF4


class GeoGrid(object):
    """geocode's located-array class"""

    def __init__(self, srs_string=None,
                 array=None, x=None, y=None):
        """GeoGrid is an array with projection"""
        self._srs = osr.SpatialReference()

        #if fname is not None:

        if srs_string is not None:
            try:
                if self._srs.ExportToWkt() == "" and \
                   isinstance(srs_string, int):
                    self._srs.ImportFromEPSG(srs_string)
            except TypeError:
                print('Ignoring TypeError when importing from epsg')

            try:
                if self._srs.ExportToWkt() == "" and \
                   isinstance(srs_string, str) and \
                   'PROJ' in srs_string:
                    self._srs.ImportFromWkt(srs_string)
            except TypeError:
                print('Ignoring TypeError when importing from WKT')

            try:
                if self._srs.ExportToWkt() == "" and \
                   isinstance(srs_string, str) and \
                   '+proj' in srs_string:
                    self._srs.ImportFromProj4(srs_string)
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


'''
    def isValid(self):
        """Check if all variables have a value"""
        return self._srs is not None and \
               self._value is not None and \
               self._x is not None and \
               self._y is not None


def transform_point(src_gp, dst_srs):
    """Transform one gp into another gp with specified SRS"""
    transform = osr.CoordinateTransformation(src_gp.srs, dst_srs)
    # Note: ignoring dst_z
    dst_x, dst_y, dst_z = transform.TransformPoint(src_gp.x, src_gp.y)

    return GeoPoint(dst_srs, dst_x, dst_y, src_gp.value)
'''
