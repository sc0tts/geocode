"""
geopoint.py

Provicdes GeoPoint class
"""
from osgeo import osr


class GeoPoint(object):
    """geocode's single location value"""
    def __init__(self, srs=None, x=None, y=None, value=None):
        """GeoPoint has value at (x, y) in projection srs"""
        self._srs = srs
        self._value = value
        self._x = x  # Floating point value of k
        self._y = y


    @property
    def x(self):
        """Return the x coordinate"""
        return self._x


    @property
    def y(self):
        """Return the y coordinate"""
        return self._y


    @property
    def value(self):
        """Return the value at the GeoPoint"""
        return self._value


    @property
    def srs(self):
        """Return the srs"""
        return self._srs


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
