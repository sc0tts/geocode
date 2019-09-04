"""
test_initial.py

Initial tests to ensure project layout is adequate
"""

#import numpy as np
from osgeo import osr
#from osgeo import gdal
import geogrid

def test_imports():
    """Test that standard imports raise no exceptions"""
    pass


def test_geopoint_exists():
    """Test can declare an empty GeoPoint class"""
    gp_empty = geogrid.geopoint.GeoPoint()


def test_geopoint_can_be_declared():
    """Test can declare a proper GeoPoint class"""
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(3411)
    xval=-3800000
    yval=5800000
    init_value=100
    gp = geogrid.geopoint.GeoPoint(
        srs=srs,
        x=-3800000,
        y=5800000,
        value=100)
    assert gp.x == xval
    assert gp.y == yval
    assert gp.value == init_value
    assert gp.srs == srs
    print(gp.srs)


