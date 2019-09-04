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
    gp_empty = geogrid.GeoPoint()


def test_can_declare_geopoint():
    """Test can declare a proper GeoPoint class"""
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(3411)
    xval=-3800000
    yval=5800000
    init_value=100
    gp = geogrid.GeoPoint(
        srs=srs,
        x=-3800000,
        y=5800000,
        value=100
    )
    assert gp.x == xval
    assert gp.y == yval
    assert gp.value == init_value
    assert gp.srs == srs


def test_geopoint_isValid():
    gp_empty = geogrid.GeoPoint()
    assert gp_empty.isValid() == False

    gp_minimal = geogrid.GeoPoint(
        srs=osr.SpatialReference(),
        x=0,
        y=0,
        value=0)
    assert gp_minimal.isValid()
