"""
test_initial.py

Initial tests to ensure project layout is adequate
"""

from osgeo import osr
import georef
from numpy.testing import assert_almost_equal
#from osgeo import gdal
#import numpy as np


# Helper functions for tests
polar_stereo_north_srs = osr.SpatialReference()
polar_stereo_north_srs.ImportFromEPSG(3411)

latlon_srs = osr.SpatialReference()
latlon_srs.ImportFromEPSG(4326)


def test_imports():
    """Test that standard imports raise no exceptions"""
    pass


def test_geopoint_exists():
    """Test can declare an empty GeoPoint class"""
    gp_empty = georef.GeoPoint()


def test_can_declare_geopoint():
    """Test can declare a proper GeoPoint class"""
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(3411)
    xval=-3800000
    yval=5800000
    init_value=100
    gp = georef.GeoPoint(
        srs=srs,
        x=xval,
        y=yval,
        value=init_value
    )
    assert gp.x == xval
    assert gp.y == yval
    assert gp.value == init_value
    assert gp.srs == srs


def create_geopoint_on_psn(
        srs=polar_stereo_north_srs,
        xval=-38000,
        yval=5800000,
        init_value=100):
    """Create a sample geopoint on North Polar Stereo projection"""
    return georef.GeoPoint(srs=srs, x=xval, y=yval, value=init_value)


def test_geopoint_isValid():
    gp_empty = georef.GeoPoint()
    assert gp_empty.isValid() == False

    gp_minimal = create_geopoint_on_psn()
    assert gp_minimal.isValid()


def test_can_transform_geopoint():
    """Test conversion of geopoint from one SRS to another SRS"""
    gp_on_psn = create_geopoint_on_psn()
    gp_transformed = georef.geopoint.transform_point(gp_on_psn, latlon_srs)
    assert gp_transformed.srs == latlon_srs
    assert gp_transformed.value == gp_on_psn.value
    assert_almost_equal(gp_transformed.x, 39.80610930, decimal=5)
    assert_almost_equal(gp_transformed.y, 135.3753808, decimal=5)


