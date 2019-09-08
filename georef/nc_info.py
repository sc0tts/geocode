"""
nc_info.py

Support routines for georef netCDF routines
"""

import os
import netCDF4
import datetime as dt
import numpy as np


def write_simple_netCDF_file(nc_fname, overwrite=False):
    """An example of writing a simply-referenced netCDF file

    Sample panoply settings:
        Stereographic
        -160E  63N
        Edge Angle: 5.0
        Color: 29 to 170
    """

    try:
        assert not os.path.isfile(nc_fname)
    except AssertionError:
        if overwrite:
            os.remove(nc_fname)
        else:
            print('nc_fname exists, skipping creation: {}'.format(nc_fname))

    # Data extracted from 'nt_20180101_f18_nrt_n.bin'
    #  at data[177:185+1, 44:53+1]
    data_values = np.array(\
        [[254, 254, 254, 254, 253, 173, 155, 139, 133, 142],
         [254, 254, 254, 254, 253, 162, 141, 138, 145, 152],
         [254, 254, 254, 254, 254, 253, 141, 143, 151, 150],
         [254, 254, 254, 254, 254, 253, 153, 149, 154, 154],
         [254, 254, 254, 254, 253, 160, 142, 144, 153, 154],
         [254, 254, 254, 254, 253, 158, 131, 129, 138, 144],
         [254, 254, 254, 254, 254, 253, 153, 147, 148, 156],
         [254, 254, 254, 254, 254, 254, 253, 253, 253, 174],
         [254, 254, 254, 254, 254, 254, 254, 254, 254, 253]],
        dtype=np.float32)

    # This example uses a subset of NSIDC's North Polar Stereo grid
    x_left_map_edge = -3850000
    x_grid_step = 25000
    y_top_map_edge = 5850000
    y_grid_step = -25000

    i0_index = 44
    j0_index = 177
    i1_index = 53
    j1_index = 185

    xdim = i1_index - i0_index + 1
    ydim = j1_index - j0_index + 1

    x_left_edge = x_left_map_edge + x_grid_step * i0_index
    y_top_edge = y_top_map_edge + y_grid_step * j0_index

    # Open the netCDF dataset on-disk
    if not overwrite:
        if os.path.isfile(nc_fname):
            print('test file exists, skipping: {}'.format(nc_fname))
            return

    ds = netCDF4.Dataset(nc_fname, 'w')

    # Define dimensions
    x = ds.createDimension('x', xdim)
    y = ds.createDimension('y', ydim)
    time = ds.createDimension('time', None)

    # Assign variables
    xs = ds.createVariable('x', np.float64, ('x',))
    xs_values = np.linspace(x_left_edge + x_grid_step / 2.,
                            x_left_edge + x_grid_step * (xdim - 0.5),
                            num=xdim,
                            dtype=np.float64)
    xs[:] = xs_values
    xs.standard_name = 'projection_x_coordinate'
    xs.coverage_content_type = 'coordinate'
    xs.long_name = 'x'
    xs.units = 'meters'
    xs.axis = 'X'
    # Note: because y-index is at top, min and max are reversed
    xs.valid_min = x_left_edge
    xs.valid_max = x_left_edge + x_grid_step * (xdim + 1)

    ys = ds.createVariable('y', np.float64, ('y',))
    ys_values = np.linspace(y_top_edge + y_grid_step / 2.,
                            y_top_edge + y_grid_step * (ydim - 0.5),
                            num=ydim,
                            dtype=np.float64)
    ys[:] = ys_values
    ys.standard_name = 'projection_y_coordinate'
    ys.coverage_content_type = 'coordinate'
    ys.long_name = 'y'
    ys.units = 'meters'
    ys.axis = 'Y'
    ys.valid_min = y_top_edge + y_grid_step * (ydim + 1)
    ys.valid_max = y_top_edge


    times = ds.createVariable('time', np.float64, ('time',))
    times.standard_name = 'time'
    times.units = 'days since 1970-01-01'
    times.calendar = 'gregorian'
    times.axis = 'T'

    n_times = 10

    first_date = dt.datetime(2018, 1, 1, 12, 0)
    date_increment = dt.timedelta(days=7)
    time_values = []
    for t_index in range(n_times):
        time_values.append(
            first_date + t_index * date_increment)

    times[:] = netCDF4.date2num(time_values,
                                units=times.units,
                                calendar=times.calendar
                               )

    # Define the CRS in GDAL-Create-like manner
    crs_var = ds.createVariable('polar_stereographic', str)
    crs_var.grid_mapping_name = 'polar_stereographic'
    crs_var.straight_vertical_longitude_from_pole = -45.
    crs_var.false_easting = 0.
    crs_var.false_northing = 0.
    crs_var.latitude_of_projection_origin = 90.
    crs_var.standard_parallel = 70.
    crs_var.long_name = 'CRS definition'
    crs_var.longitude_of_prime_meridian = 0.
    crs_var.semi_major_axis = 6378273.
    crs_var.inverse_flattening = 298.279411123064
    crs_var.spatial_ref = 'PROJCS[\"NSIDC Sea Ice Polar Stereographic North\",GEOGCS[\"Unspecified datum based upon the Hughes 1980 ellipsoid\",DATUM[\"Not_specified_based_on_Hughes_1980_ellipsoid\",SPHEROID[\"Hughes 1980\",6378273,298.279411123064,AUTHORITY[\"EPSG\",\"7058\"]],AUTHORITY[\"EPSG\",\"6054\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AUTHORITY[\"EPSG\",\"4054\"]],PROJECTION[\"Polar_Stereographic\"],PARAMETER[\"latitude_of_origin\",70],PARAMETER[\"central_meridian\",-45],PARAMETER[\"scale_factor\",1],PARAMETER[\"false_easting\",0],PARAMETER[\"false_northing\",0],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]],AXIS[\"X\",EAST],AXIS[\"Y\",NORTH],AUTHORITY[\"EPSG\",\"3411\"]]'
    crs_var.GeoTransform = '-3850000 25000 0 5850000 0 -25000'

    conc = ds.createVariable('conc',
                              np.float32,
                              ('time', 'y', 'x',),
                             zlib=True,
                             fill_value=np.nan)
    conc.short_name = 'conc'
    conc.long_name = 'sea_ice_concentration'
    conc.standard_name = 'sea_ice_area_fraction'
    conc.units = 1;
    conc.grid_mapping = 'polar_stereographic'
    conc.flag_values = np.array((251, 253, 254, 255), dtype=np.float32)
    conc.flag_meaning = 'pole_hole coastline land_mask missing_data'

    ydim_data, xdim_data = data_values.shape
    assert xdim == xdim_data
    assert ydim == ydim_data

    for time_index in range(n_times):
        temp_array = np.zeros(data_values.shape, dtype=np.float32)
        temp_array[:] = \
            np.subtract(data_values[:], 10 * (n_times - time_index))

        temp_array[data_values >= 250] = data_values[data_values >= 250]

        conc[time_index, :] = temp_array[:]

    ds.close()
