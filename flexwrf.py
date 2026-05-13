# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 10:19:11 2019
This library contains functions for writing, processing, and analyzing flexpart-wrf simulations.
@author: corn062
"""


def modify_lasso(fname, gridcell_width, numgridcells):
    """
    This function reads in a LASSO output file and modifies the latitude/longitude
    values such that they are not all the same. Without this, the simulations
    cannot run properly.
    
    fname is the filename of the ouptut to be modified
    gridcell_width is the width of each grid cell in meters
    numgridcells is the number of gridcells in the simulation
    """
    import numpy as np, netCDF4 as nc, os
    files = os.listdir('./')
    files_nc = [i for i in files if i.startswith('wrfout')]
    for k in np.arange(0, len(files_nc)):
        x = nc.Dataset(fname, 'r+')
        # formula to calculate the distance between degrees of longitude
        # 1 degree of longitude = cosein(lat-radians) * length of degree at equator
        long_rad = (36.6 / (180 / np.pi))  # latitude in radians
        long_1deg = np.cos(long_rad) * 111.32  # km
        long_1deg = long_1deg * 1000  # m
        long_spacing = gridcell_width / long_1deg
        new_longitude = -97.5 + np.round(np.arange(0, numgridcells) * long_spacing, 3)
        lat_1deg = 111  # km
        lat_spacing = gridcell_width / (lat_1deg * 1000)
        new_latitude = 36.6 + np.round(np.arange(0, numgridcells) * lat_spacing, 3)
        lat = x['XLAT']
        lon = x['XLONG']
        for j in np.arange(0, 6):
            for i in np.arange(0, 250):
                lat[j, :, i] = new_latitude
                lon[j, i, :] = new_longitude
        x.close()
        print(fname + ' finished processing, chief')


def read_turbulence(filename, numpart):
    """
    This function reads in a turboutput file from the modified FLEXWRF code
    and returns a matrix of values for each point.
    
    output matrix columns are structured the following way:
        (1) simulation time advancef.90 is entered (in minutes); (2) cycle # through advance.f90;
        (3) delta z (m); (4) particle altitude (m)
    """
    import numpy as np
    xout = np.loadtxt(filename)  # read in data
    xout_len = np.shape(xout)[0]
    # next several lines are to separate out the data for separate particles
    dff = np.diff(xout[:, 1])  #
    tempidx = np.where(dff < 0)
    dff_idx = np.zeros(len(tempidx[0]) + 2, dtype=int)
    dff_idx[1:-1] = tempidx[0] + 1
    dff_idx[-1] = xout_len
    particledataout = [[]] * numpart  # pre-declare variable
    for i in np.arange(0, numpart):
        particledataout[i] = np.empty([0, 5], dtype=float)
    count = 0  # counter for particles
    for i in np.arange(0, len(dff_idx) - 1):
        blocksize = dff_idx[i + 1] - dff_idx[i]
        temp = np.empty([blocksize, 5])
        temp = xout[dff_idx[i]:dff_idx[i + 1], [1, 2, 3, 5, 6]]
        temp = xout[dff_idx[i]:dff_idx[i + 1], [1, 2, 3, 5, 6]]
        particledataout[count] = np.vstack((particledataout[count], temp))
        count = count + 1
        if count == numpart:
            count = 0
    return particledataout


def ll_to_wrf_xy_HRRR(new_lon, new_lat):
    # script to calculate the WRF coordinates for a WRF simulation, for a given lon/lat point
    from pyproj import Proj
    import numpy as np
    
    # new_lon, new_lat = -97.485, 36.605
    
    # hardcoded parameters and projection information
    lcc_proj = Proj(proj='lcc', lat_1=38.5, lat_2=38.5, lat_0=38.5, lon_0=262.5,
                    R=6371229)  # projection from the NAM12 grid, information taken from grib files
    lon1 = 237.280472  # from grib files
    lat1 = 21.138123  # from grib files
    
    # first makes arrays that correspond to the bin edges, which will be used when binning
    # this is needed because WRF coordinates correspond to the lower left corner of the grid cell
    llcrnrx, llcrnry = lcc_proj(lon1, lat1)
    new_x, new_y = lcc_proj(new_lon, new_lat)
    new_x = int(np.round(new_x + np.abs(llcrnrx)))
    new_y = int(np.round(new_y + np.abs(llcrnry)))
    return new_x, new_y

