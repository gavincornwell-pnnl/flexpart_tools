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


def read_partposition(fname):
    """
    This function reads in a partposition output file from a flexwrf simulation
    and returns a matrix of values for each point.
    
    partout columns are structured the following way:
        (1) release; (2) timepoint; (3) longitude; (4) latitude; (5) altitude;
        (6) topography; (7) potential vorticity; (8) water mixing ratio; (9) air density;
        (10) PBL height, m agl; (11) tropopause height, m agl; (12) temperature, K;
        (13) mass for each particle

    """
    import numpy as np
    f = open(fname)
    x = f.readline().split()
    numpart = int(x[1])
    time = int(x[0])
    #    print(time,numpart)
    f.close()
    partout = np.loadtxt(fname, skiprows=1)  # read in file and skip the first line
    if np.size(partout[0]) > 1:  # exception for if file is empty
        partout = partout[:-1, :]  # remove the last column, which is just empty data
    else:
        partout = np.empty([0, 13])
    return partout, time, numpart


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


def calc_srs(fname, grid, delta_t, utot, *args):
    """
    This function calcualtes a source-receptor influence footprint for a given partouput file.
    fname is the full file name for which a SRS should be generated for.
    grid is a two dimensional grid specifying the lon/lat
    delta_t is the time interval that should be used in the calculation of the SRS (units in seconds)
    utot is the total mass released (units in kg)
    *args is an optional input for whether or not to turn the flexpart kernel on (see above function)
    """
    alt1 = 0
    alt2 = 300
    import numpy as np
    if len(args) > 0:
        kernel = args[0]
    else:
        kernel = 0
    temp_part, time, numpart = read_partposition_flexpart(fname)
    srs = np.zeros((len(grid[0]), len(grid[1])))  # set up grid
    grid_diff1 = np.diff(grid[0])[0]
    grid_diff2 = np.diff(grid[1])[0]
    for jidx, j in enumerate(temp_part):
        temp_alt = j[3] + j[4]
        if np.logical_and(temp_alt >= alt1, temp_alt <= alt2):
            q = j[-1]  # for mass concentration output
            if kernel == 1:
                weights = calc_kernel(j[2], j[3], grid)
                if weights[0][0] >= 0:
                    srs[weights[1][0]][weights[2][0]] = srs[weights[1][0]][weights[2][0]] + (
                            q * weights[0][0] * np.abs(delta_t)) / utot
                    srs[weights[1][1]][weights[2][1]] = srs[weights[1][1]][weights[2][1]] + (
                            q * weights[0][1] * np.abs(delta_t)) / utot
                    srs[weights[1][2]][weights[2][2]] = srs[weights[1][2]][weights[2][2]] + (
                            q * weights[0][2] * np.abs(delta_t)) / utot
                    srs[weights[1][3]][weights[2][3]] = srs[weights[1][3]][weights[2][3]] + (
                            q * weights[0][3] * np.abs(delta_t)) / utot
            else:
                idx1 = np.where(np.logical_and(j[2] >= grid[0], j[2] - grid_diff1 <= grid[0]))[0]
                idx2 = np.where(np.logical_and(j[3] >= grid[1], j[3] - grid_diff2 <= grid[1]))[0]
                if np.logical_and(len(idx1) > 0, len(idx2) > 0):
                    idx1 = idx1[0]
                    idx2 = idx2[0]
                    #					print(j[2],grid[0][idx1],j[3],grid[1][idx2])
                    srs[idx1][idx2] = srs[idx1][idx2] + (q * np.abs(delta_t) / utot)
    #				else:
    #					print(j[2],j[3])

    return srs, time



def calc_srs_hrrr(fname, delta_t=3600, utot=1000):
    """
    This function calculates a source-receptor influence footprint for a given partouput file from a flexpart 11
    simulation. The grid corresponds to the NAM12 grid https://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID218).
    fname is the full file name for which a SRS should be generated for.
    delta_t is the time interval that should be used in the calculation of the SRS (units in seconds)
    utot is the total mass released (units in kg)

    returns a srs, which is
    """
    import flexpart as fp
    import numpy as np
    srs = np.zeros((24, 1799, 1059))
    hmix = np.zeros((24, 1799, 1059))
    part_lon, part_lat, part_mass, part_z, part_hmix, time, seconds = read_partposition_fp11(fname)
    for i in np.arange(24):
        tmp_lon = part_lon[:, i]
        tmp_lat = part_lat[:, i]
        tmp_mass = part_mass[:, i]
        tmp_z = part_z[:, i]
        tmp_hmix = part_hmix[:, i]
        mass_grid, hmix_grid, lon, lat, lon_e, lat_e, indices = bin_particles_NAM12(tmp_lon, tmp_lat, tmp_z,
                                                                                    tmp_hmix, tmp_mass)
        srs_tmp = (mass_grid * delta_t / utot)
        srs[i] = srs_tmp
        hmix[i] = hmix_grid
    return srs, hmix, seconds



def calc_srs_nam12(fname, delta_t=3600, utot=1000):
    """
    This function calculates a source-receptor influence footprint for a given partouput file from a flexpart 11
    simulation. The grid corresponds to the NAM12 grid https://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID218).
    fname is the full file name for which a SRS should be generated for.
    delta_t is the time interval that should be used in the calculation of the SRS (units in seconds)
    utot is the total mass released (units in kg)

    returns a srs, which is
    """
    import flexpart as fp
    import numpy as np
    srs = np.zeros((24, 614, 428))
    hmix = np.zeros((24, 614, 428))
    part_lon, part_lat, part_mass, part_z, part_hmix, time, seconds = read_partposition_fp11(fname)
    for i in np.arange(24):
        tmp_lon = part_lon[:, i]
        tmp_lat = part_lat[:, i]
        tmp_mass = part_mass[:, i]
        tmp_z = part_z[:, i]
        tmp_hmix = part_hmix[:, i]
        mass_grid, hmix_grid, lon, lat, lon_e, lat_e, indices = bin_particles_NAM12(tmp_lon, tmp_lat, tmp_z,
                                                                                    tmp_hmix, tmp_mass)
        srs_tmp = (mass_grid * delta_t / utot)
        srs[i] = srs_tmp
        hmix[i] = hmix_grid
    return srs, hmix, seconds


def bin_particles_NAM12(part_lon, part_lat, part_z, hmix, part_mass):
    '''
    This function is intended to take the geodetic data from particle from FLEXPART simulations and transpose it on to the NAM12 grid.
    Particle locations are transformed using pyproj and then the distance to the grid cell corners are minimized to determine what cell the particles are in.
    '''
    from pyproj import Proj
    import numpy as np
    from scipy.stats import binned_statistic_2d as bs2d

    ## set up output grid
    lon, lat, lon_e, lat_e = NAM12_grid() # XX_e denotes edges

    # only select particles that are within the boundary layer
    idx = np.where(part_z < hmix)

    # get particle grid coordinate distances
    part_x, part_y = lcc_proj(part_lon[idx], part_lat[idx], inverse=False)

    # bin particles by their location, and sum their mass
    mass_grid, x_edges, y_edges, indices = bs2d(x=part_x, y=part_y, values=part_mass[idx], statistic='sum',
                                                bins=(x_e, y_e))

    # bin particles by their location, and get their average mixing height
    hmix_grid, _, _, _ = bs2d(x=part_x, y=part_y, values=hmix[idx], statistic='mean',
                              bins=(x_e, y_e))
    return mass_grid, hmix_grid, lon, lat, lon_e, lat_e, indices


def bin_particles_HRRR(part_lon, part_lat, part_z, hmix, part_mass):
    '''
    This function is intended to take the geodetic data from particle from FLEXPART simulations and transpose it on to the NAM12 grid.
    Particle locations are transformed using pyproj and then the distance to the grid cell corners are minimized to determine what cell the particles are in.
    '''
    from pyproj import Proj
    import numpy as np
    from scipy.stats import binned_statistic_2d as bs2d

    ## set up output grid
    lon, lat, lon_e, lat_e = HRRR_grid() # XX_e denotes edges

    # only select particles that are within the boundary layer
    idx = np.where(part_z < hmix)

    # get particle grid coordinate distances
    part_x, part_y = lcc_proj(part_lon[idx], part_lat[idx], inverse=False)

    # bin particles by their location, and sum their mass
    mass_grid, x_edges, y_edges, indices = bs2d(x=part_x, y=part_y, values=part_mass[idx], statistic='sum',
                                                bins=(x_e, y_e))

    # bin particles by their location, and get their average mixing height
    hmix_grid, _, _, _ = bs2d(x=part_x, y=part_y, values=hmix[idx], statistic='mean',
                              bins=(x_e, y_e))


def NAM12_grid():
    from pyproj import Proj
    import numpy as np

    ## set up output grids
    dx, dy = (12191.0, 12191.0)
    nx, ny = (614, 428)

    # hardcode parameters and projection information
    lcc_proj = Proj(proj='lcc', lat_1=25, lat_2=25, lat_0=25, lon_0=265,
                    R=6371229)  # projection from the NAM12 grid, information taken from grib files
    lon1 = 226.541  # from grib files
    lat1 = 12.19  # from grib files

    # first makes arrays that correspond to the bin edges, which will be used when binning
    # this is needed because WRF coordinates correspond to the lower left corner of the grid cell
    llcrnrx, llcrnry = lcc_proj(lon1, lat1)
    x_e = llcrnrx + np.arange(nx + 1) * (dx)  # has to extend out an extra grid cell because it is the edge
    y_e = llcrnry + np.arange(ny + 1) * (dy)

    x_grid, y_grid = np.meshgrid(x_e, y_e)
    if x_e.shape != x_grid.shape[0]:
        x_grid = x_grid.T
        y_grid = y_grid.T
    lon_e, lat_e = lcc_proj(x_grid, y_grid, inverse=True)

    # second make a lon/lat set of arrays that corresponds to the lower left corners of the grid cells
    x = x_e[0:-1]
    y = y_e[0:-1]
    x_grid, y_grid = np.meshgrid(x, y)
    if x.shape != x_grid.shape[0]:
        x_grid = x_grid.T
        y_grid = y_grid.T
    lon, lat = lcc_proj(x_grid, y_grid, inverse=True)
    return lon, lat, lon_e, lat_e


def HRRR_grid():
    from pyproj import Proj
    import numpy as np

    ## set up output grids
    dx, dy = (3000.0, 3000.0)
    nx, ny = (1799, 1059)

    # hardcode parameters and projection information
    lcc_proj = Proj(proj='lcc', lat_1=38.5, lat_2=38.5, lat_0=38.5, lon_0=262.5,
                    R=6371229)  # projection from the NAM12 grid, information taken from grib files
    lon1 = 237.280472  # from grib files
    lat1 = 21.138123  # from grib files

    # first makes arrays that correspond to the bin edges, which will be used when binning
    # this is needed because WRF coordinates correspond to the lower left corner of the grid cell
    llcrnrx, llcrnry = lcc_proj(lon1, lat1)
    x_e = llcrnrx + np.arange(nx + 1) * (dx)  # has to extend out an extra grid cell because it is the edge
    y_e = llcrnry + np.arange(ny + 1) * (dy)

    x_grid, y_grid = np.meshgrid(x_e, y_e)
    if x_e.shape != x_grid.shape[0]:
        x_grid = x_grid.T
        y_grid = y_grid.T
    lon_e, lat_e = lcc_proj(x_grid, y_grid, inverse=True)

    # second make a lon/lat set of arrays that corresponds to the lower left corners of the grid cells
    x = x_e[0:-1]
    y = y_e[0:-1]
    x_grid, y_grid = np.meshgrid(x, y)
    if x.shape != x_grid.shape[0]:
        x_grid = x_grid.T
        y_grid = y_grid.T
    lon, lat = lcc_proj(x_grid, y_grid, inverse=True)
    return lon, lat, lon_e, lat_e
