# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 10:19:11 2019
This library contains functions for the analysis of simulations from FLEXWRF.
Additional scripts are included for pre-processing of data as well.
@author: corn062
"""


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


def read_partposition_flexpart(fname):
    """
    This function reads in a partposition output file from a flexwrf simulation
    and returns a matrix of values for each point.
    
    partout columns are structured the following way:
        (1) release; (2) timepoint particle was released; (3) longitude; (4) latitude; (5) altitude;
        (6) topography; (7) potential vorticity; (8) water mixing ratio; (9) air density;
        (10) PBL height, m agl; (11) tropopause height, m agl; (12) temperature, K;
        (13) mass for each particle

    """
    import numpy as np
    import matplotlib.dates as dt
    f = open(fname)
    x = f.readline().split()
    numpart = int(x[0])
    time = dt.num2date(dt.datestr2num(fname.split('_')[-1]))
    #    print(time,numpart)
    f.close()
    partout = np.loadtxt(fname, skiprows=1)  # read in file and skip the first line
    if np.size(partout[0]) > 1:  # exception for if file is empty
        partout = partout[:-1, :]  # remove the last column, which is just empty data
    else:
        partout = np.empty([0, 13])
    return partout, time, numpart


def read_partposition_fp11(fname):
    '''
    :param fname: name of file to be read
    :return: part_lon
    :return: part_lat
    :return: part_mass
    :return: part_z
    :return: part_hmix
    '''
    import netCDF4 as nc
    import numpy as np
    import matplotlib.dates as dt

    ds = nc.Dataset(fname)
    part_lon = ds.variables['lon_av'][:, :]
    part_lat = ds.variables['lat_av'][:, :]
    part_mass = ds.variables['m_av001'][:, :]
    part_z = ds.variables['z_av'][:, :]
    part_hmix = ds.variables['hmix_av'][:, :]
    time = ds.variables['time'][:]
    delta_t = 3600
    if (len(time)) > 1:
        delta_t = np.diff(time)[0]
    tmp = fname.split('_')[-1].split('.')[0]
    if tmp == 'init':
        tmp = fname.split('_')[-2]
    year = int(tmp[0:4])
    month = int(tmp[4:6])
    day = int(tmp[6:8])
    hour = int(tmp[8:10])
    minute = int(tmp[10:12])
    second = int(tmp[12:14])
    basetime = dt.date2num(dt.datetime.datetime(year, month, day, hour, minute, second)) * 24 * 3600
    matlab_times = basetime + np.arange(0, 24) * delta_t
    ds.close()
    return part_lon, part_lat, part_mass, part_z, part_hmix, time, matlab_times


def calc_kernel(lon, lat, grid):
    """
    This script replicates the kernel used by FLEXPART/FLEXWRF in order to spread particle mass out
    over multiple cells.
    """
    import numpy as np
    lonbins = grid[0]  # longitude grid
    latbins = grid[1]  # latitude grid
    grid_diff1 = np.diff(grid[0])[0]
    grid_diff2 = np.diff(grid[1])[0]
    # lower- and left-most coor grid cell coordinates for lat/lon
    idx1 = np.where(np.logical_and(lon >= lonbins, lon - grid_diff1 <= lonbins))[0]
    idx2 = np.where(np.logical_and(lat >= latbins, lat - grid_diff2 <= latbins))[0]
    if np.logical_and(idx1.size > 0, idx2.size > 0):
        lon_idx1 = idx1[0]
        lat_idx1 = idx2[0]
        # determine distance to edge of grid cell so that direction of weighting can be determined
        dx = np.round(np.diff(grid[0])[0], 6)
        dy = np.round(np.diff(grid[1])[0], 6)
        ddx = round((lon - lonbins[lon_idx1]) / dx, 6)
        ddy = round((lat - latbins[lat_idx1]) / dy, 6)
        # get weighting for different cells
        if ddx > 0.5:  # if closer to rightside of cell, shift upward
            lon_idx2 = lon_idx1 + 1
            if lon_idx2 >= len(lonbins):
                lon_idx2 = 0
            wx = 1.5 - ddx
        else:  # closer to leftside of cell, shift down
            lon_idx2 = lon_idx1 - 1
            wx = 0.5 + ddx
        if ddy > 0.5:  # if closer to top of cell, shiftupward
            lat_idx2 = lat_idx1 + 1
            if lat_idx2 >= len(latbins):
                lat_idx2 = 0
            wy = 1.5 - ddy
        else:  # closer to bottom of cell, shift down
            lat_idx2 = lat_idx1 - 1
            wy = 0.5 + ddy
        # graph out
        weights = [round(wx * wy, 6), round(wx * (1 - wy), 6), round((1 - wx) * wy, 6), round((1 - wx) * (1 - wy), 6)]
        lonidx_out = [lon_idx1, lon_idx1, lon_idx2, lon_idx2]
        latidx_out = [lat_idx1, lat_idx2, lat_idx1, lat_idx2]
    else:
        weights = np.ones(4) * -1
        lonidx_out = np.ones(4) * -1
        latidx_out = np.ones(4) * -1
    return weights, lonidx_out, latidx_out


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


def calc_srs_nam12(fname, delta_t, utot):
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
    utot = 1000000
    delta_t = 3600
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


def calc_srs_wrf(fname, x_grid, y_grid, delta_t, utot):
    """
    This function calculates a source-receptor influence footprint for a given partouput file from a flexpart-wrf
    simulation. The grid corresponds to the HRRR grid
    fname is the full file name for which a SRS should be generated for.
    delta_t is the time interval that should be used in the calculation of the SRS (units in seconds)
    utot is the total mass released (units in kg)
    """
    import flexpart as fp
    import numpy as np

    srs = np.zeros((x_grid.shape[0], y_grid.shape[0]))
    hmix = np.zeros(np.shape(srs))
    partout, seconds, numpart = read_partposition(fname)
    part_x = partout[:, 2]
    part_y = partout[:, 3]
    part_z = partout[:, 4]
    part_hmix = partout[:, 5]
    part_mass = partout[:, -1]
    mass_grid, hmix_grid, indices = bin_particles_WRF(part_x, part_y, part_z, part_hmix, part_mass)
    srs_tmp = (mass_grid * delta_t / utot)
    srs = srs_tmp
    hmix = hmix_grid
    return srs, hmix, seconds


def calc_srs_fp11(fname, x_grid, y_grid, delta_t, utot):
    """
    This function calcualtes a source-receptor influence footprint for a given partouput.nc file
    fname is the full file name for which a SRS should be generated for.
    x_grid and y_grid are both 1-d arrays specifying the lon/lat
    delta_t is the time interval that should be used in the calculation of the SRS (units in seconds)
    utot is the total mass released (units in kg)
    *args is an optional input for whether or not to turn the flexpart kernel on (see above function)
    """
    import numpy as np

    part_lon, part_lat, part_mass, part_z, part_hmix, time, matlab_times = read_partposition_fp11(fname)
    srs = np.zeros((len(matlab_times), len(x_grid), len(y_grid)))  # set up grid
    hmix = np.zeros((len(matlab_times), len(x_grid), len(y_grid)))  # set up grid
    for iidx, i in enumerate(matlab_times):
        srs_tmp = np.zeros((len(x_grid), len(y_grid)))
        temp_lon = part_lon[:, iidx]
        temp_lat = part_lat[:, iidx]
        temp_alt = part_z[:, iidx]
        temp_hmix = part_hmix[:, iidx]
        temp_mass = part_mass[:, iidx]
        mass_grid, hmix_grid, indices = bin_particles_FP11(temp_lon, temp_lat, temp_alt, temp_hmix, temp_mass, x_grid,
                                                           y_grid)
        srs_tmp = (mass_grid * delta_t / utot)
        srs[iidx] = srs_tmp
        hmix[iidx] = hmix_grid


def bin_particles_NAM12(part_lon, part_lat, part_z, hmix, part_mass):
    '''
    This function is intended to take the geodetic data from particle from FLEXPART simulations and transpose it on to the NAM12 grid.
    Particle locations are transformed using pyproj, and the
    '''
    from pyproj import Proj
    import numpy as np
    from scipy.stats import binned_statistic_2d as bs2d

    # hardcode parameters and projection information
    lcc_proj = Proj(proj='lcc', lat_1=25, lat_2=25, lat_0=25, lon_0=265,
                    R=6371229)  # projection from the NAM12 grid, information taken from grib files
    lon, lat, lon_e, lat_e = NAM12_grid()  # lon/lat values

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
    Particle locations are transformed using pyproj, and the
    '''
    from pyproj import Proj
    import numpy as np
    from scipy.stats import binned_statistic_2d as bs2d

    # hardcode parameters and projection information
    lcc_proj = Proj(proj='lcc', lat_1=25, lat_2=25, lat_0=25, lon_0=265,
                    R=6371229)  # projection from the NAM12 grid, information taken from grib files
    lon, lat, lon_e, lat_e = HRRR_grid()  # lon/lat values

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


def bin_particles_WRF(part_lon, part_lat, part_z, hmix, part_mass):
    '''
    This function is intended to take the geodetic data from particle from FLEXPART simulations and transpose it on to WRF grid.
    '''
    import numpy as np
    from scipy.stats import binned_statistic_2d as bs2d

    # grid information
    dx, dy = (3000, 3000)
    nx, ny = (1799, 1059)
    x_e = np.arange(nx + 1) * (dx)  # has to extend out an extra grid cell because it is the edge
    y_e = np.arange(ny + 1) * (dy)

    # only select particles that are within the boundary layer
    idx = np.where(part_z < hmix)

    # get particle grid coordinate distances
    part_x, part_y = (part_lon[idx], part_lat[idx])

    # bin particles by their location, and sum their mass
    mass_grid, x_edges, y_edges, indices = bs2d(x=part_x, y=part_y, values=part_mass[idx], statistic='sum',
                                                bins=(x_e, y_e))

    # bin particles by their location, and get their average mixing height
    hmix_grid, _, _, _ = bs2d(x=part_x, y=part_y, values=hmix[idx], statistic='mean',
                              bins=(x_e, y_e))
    return mass_grid, hmix_grid, indices


def bin_particles_FP11(part_lon, part_lat, part_z, hmix, part_mass, x_grid, y_grid):
    '''
    This function is intended to take bin particle data from FLEXPART simulations onto a user-specified grid.
    Data coordinates are in lon,lat
    part_z should be m AGL
    x_grid and y_grid are both 1-d arrays of lon/lat
    '''
    import numpy as np
    from scipy.stats import binned_statistic_2d as bs2d

    # only select particles that are within the boundary layer
    idx = np.where(part_z < hmix)

    # get particle grid coordinate distances
    part_x, part_y = (part_lon[idx], part_lat[idx])

    # bin particles by their location, and sum their mass
    mass_grid, x_edges, y_edges, indices = bs2d(x=part_x, y=part_y, values=part_mass[idx], statistic='sum',
                                                bins=(x_grid, y_grid))

    # bin particles by their location, and get their average mixing height
    hmix_grid, _, _, _ = bs2d(x=part_x, y=part_y, values=hmix[idx], statistic='mean',
                              bins=(x_grid, y_grid))
    return mass_grid, hmix_grid, indices


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

    # hardcoded parameters and projection information
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


def calculate_point_HRRR(lon, lat):
    from pyproj import Proj
    import numpy as np

    ## set up output grids
    dx, dy = (3000.0, 3000.0)
    nx, ny = (1799, 1059)

    # hardcoded parameters and projection information
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
