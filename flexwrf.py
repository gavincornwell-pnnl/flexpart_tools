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
    new_x = np.round(new_x + np.abs(llcrnrx))
    new_y = np.round(new_y + np.abs(llcrnry))
    return new_x, new_y


def extract_3d_var(grbs, varname):
    '''
    This function extracts 3D data from a grib file.
    :param fname: filename of grib file
    :param varname: variable to be extracted, must be a 3-d variable
    :return data_corr: data arranged from 1000 to 50 hPa
    :return hPa: hPa of data
    :return lon: longitude
    :return lat: latitude
    '''
    import numpy as np
    import pygrib as pg
    
    grbs.seek(0)
    msg = grbs.read(1)[0]
    nx, ny = msg.Nx, msg.Ny
    msg = grbs.select(shortName=varname, typeOfLevel='hybrid')
    datacat = np.zeros((len(msg), ny, nx))
    levelcat = np.zeros(len(msg))
    for iidx, item in enumerate(msg):
        data, lat, lon = item.data()
        data = (data  * item.scaleValuesBy) + item.offsetValuesBy
        level = item.level
        datacat[iidx, :, :] = data
        levelcat[iidx] = level
    levels = levelcat
    data = datacat[:]
    return data, levels, lon, lat


def bilinear_rect(lon, lat, F, xp, yp, clip=True):
    """
    Bilinear interpolation of F on a rectilinear grid (lat, lon).

    lon: (nx,) increasing
    lat: (ny,) increasing
    F:   (ny, nx) or (..., ny, nx)
    xp:  (npart,) lon of query points
    yp:  (npart,) lat of query points
    """
    import numpy as np
    nx = lon.size
    ny = lat.size
    npart = xp.size

    # 1) bracket indices (hi is first index where coord >= point)
    ix_hi = np.searchsorted(lon, xp, side="left")
    iy_hi = np.searchsorted(lat, yp, side="left")
    ix_lo = ix_hi - 1
    iy_lo = iy_hi - 1

    # 2) handle boundaries
    if clip:
        ix_lo = np.clip(ix_lo, 0, nx - 2)
        iy_lo = np.clip(iy_lo, 0, ny - 2)
        ix_hi = ix_lo + 1
        iy_hi = iy_lo + 1
    else:
        # you can decide your own behavior for out-of-domain points
        pass

    # 3) fractional distance in cell [lo, hi]
    x0 = lon[ix_lo]; x1 = lon[ix_hi]
    y0 = lat[iy_lo]; y1 = lat[iy_hi]

    tx = (xp - x0) / (x1 - x0)
    ty = (yp - y0) / (y1 - y0)

    # keep weights sane if you clipped indices
    tx = np.clip(tx, 0.0, 1.0)
    ty = np.clip(ty, 0.0, 1.0)

    # 4) gather four corners
    # F00 = (y0, x0), F10 = (y0, x1), F01 = (y1, x0), F11 = (y1, x1)
    # Works for F of shape (ny, nx) and also (..., ny, nx)
    F00 = F[..., iy_lo, ix_lo]
    F10 = F[..., iy_lo, ix_hi]
    F01 = F[..., iy_hi, ix_lo]
    F11 = F[..., iy_hi, ix_hi]

    # 5) bilinear combination
    w00 = (1 - tx) * (1 - ty)
    w10 = tx * (1 - ty)
    w01 = (1 - tx) * ty
    w11 = tx * ty

    # If F has leading dims, broadcast weights to match
    # (npart,) -> (1,...,1,npart) would be needed only if you keep npart as last dim.
    # Here F[..., iy, ix] returns shape (..., npart), so weights broadcast fine.
    out = w00 * F00 + w10 * F10 + w01 * F01 + w11 * F11

    return out, (ix_lo, ix_hi, iy_lo, iy_hi), (tx, ty)

