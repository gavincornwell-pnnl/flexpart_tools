import numpy as np
import netCDF4 as nc
import flexpart as fp
import flexwrf as fw
import os
import sys
import matplotlib.dates as dt
from scipy.stats import binned_statistic_2d as bs2d
import pygrib as pg
import matplotlib.pyplot as plt

def HRRR_trajectory_processing(folder):
    '''
    1. identify the variables from grib to correlate to
    2. create placeholder variables
    3. read in particle data, and extract location information, as well as mass
    4. locate relevant metfile in grib data
    5. get variable data at each particle location
    6. open save file
    7. save data
    '''
    ### read in variables to be assigned to
    var_name = r'/people/corn062/flexpart_tools/Surface_variables.txt'
    fid = open(var_name, 'r')
    lines_s = np.sort(fid.readlines())
    fid.close()
    var_name = r'/people/corn062/flexpart_tools/AGL_variables.txt'
    fid = open(var_name, 'r')
    lines_a = np.sort(fid.readlines())
    fid.close()
    var_name = r'/people/corn062/flexpart_tools/Hybrid_variables.txt'
    fid = open(var_name, 'r')
    lines_h = np.sort(fid.readlines())
    fid.close()

    ### create placeholder variables for metadata writing
    units_s = [[]] * len(lines_s)
    shortnames_s = [[]] * len(lines_s)
    names_s = [[]] * len(lines_s)
    units_a = [[]] * len(lines_a)
    shortnames_a = [[]] * len(lines_a)
    names_a = [[]] * len(lines_a)
    units_h = [[]] * len(lines_h)
    shortnames_h = [[]] * len(lines_h)
    names_h = [[]] * len(lines_h)

    #### get filelist from the folder to be analyzed
    # folder = '/rcfs/scratch/corn062/AGINSGP/ncep/20220410_010000_050'
    # folder = '/rcfs/scratch/corn062/AGINSGP/HRRR/AGINSGP.220408000000'
    print(folder)
    fnamelist = os.listdir(folder)
    # fnamelist = [i for i in fnamelist if i.endswith('.nc')]
    print(fnamelist)
    fnamelist = np.sort([i for i in fnamelist if i.startswith('part')])[::-1]
    gribfolder = '/rcfs/scratch/corn062/BNF/met/HRRR/grib'

    ### variables that are used in defining spatial and particle arrays
    numPart = 100000
    partarray = np.arange(1, numPart + 1)
    dx = dy = 3000
    nx = 1799
    ny = 1059
    x_arr = np.arange(0, nx) * dx
    y_arr = np.arange(0, ny) * dy

    ### Create variables that will be written
    part_surface = np.ones((len(lines_s), numPart, len(fnamelist))) * -9999
    part_hybrid = np.ones((len(lines_h), numPart, len(fnamelist))) * -9999
    part_agl = np.ones((len(lines_a), numPart, len(fnamelist))) * -9999
    part_mass = np.ones((numPart, len(fnamelist))) * -9999
    part_hmix = np.ones((numPart, len(fnamelist))) * -9999
    part_alt = np.ones((numPart, len(fnamelist))) * -9999
    times = np.ones(len(fnamelist)) * -9999

    ## get individual particle output file and read in its data
    for fidx, f in enumerate(fnamelist[:]):
        print(fidx / len(fnamelist))
        fid = open(os.path.join(folder, f), 'r')
        fname = os.path.join(folder, fnamelist[fidx])
        # part_lon, part_lat, part_mass, part_z, part_hmix, time, matlab_times = fp.read_partposition_fp11(fname)
        partout, time, _ = fp.read_partposition_flexpart(fname)
        pidx = np.where(np.isin(partout[:, 0], partarray))[0]
        part_x_ = partout[:, 2]
        part_y_ = partout[:, 3]
        part_z_ = partout[:, 4]
        part_hmix_ = partout[:, 5]
        part_mass_ = partout[:, 12]
        date = time
        times[fidx] = dt.date2num(date)
        part_mass[:,fidx] = part_mass_
        part_hmix[:,fidx] = part_hmix_
        part_alt[:,fidx] = part_z_

        ## get grib file and read in data from it
        grib_fname = os.path.join(gribfolder, date.strftime('hrrr_%Y%m%d_%H00.grib2'))
        grib = pg.open(grib_fname)

        ## get altitude of levels by taking the geopotential and subtracting the terrain (orography)
        gh, _, _, _ = fw.extract_3d_var(grib, 'gh')
        m = grib.select(shortName='orog')[0]
        oro, _, _ = m.data()
        gh = gh - oro  # turn GH to above AGL

        ### bilinearly interpolate altitude to particle position
        part_gh, _, _ = fw.bilinear_rect(x_arr, y_arr, gh, part_x_, part_y_)
        idx = np.arange(part_z_.size)  # place holder index variable
        part_gh_T = part_gh.T  # (npart, nlev)

        # index of first level >= parcel altitude
        k_hi = np.array([np.searchsorted(part_gh_T[i], part_z_[i], side="left") for i in range(part_gh_T.shape[0])])
        k_lo = k_hi - 1  # convert to "below" index
        # clamp to valid range
        k_lo = np.clip(k_lo, 0, part_gh_T.shape[1] - 1)
        k_hi = np.clip(k_hi, 0, part_gh_T.shape[1] - 1)
        ## weight altitude and find weights for low and high indices
        z_lo = part_gh[k_lo, idx]
        z_hi = part_gh[k_hi, idx]
        den = (z_hi - z_lo)
        w = (part_z_ - z_lo) / den  ## weight
        w = np.where(den != 0, w, 0.0)

        ### loop through all the 3D variables
        for hidx, h in enumerate(lines_h[:]):
            shortnames_h[hidx] = shortName = h.split(',')[1].strip()
            print('Hybrid, %s' % shortName)
            msg = grib.select(shortName=shortName, typeOfLevel='hybrid')[0]
            units_h[hidx] = msg.units
            names_h[hidx] = msg.name
            data, _, _, _ = fw.extract_3d_var(grib, shortName)
            #    bilin_data = fw.bilinear_rect(x_arr, y_arr, data, part_x_, part_y_)
            shortnames_h[hidx] = '3d_%s' % shortName
            bilin_data, _, _ = fw.bilinear_rect(x_arr, y_arr, data, part_x_, part_y_)
            data_lo = bilin_data[k_lo, idx]
            data_hi = bilin_data[k_hi, idx]
            data_interp = data_lo + w * (data_hi - data_lo)  # same as (1-w)*v_lo + w*v_hi
            part_hybrid[hidx, :, fidx] = data_interp

        # loop through surface variables
        for sidx, s in enumerate(lines_s[:]):
            shortName = s.split(',')[1].strip()
            msg = grib.select(shortName=shortName, typeOfLevel='surface')[0]
            data, _, _ = msg.data()
            data = (data * msg.scaleValuesBy) + msg.offsetValuesBy
            units_s[sidx] = msg.units
            names_s[sidx] = msg.name
            shortnames_s[sidx] = 'surface_%s' % shortName
            bilin_data, _, _ = fw.bilinear_rect(x_arr, y_arr, data, part_x_, part_y_)
            part_surface[sidx, :, fidx] = bilin_data

        # loop through above ground levels
        for aidx, a in enumerate(lines_a[:]):
            shortName = a.split(',')[1].strip()
            msg = grib.select(shortName=shortName, typeOfLevel='heightAboveGround')[0]
            data, _, _ = msg.data()
            data = (data * msg.scaleValuesBy) + msg.offsetValuesBy
            units_a[aidx] = msg.units
            names_a[aidx] = msg.name
            shortnames_a[aidx] = shortName
            bilin_data, _, _ = fw.bilinear_rect(x_arr, y_arr, data, part_x_, part_y_)
            part_agl[aidx, :, fidx] = bilin_data

    ######### create output file
    basename = os.path.basename(folder)
    if basename.endswith('.nc'):
        basename = basename.split('_')[1]
    dirname = os.path.dirname(fname)
    sname = os.path.join(dirname, 'HRRR_var_tracking.%s.nc' % basename)
    ## open file and create dimensions
    ds_out = nc.Dataset(sname, 'w')
    time_dim = ds_out.createDimension(dimname='time', size=np.shape(part_agl)[2])
    particle_dim = ds_out.createDimension(dimname='particle', size=np.shape(part_agl)[1])

    ##### first create standard particle variables and load attributes and data
    ## create variables
    time_var = ds_out.createVariable(varname='time', datatype='f8', dimensions='time')
    part_mass_var = ds_out.createVariable(varname='part_mass', datatype='f4', dimensions=('particle', 'time'))
    part_alt_var = ds_out.createVariable(varname='part_alt', datatype='f4', dimensions=('particle', 'time'))
    hmix_var = ds_out.createVariable(varname='part_hmix', datatype='f4', dimensions=('particle', 'time'))
    # write attributes
    time_var.units = 'matplotlib date number'
    time_var.setncattr_string('name', 'datenumber')
    part_mass_var.units = 'kg'
    part_mass_var.setncattr_string('name', 'Particle mass')
    part_alt_var.units = 'm'
    part_alt_var.setncattr_string('name', 'Particle altitude AGL')
    hmix_var.units = 'm'
    hmix_var.setncattr_string('name', 'Mixing layer height AGL')
    # write data
    time_var[:] = times
    part_mass_var[:] = part_mass
    part_alt_var[:] = part_alt
    hmix_var[:] = part_hmix

    ##### next step through all surface variables: create variable, write attributes and load data
    for sidx, s in enumerate(lines_s[:]):
        tmp_var = ds_out.createVariable(varname=shortnames_s[sidx], datatype='f4', dimensions=('particle', 'time'))
        tmp_var.units = units_s[sidx]
        tmp_var.setncattr_string('name', names_s[sidx])
        tmp_var[:] = part_surface[sidx]

    ##### next step through all AGL variables: create variable, write attributes and load data
    for aidx, a in enumerate(lines_a[:]):
        tmp_var = ds_out.createVariable(varname=shortnames_a[aidx], datatype='f4', dimensions=('particle', 'time'))
        tmp_var.units = units_a[aidx]
        tmp_var.setncattr_string('name', names_a[aidx])
        tmp_var[:] = part_agl[aidx]

    ##### next step through all 3d variables: create variable, write attributes and load data
    for hidx, a in enumerate(lines_h[:]):
        tmp_var = ds_out.createVariable(varname=shortnames_h[hidx], datatype='f4', dimensions=('particle', 'time'))
        tmp_var.units = units_h[hidx]
        tmp_var.setncattr_string('name', names_h[hidx])
        tmp_var[:] = part_hybrid[hidx]

    ds_out.close()
if __name__ == "__main__":
    folder = sys.argv[1]
    print(folder)
    HRRR_trajectory_processing(folder)
