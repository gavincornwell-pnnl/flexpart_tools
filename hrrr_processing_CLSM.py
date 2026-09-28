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


# def trajectory_processing(folder):
'''
1. identify the variables from grib to correlate to
2. create placeholder variables
3. read in particle data, and extract location information, as well as mass
4. locate relevant metfile in grib data
5. get variable data at each particle location
6. open save file
7. save data
'''

var_name = r'/people/corn062/flexpart_tools/surface_variables.txt'
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


#### get filelist from the folder to be analyzed
folder = '/rcfs/scratch/corn062/AGINSGP/ncep/20220410_010000_050'
fnamelist = os.listdir(folder)
fnamelist = [i for i in fnamelist if i.endswith('.nc')]
fnamelist = np.sort([i for i in fnamelist if i.startswith('partoutput')])[::-1]
gribfolder = '/rcfs/scratch/corn062/BNF/met/HRRR/grib'

## get individual particle output file and read in its data
fidx = 0
fname = os.path.join(folder, fnamelist[fidx])
part_lon, part_lat, part_mass, part_z, part_hmix, time, matlab_times = fp.read_partposition_fp11(fname)
part_x, part_y = fw.ll_to_wrf_xy_HRRR(part_lon,part_lat)
days = np.round(matlab_times / (24 * 3600))
dates = dt.num2date(days)

### Create variables that will be written
numPart = 100000
dx = dy = 3000
nx = 1799
ny = 1059
x_arr = np.arange(0, nx) * dx
y_arr = np.arange(0, ny) * dy
part_surface = np.ones((len(lines_s), numPart, len(days))) * -9999
part_hybrid = np.ones((len(lines_h), numPart, len(days))) * -9999
part_agl = np.ones((len(lines_a), numPart, len(days))) * -9999

didx = 0
date = dates[didx]
grib_fname = os.path.join(gribfolder,date.strftime('hrrr_%Y%m%d_%H%M.grib2'))
grib = pg.open(grib_fname)

part_x_ = part_x[:,didx]
part_y_ = part_y[:,didx]
part_z_ = part_z[:,didx]
part_mass_ = part_mass[:,didx]
part_hmix_ = part_hmix[:,didx]

idx_x = np.digitize(part_x_, x_arr) - 1
idx_x[idx_x < 0] = 0
idx_y = np.digitize(part_y, y_arr) - 1
idx_y[idx_y < 0] = 0

units_s = [[]] * len(lines_s)
shortnames_s = [[]] * len(lines_s)
names_s = [[]] * len(lines_s)

units_a = [[]] * len(lines_a)
shortnames_a = [[]] * len(lines_a)
names_a = [[]] * len(lines_a)

units_h = [[]] * len(lines_h)
shortnames_h = [[]] * len(lines_h)
names_h = [[]] * len(lines_h)

gh, _, _, _ = fw.extract_3d_var(grib, 'gh')
oro, _, _, _ = fw.extract_3d_var(grib, 'orog')

part_gh = gh[:,idx_y, idx_x]

Z = part_gh.T  # (npart, nlev)

# index of first level >= parcel altitude
k_hi = np.array([np.searchsorted(Z[i], part_z_[i], side="left") for i in range(Z.shape[0])])

# convert to "below" index
k_lo = k_hi - 1

# clamp to valid range
k_hi = np.clip(k_hi, 0, Z.shape[1]-1)
k_lo = np.clip(k_lo, 0, Z.shape[1]-1)

asdfasdf

for hidx, h in enumerate(lines_h[:]):
    shortName = h.split(',')[1].strip()
    msg = grib.select(shortName = shortName, typeOfLevel = 'hybrid')[0]
    data, _, _ = msg.data()
    units_h[hidx] = msg.units
    names_h[hidx] = msg.name
    shortnames_s[hidx] = shortName
    part_surface[hidx, :, didx] = (data[idx_y, idx_x] * msg.scaleValuesBy) + msg.offsetValuesBy
'''
for sidx, s in enumerate(lines_s[:]):
    shortName = s.split(',')[1].strip()
    msg = grib.select(shortName = shortName, typeOfLevel = 'surface')[0]
    data, _, _ = msg.data()
    units_s[sidx] = msg.units
    names_s[sidx] = msg.name
    shortnames_s[sidx] = shortName
    part_surface[sidx, :, didx] = (data[idx_y, idx_x] * msg.scaleValuesBy) + msg.offsetValuesBy

for aidx, a in enumerate(lines_a[:]):
    shortName = a.split(',')[1].strip()
    msg = grib.select(shortName = shortName, typeOfLevel = 'aboveGroundLevel')[0]
    data, _, _ = msg.data()
    units_a[aidx] = msg.units
    names_a[aidx] = msg.name
    shortnames_a[aidx] = shortName
    part_agl[aidx, :, didx] = (data[idx_y, idx_x] * msg.scaleValuesBy) + msg.offsetValuesBy

######### create output file
basename = os.path.basename(folder)
dirname = os.path.dirname(folder)
sname = os.path.join(dirname, 'CLMS_flexwrf_HRRR.%s.nc' % basename)
## open file and create dimensions
ds_out = nc.Dataset(sname, 'w')
time_dim = ds_out.createDimension(dimname='time', size=np.shape(part_lon)[1])
particle_dim = ds_out.createDimension(dimname='particle', size=np.shape(part_lon)[0])

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
time_var[:] = matlab_times
part_mass_var[:] = part_mass
part_alt_var[:] = part_alt
hmix_var[:] = part_hmix

##### next step through all surface variables: create variable, write attributes and load data
for sidx, s in enumerate(lines_s[:]):
    tmp_var = ds_out.createVariable(varname= shortnames_s[sidx],datatype='f4',dimensions=('particle','time'))
    tmp_var.units = units_s[sidx]
    tmp_var.setncattr_string('name',names_s[sidx])
    tmp_var[:] = part_surface[sidx]

##### next step through all AGL variables: create variable, write attributes and load data
for aidx, a in enumerate(lines_a[:]):
    tmp_var = ds_out.createVariable(varname= shortnames_a[aidx],datatype='f4',dimensions=('particle','time'))
    tmp_var.units = units_a[aidx]
    tmp_var.setncattr_string('name',names_a[aidx])
    tmp_var[:] = part_agl[aidx]

##### next step through all 3d variables: create variable, write attributes and load data
for hidx, a in enumerate(lines_h[:]):
    tmp_var = ds_out.createVariable(varname= shortnames_h[hidx],datatype='f4',dimensions=('particle','time'))
    tmp_var.units = units_h[hidx]
    tmp_var.setncattr_string('name',names_h[hidx])
    tmp_var[:] = part_hybrid[hidx]

''
### loop over matlab times
#for iidx, d in enumerate(dates):
#    date = dates[iidx]
#    year = str(date.year)
#    month = str(date.month).zfill(2)
#    day = str(date.day).zfill(2)
    # for
#    part_lon_ = part_lon[:, iidx]
#    part_lat_ = part_lat[:, iidx]
#    part_hmix_ = part_hmix[:, iidx]

    ######## following block does the following for each variable, for a given time:
    ### 1. finds the file closest in time to the timepoint
    ### 2. opens file and reads in variables of interest
    ### 3. finds indices that correspond from particle lat/lon to VOI
    ### 4. concatenate variables

    lon = lon[idx1]
    lat = lat[idx2]
    BA = ds.variables['burned_fraction'][0, idx2, idx1] * ds.variables['burned_fraction'].scale_factor
    ### get indices for each particle that corresponds to lon/lat variables
    lon_idx = np.digitize(part_lon_, lon)
    lat_idx = np.digitize(part_lat_, lat)
    part_ba[:, iidx] = BA[lat_idx, lon_idx]
    ds.close()

### create output file and save it
basename = os.path.basename(folder)
dirname = os.path.dirname(folder)
sname = os.path.join(dirname, 'CLMS_flexwrf_HRRR.%s.nc' % basename)

## open file and create dimensions
ds_out = nc.Dataset(sname, 'w')
time_dim = ds_out.createDimension(dimname='time', size=np.shape(part_lon)[1])
particle_dim = ds_out.createDimension(dimname='particle', size=np.shape(part_lon)[0])

## create variables
time_var = ds_out.createVariable(varname='time', datatype='f8', dimensions='time')
part_mass_var = ds_out.createVariable(varname='part_mass', datatype='f4', dimensions=('particle', 'time'))
part_alt_var = ds_out.createVariable(varname='part_alt', datatype='f4', dimensions=('particle', 'time'))
hmix_var = ds_out.createVariable(varname='part_hmix', datatype='f4', dimensions=('particle', 'time'))
ba_var = ds_out.createVariable(varname='burnt_area', datatype='f4', dimensions=('particle', 'time'))
npp_var = ds_out.createVariable(varname='npp', datatype='f4', dimensions=('particle', 'time'))
ndvi_var = ds_out.createVariable(varname='ndvi', datatype='f4', dimensions=('particle', 'time'))
fapar_var = ds_out.createVariable(varname='fapar', datatype='f4', dimensions=('particle', 'time'))
lai_var = ds_out.createVariable(varname='lai', datatype='f4', dimensions=('particle', 'time'))
fcover_var = ds_out.createVariable(varname='fcover', datatype='f4', dimensions=('particle', 'time'))
swi_001_var = ds_out.createVariable(varname='swi_001', datatype='f4', dimensions=('particle', 'time'))
swi_005_var = ds_out.createVariable(varname='swi_005', datatype='f4', dimensions=('particle', 'time'))
swi_010_var = ds_out.createVariable(varname='swi_010', datatype='f4', dimensions=('particle', 'time'))

# write attributes
time_var.units = 'matplotlib date number'
time_var.setncattr_string('name', 'datenumber')
part_mass_var.units = 'kg'
part_mass_var.setncattr_string('name', 'Particle mass')
part_alt_var.units = 'm'
part_alt_var.setncattr_string('name', 'Particle altitude AGL')
hmix_var.units = 'm'
hmix_var.setncattr_string('name', 'Mixing layer height AGL')
ba_var.units = 'fraction'
ba_var.setncattr_string('name', 'Burnt fraction')
ba_var.setncattr_string('long_name', 'Fraction of pixel surface affected by fire at the day of the burn detection')
npp_var.units = 'gC / m2/ day'
npp_var.setncattr_string('name', 'Net primary production')
npp_var.setncattr_string('long_name', 'Net Primary Production 333 m')
ndvi_var.units = ''
ndvi_var.setncattr_string('name', 'normalized_difference_vegetation_index')
ndvi_var.setncattr_string('long_name', 'Normalized Difference Vegetation Index 333m')
fapar_var.units = ''
fapar_var.setncattr_string('name',
                           'fraction_of_surface_downwelling_photosynthetic_radiative_flux_absorbed_by_vegetation')
fapar_var.setncattr_string('long_name', 'Fraction of Absorbed Photosynthetically Active Radiation 333m')
lai_var.units = 'm2 / m2'
lai_var.setncattr_string('name', 'leaf_area_index')
lai_var.setncattr_string('long_name', 'Leaf Area Index 333m')
fcover_var.units = ''
fcover_var.setncattr_string('name', 'vegetation_area_fraction')
fcover_var.setncattr_string('long_name', 'Fraction of green Vegetation Cover 333m')
swi_001_var.units = '%'
swi_001_var.setncattr_string('name', 'soil_water_index_001')
swi_001_var.setncattr_string('long_name', 'Soil Water Index with T=1')
swi_005_var.units = '%'
swi_005_var.setncattr_string('name', 'soil_water_index_005')
swi_005_var.setncattr_string('long_name', 'Soil Water Index with T=5')
swi_010_var.units = '%'
swi_010_var.setncattr_string('name', 'soil_water_index_010')
swi_010_var.setncattr_string('long_name', 'Soil Water Index with T=10')

# write data
hmix_var[:] = part_hmix
time_var[:] = matlab_times
ba_var[:] = part_ba
npp_var[:] = part_npp
ndvi_var[:] = part_ndvi
fapar_var[:] = part_fapar
lai_var[:] = part_lai
fcover_var[:] = part_fcover
swi_001_var[:] = part_swi_01
swi_005_var[:] = part_swi_05
swi_010_var[:] = part_swi_10
# close file
ds_out.close()


# if __name__ == "__main__":
#     folder = sys.argv[1]
#     print(folder)
#     trajectory_processing(folder)
'''
