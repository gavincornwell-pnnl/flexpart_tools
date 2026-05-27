import numpy as np
import netCDF4 as nc
import flexpart as fp
import os
import sys
import matplotlib.dates as dt
from scipy.stats import binned_statistic_2d as bs2d
import matplotlib.pyplot as plt

def trajectory_processing(foldername):
    ### get all BURNT AREA files and then find the date that they correspond to
    ba_folder = os.path.join('/rcfs/scratch/corn062/CLMS/bio-geophysical/burnt_area/ba_global_300m_monthly_v4')
    #ba_fnames = os.listdir(ba_folder)
    fnames = np.zeros(0,dtype=str)
    for root, dirs, files in os.walk(ba_folder):
        for name in files:
            fnames = np.append(fnames,os.path.join(root,name))
    ba_fnames =np.sort([i for i in fnames if 'BA' in i])
    ba_dates = np.zeros(len(ba_fnames),dtype=dt.datetime.datetime)
    ## times for each file
    for tidx,tmp in enumerate(ba_fnames):
        tmp2 = tmp.split('_')
        tmp2 = [i for i in tmp2 if i.startswith('20')]
        if len(tmp2) > 0:
            tmp2 = tmp2[0]
            ba_dates[tidx] = dt.num2date(dt.datestr2num(tmp2))

    ### get all NET PRIMARY PRODUCTION files and then find the date that they correspond to
    npp_folder = os.path.join('/rcfs/scratch/corn062/CLMS/bio-geophysical/net-gross_primary_production/npp_global_300m_10daily_v2')
    #npp_fnames = os.listdir(npp_folder)
    fnames = np.zeros(0,dtype=str)
    for root, dirs, files in os.walk(npp_folder):
        for name in files:
            fnames = np.append(fnames,os.path.join(root,name))
    npp_fnames =np.sort([i for i in fnames if 'V2.0.2' in i])
    npp_dates = np.zeros(len(npp_fnames),dtype=dt.datetime.datetime)
    ## times for each file
    for tidx,tmp in enumerate(npp_fnames):
        tmp2 = tmp.split('_')
        tmp2 = [i for i in tmp2 if i.startswith('20')]
        if len(tmp2) > 0:
            tmp2 = tmp2[0]
            npp_dates[tidx] = dt.num2date(dt.datestr2num(tmp2))

    ### get all NDVI files and then find the date that they correspond to
    ndvi_folder = os.path.join('/rcfs/scratch/corn062/CLMS/bio-geophysical/vegetation_indices/ndvi_global_300m_10daily_v3')
    fnames = np.zeros(0,dtype=str)
    for root, dirs, files in os.walk(ndvi_folder):
        for name in files:
            fnames = np.append(fnames,os.path.join(root,name))
    ndvi_fnames =np.sort([i for i in fnames if 'NDVI' in i])
    ndvi_dates = np.zeros(len(ndvi_fnames),dtype=dt.datetime.datetime)
    ## times for each file
    for tidx,tmp in enumerate(ndvi_fnames):
        tmp2 = tmp.split('_')
        tmp2 = [i for i in tmp2 if i.startswith('20')]
        if len(tmp2) > 0:
            tmp2 = tmp2[0]
            ndvi_dates[tidx] = dt.num2date(dt.datestr2num(tmp2))

    ### get all FAPAR files and then find the date that they correspond to
    fapar_folder = os.path.join('/rcfs/scratch/corn062/CLMS/bio-geophysical/vegetation_properties/fapar_global_300m_10daily_v2')
    fnames = np.zeros(0,dtype=str)
    for root, dirs, files in os.walk(fapar_folder):
        for name in files:
            fnames = np.append(fnames,os.path.join(root,name))
    fapar_fnames =np.sort([i for i in fnames if 'RT6_' in i])
    fapar_dates = np.zeros(len(fapar_fnames),dtype=dt.datetime.datetime)
    ## times for each file
    for tidx,tmp in enumerate(fapar_fnames):
        tmp2 = tmp.split('_')
        tmp2 = [i for i in tmp2 if i.startswith('20')]
        if len(tmp2) > 0:
            tmp2 = tmp2[0]
            fapar_dates[tidx] = dt.num2date(dt.datestr2num(tmp2))

    ### get all LEAF AREA INDEX files and then find the date that they correspond to
    lai_folder = os.path.join('/rcfs/scratch/corn062/CLMS/bio-geophysical/vegetation_properties/lai_global_300m_10daily_v2')
    fnames = np.zeros(0,dtype=str)
    for root, dirs, files in os.walk(lai_folder):
        for name in files:
            fnames = np.append(fnames,os.path.join(root,name))
    lai_fnames =np.sort([i for i in fnames if 'RT6_' in i])
    lai_dates = np.zeros(len(lai_fnames),dtype=dt.datetime.datetime)
    ## times for each file
    for tidx,tmp in enumerate(lai_fnames):
        tmp2 = tmp.split('_')
        tmp2 = [i for i in tmp2 if i.startswith('20')]
        if len(tmp2) > 0:
            tmp2 = tmp2[0]
            lai_dates[tidx] = dt.num2date(dt.datestr2num(tmp2))

    ### get all FCOVER files and then find the date that they correspond to
    fc_folder = os.path.join('/rcfs/scratch/corn062/CLMS/bio-geophysical/vegetation_properties/fcover_global_300m_10daily_v2')
    fnames = np.zeros(0,dtype=str)
    for root, dirs, files in os.walk(fc_folder):
        for name in files:
            fnames = np.append(fnames,os.path.join(root,name))
    fcover_fnames =np.sort([i for i in fnames if 'RT6_' in i])
    fcover_dates = np.zeros(len(fcover_fnames),dtype=dt.datetime.datetime)
    ## times for each file
    for tidx,tmp in enumerate(fcover_fnames):
        tmp2 = tmp.split('_')
        tmp2 = [i for i in tmp2 if i.startswith('20')]
        if len(tmp2) > 0:
            tmp2 = tmp2[0]
            fcover_dates[tidx] = dt.num2date(dt.datestr2num(tmp2))

    ### get all SOIL WATER INDEX files and then find the date that they correspond to
    swi_folder = os.path.join('/rcfs/scratch/corn062/CLMS/bio-geophysical/soil_water_index/swi_global_12.5km_10daily_v4')
    fnames = np.zeros(0,dtype=str)
    for root, dirs, files in os.walk(swi_folder):
        for name in files:
            fnames = np.append(fnames,os.path.join(root,name))
    swi_fnames =np.sort([i for i in fnames if 'SWI' in i])
    swi_dates = np.zeros(len(swi_fnames),dtype=dt.datetime.datetime)
    ## times for each file
    for tidx,tmp in enumerate(swi_fnames):
        tmp2 = tmp.split('_')
        tmp2 = [i for i in tmp2 if i.startswith('20')]
        if len(tmp2) > 0:
            tmp2 = tmp2[0]
            swi_dates[tidx] = dt.num2date(dt.datestr2num(tmp2))

    #### get filelist from the folder to be analyzed
    # folder = '/rcfs/scratch/corn062/AGINSGP/ncep/20220410_010000_050'
    fnamelist = os.listdir(folder)
    fnamelist = [i for i in fnamelist if i.endswith('.nc')]
    fnamelist = np.sort([i for i in fnamelist if i.startswith('partoutput')])[::-1]
    fname = os.path.join(folder,fnamelist[0])
    part_lon, _, _, _, _, _, _ = fp.read_partposition_fp11(fname)

    ### Create variables that will be written
    part_ba = np.ones(np.shape(part_lon)) * -9999
    part_npp = np.ones(np.shape(part_lon)) * -9999
    part_ndvi = np.ones(np.shape(part_lon)) * -9999
    part_fapar = np.ones(np.shape(part_lon)) * -9999
    part_lai = np.ones(np.shape(part_lon)) * -9999
    part_fcover = np.ones(np.shape(part_lon)) * -9999
    part_swi_01 = np.ones(np.shape(part_lon)) * -9999
    part_swi_05 = np.ones(np.shape(part_lon)) * -9999
    part_swi_10 = np.ones(np.shape(part_lon)) * -9999

    ## get individual particle output file and read in its data
    fidx = 0
    fname = os.path.join(folder,fnamelist[fidx])
    part_lon, part_lat, part_mass, part_z, part_hmix, time, matlab_times = fp.read_partposition_fp11(fname)
    days = np.round(matlab_times / (24 * 3600))
    dates = dt.num2date(days)

    ### loop over matlab times
    for iidx, d in enumerate(dates):
        date = dates[iidx]
        year = str(date.year)
        month = str(date.month).zfill(2)
        day = str(date.day).zfill(2)
        # for
        part_lon_ = part_lon[:,iidx]
        part_lat_ = part_lat[:,iidx]
        part_hmix_ = part_hmix[:,iidx]

        ######## following block does the following for each variable, for a given time:
        ### 1. finds the file closest in time to the timepoint
        ### 2. opens file and reads in variables of interest
        ### 3. finds indices that correspond from particle lat/lon to VOI
        ### 4. concatenate variables

        #### BURNT AREA
        ## find the file that is closest to the current date
        idx = np.argmin(np.abs(date - ba_dates))
        ## open file and get variables within a certain box
        print(ba_fnames[idx])
        ds = nc.Dataset(ba_fnames[idx])
        lon = ds.variables['lon'][:]
        lat = ds.variables['lat'][:]
        idx1 = np.where((lon >= -135) & (lon < -60))[0]
        idx2 = np.where((lat >=20) & (lat <55))[0]
        lon = lon[idx1]
        lat = lat[idx2]
        BA = ds.variables['burned_fraction'][0,idx2,idx1] * ds.variables['burned_fraction'].scale_factor
        ### get indices for each particle that corresponds to lon/lat variables
        lon_idx = np.digitize(part_lon_,lon)
        lat_idx = np.digitize(part_lat_,lat)
        part_ba[:,iidx] = BA[lat_idx,lon_idx]
        ds.close()

        #### NPP
        ## find the file that is closest to the current date
        idx = np.argmin(np.abs(date - npp_dates))
        ## open file and get variables within a certain box
        print(npp_fnames[idx])
        ds = nc.Dataset(npp_fnames[idx])
        lon = ds.variables['lon'][:]
        lat = ds.variables['lat'][:]
        idx1 = np.where((lon >= -135) & (lon < -60))[0]
        idx2 = np.where((lat >=20) & (lat <55))[0]
        lon = lon[idx1]
        lat = lat[idx2]
        NPP = ds.variables['NPP'][0,idx2,idx1] * ds.variables['NPP'].scale_factor
        ### get indices for each particle that corresponds to lon/lat variables
        lon_idx = np.digitize(part_lon_,lon)
        lat_idx = np.digitize(part_lat_,lat)
        part_npp[:,iidx] = NPP[lat_idx,lon_idx]
        ds.close()

        #### NDVI
        ## find the file that is closest to the current date
        idx = np.argmin(np.abs(date - ndvi_dates))
        ## open file and get variables within a certain box
        print(ndvi_fnames[idx])
        ds = nc.Dataset(ndvi_fnames[idx])
        lon = ds.variables['lon'][:]
        lat = ds.variables['lat'][:]
        idx1 = np.where((lon >= -135) & (lon < -60))[0]
        idx2 = np.where((lat >=20) & (lat <55))[0]
        lon = lon[idx1]
        lat = lat[idx2]
        NDVI = ds.variables['NDVI'][0,idx2,idx1] * ds.variables['NDVI'].scale_factor + ds.variables['NDVI'].add_offset
        ### get indices for each particle that corresponds to lon/lat variables
        lon_idx = np.digitize(part_lon_,lon)
        lat_idx = np.digitize(part_lat_,lat)
        part_ndvi[:,iidx] = NDVI[lat_idx,lon_idx]
        ds.close()

        #### LAI
        ## find the file that is closest to the current date
        idx = np.argmin(np.abs(date - lai_dates))
        ## open file and get variables within a certain box
        print(lai_fnames[idx])
        ds = nc.Dataset(lai_fnames[idx])
        lon = ds.variables['lon'][:]
        lat = ds.variables['lat'][:]
        idx1 = np.where((lon >= -135) & (lon < -60))[0]
        idx2 = np.where((lat >=20) & (lat <55))[0]
        lon = lon[idx1]
        lat = lat[idx2]
        LAI = ds.variables['LAI'][0,idx2,idx1] * ds.variables['LAI'].scale_factor + ds.variables['LAI'].add_offset
        ### get indices for each particle that corresponds to lon/lat variables
        lon_idx = np.digitize(part_lon_,lon)
        lat_idx = np.digitize(part_lat_,lat)
        part_lai[:,iidx] = LAI[lat_idx,lon_idx]
        ds.close()

        #### FCOVER
        ## find the file that is closest to the current date
        idx = np.argmin(np.abs(date - fcover_dates))
        ## open file and get variables within a certain box
        print(fcover_fnames[idx])
        ds = nc.Dataset(fcover_fnames[idx])
        lon = ds.variables['lon'][:]
        lat = ds.variables['lat'][:]
        idx1 = np.where((lon >= -135) & (lon < -60))[0]
        idx2 = np.where((lat >=20) & (lat <55))[0]
        lon = lon[idx1]
        lat = lat[idx2]
        FCOVER = ds.variables['FCOVER'][0,idx2,idx1] * ds.variables['FCOVER'].scale_factor + ds.variables['FCOVER'].add_offset
        ### get indices for each particle that corresponds to lon/lat variables
        lon_idx = np.digitize(part_lon_,lon)
        lat_idx = np.digitize(part_lat_,lat)
        part_fcover[:,iidx] = FCOVER[lat_idx,lon_idx]
        ds.close()

        #### SOIL WATER INDEX
        ## find the file that is closest to the current date
        idx = np.argmin(np.abs(date - swi_dates))
        ## open file and get variables within a certain box
        print(swi_fnames[idx])
        ds = nc.Dataset(swi_fnames[idx])
        lon = ds.variables['lon'][:]
        lat = ds.variables['lat'][:]
        idx1 = np.where((lon >= -135) & (lon < -60))[0]
        idx2 = np.where((lat >=20) & (lat <55))[0]
        lon = lon[idx1]
        lat = lat[idx2]
        SWI_001 = ds.variables['SWI_001'][0,idx2,idx1] * ds.variables['SWI_001'].scale_factor
        SWI_005 = ds.variables['SWI_005'][0,idx2,idx1] * ds.variables['SWI_005'].scale_factor
        SWI_010 = ds.variables['SWI_010'][0,idx2,idx1] * ds.variables['SWI_010'].scale_factor
        ### get indices for each particle that corresponds to lon/lat variables
        lon_idx = np.digitize(part_lon_,lon)
        lat_idx = np.digitize(part_lat_,lat)
        part_swi_01[:,iidx] = SWI_001[lat_idx,lon_idx]
        part_swi_05[:,iidx] = SWI_005[lat_idx,lon_idx]
        part_swi_10[:,iidx] = SWI_010[lat_idx,lon_idx]
        ds.close()

    ### create output file and save it
    sname = os.path.join(folder,'CLMS_processed_%s.nc' % tmp.split('_')[1])
    print(sname)

    ## open file and create dimensions
    ds_out=nc.Dataset(sname,'w')
    time_dim=ds_out.createDimension(dimname='time',size=np.shape(part_lon)[1])
    particle_dim=ds_out.createDimension(dimname='particle',size=np.shape(part_lon)[0])

    ## create variables
    time_var=ds_out.createVariable(varname='time',datatype='f8',dimensions='time')
    part_mass_var = ds_out.createVariable(varname='part_mass', datatype='f4', dimensions=('particle', 'time'))
    part_alt_var=ds_out.createVariable(varname='part_alt',datatype='f4',dimensions=('particle','time'))
    hmix_var=ds_out.createVariable(varname='part_hmix',datatype='f4',dimensions=('particle','time'))
    ba_var=ds_out.createVariable(varname='burnt_area',datatype='f4',dimensions=('particle','time'))
    npp_var=ds_out.createVariable(varname='npp',datatype='f4',dimensions=('particle','time'))
    ndvi_var=ds_out.createVariable(varname='ndvi',datatype='f4',dimensions=('particle','time'))
    fapar_var=ds_out.createVariable(varname='fapar',datatype='f4',dimensions=('particle','time'))
    lai_var=ds_out.createVariable(varname='lai',datatype='f4',dimensions=('particle','time'))
    fcover_var=ds_out.createVariable(varname='fcover',datatype='f4',dimensions=('particle','time'))
    swi_001_var=ds_out.createVariable(varname='swi_001',datatype='f4',dimensions=('particle','time'))
    swi_005_var=ds_out.createVariable(varname='swi_005',datatype='f4',dimensions=('particle','time'))
    swi_010_var=ds_out.createVariable(varname='swi_010',datatype='f4',dimensions=('particle','time'))

    # write attributes
    time_var.units='matplotlib date number'
    time_var.setncattr_string('name','datenumber')
    part_mass_var.units = 'kg'
    part_mass_var.setncattr_string('name', 'Particle mass')
    part_alt_var.units='m'
    part_alt_var.setncattr_string('name','Particle altitude AGL')
    hmix_var.units='m'
    hmix_var.setncattr_string('name','Mixing layer height AGL')
    ba_var.units='fraction'
    ba_var.setncattr_string('name','Burnt fraction')
    ba_var.setncattr_string('long_name','Fraction of pixel surface affected by fire at the day of the burn detection')
    npp_var.units='gC / m2/ day'
    npp_var.setncattr_string('name','Net primary production')
    npp_var.setncattr_string('long_name','Net Primary Production 333 m')
    ndvi_var.units=''
    ndvi_var.setncattr_string('name','normalized_difference_vegetation_index')
    ndvi_var.setncattr_string('long_name','Normalized Difference Vegetation Index 333m')
    fapar_var.units=''
    fapar_var.setncattr_string('name','fraction_of_surface_downwelling_photosynthetic_radiative_flux_absorbed_by_vegetation')
    fapar_var.setncattr_string('long_name','Fraction of Absorbed Photosynthetically Active Radiation 333m')
    lai_var.units='m2 / m2'
    lai_var.setncattr_string('name','leaf_area_index')
    lai_var.setncattr_string('long_name','Leaf Area Index 333m')
    fcover_var.units=''
    fcover_var.setncattr_string('name','vegetation_area_fraction')
    fcover_var.setncattr_string('long_name','Fraction of green Vegetation Cover 333m')
    swi_001_var.units='%'
    swi_001_var.setncattr_string('name','soil_water_index_001')
    swi_001_var.setncattr_string('long_name','Soil Water Index with T=1')
    swi_005_var.units='%'
    swi_005_var.setncattr_string('name','soil_water_index_005')
    swi_005_var.setncattr_string('long_name','Soil Water Index with T=5')
    swi_010_var.units='%'
    swi_010_var.setncattr_string('name','soil_water_index_010')
    swi_010_var.setncattr_string('long_name','Soil Water Index with T=10')

    # write data
    hmix_var[:]=part_hmix
    time_var[:]=matlab_times
    ba_var[:]=part_ba
    npp_var[:]=part_npp
    ndvi_var[:]=part_ndvi
    fapar_var[:]=part_fapar
    lai_var[:]=part_lai
    fcover_var[:]=part_fcover
    swi_001_var[:]=part_swi_01
    swi_005_var[:]=part_swi_05
    swi_010_var[:]=part_swi_10
    # close file
    ds_out.close()


if __name__ == "__main__":
    folder = sys.argv[1]
    print(folder)
    trajectory_processing(folder)
                                              

