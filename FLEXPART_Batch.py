def writeCOMMAND(startDate, endDate, folder):
    # Function to write COMMAND file for FLEXPART backtrajectories.
    # Start date is a 1x1 date number marking the EARLIEST time point
    # for the simulation (ie the end of the trajectory).
    # End date is a 1x1  date number marking the LATEST time point for
    # the simulation (ie the beginning of the trajectory).
    # folder is the directory where the COMMAND file should be written to
    import matplotlib.dates as dt
    import os
    format = '%Y%m%d %H%M%S'
    if type(startDate) != dt.datetime.datetime:
        startDate = dt.num2date(startDate)
    startDateStr = dt.datetime.datetime.strftime(startDate, format=format)
    if type(endDate) != dt.datetime.datetime:
        endDate = dt.num2date(endDate)
    endDateStr = dt.datetime.datetime.strftime(endDate, format=format)
    fid = open(os.path.join(folder, 'COMMAND'), 'w')
    # writing header infmtion
    fid.write('********************************************************************************\n');
    fid.write('*                                                                              *\n');
    fid.write('*      Input file for the Lagrangian particle dispersion model FLEXPART        *\n');
    fid.write('*                           Please select your options                         *\n');
    fid.write('*                                                                              *\n');
    fid.write('********************************************************************************\n');
    fid.write('&COMMAND\n');
    # writing command options
    fid.write(' LDIRECT=              -1, ! Simulation direction in time   ; 1(forward) or -1 (backward)\n')
    fid.write(' IBDATE=         ' + startDateStr[
                                    0:8] + ', ! Start date of the simulation   ; YYYYMMDD: YYYY=year, MM=month, DD=day\n')
    fid.write(' IBTIME=          ' + startDateStr[
                                     8:] + ', ! Start time of the simulation   ; HHMISS: HH=hours, MI=min, SS=sec; UTC \n')
    fid.write(' IEDATE=         ' + endDateStr[
                                    0:8] + ', ! Start date of the simulation   ; YYYYMMDD: YYYY=year, MM=month, DD=day\n')
    fid.write(' IETIME=          ' + endDateStr[
                                     8:] + ', ! Start time of the simulation   ; HHMISS: HH=hours, MI=min, SS=sec; UTC \n')
    fid.write(
        ' LOUTSTEP=           3600, ! Interval of model output; average concentrations calculated every LOUTSTEP (s)  \n')
    fid.write(' LOUTAVER=           3600, ! Interval of output averaging (s)\n')
    fid.write(
        ' LOUTSAMPLE=          900, ! Interval of output sampling  (s), higher stat. accuracy with shorter intervals\n')

    fid.write(' LOUTRESTART=       28800, ! Interval of writing restart files (s), switched off when set to -1\n')
    fid.write(' LRECOUTSTEP=        3600, ! Interval of model output at receptors (s)\n')
    fid.write(' LRECOUTAVER=        3600, ! Interval of receptor output averaging (s)\n')
    fid.write(' LRECOUTSAMPLE=      1200, ! Interval of receptor output sampling (s)\n')
    # fid.write(' ITSPLIT=        99999999, ! Interval of particle splitting (s)\n')
    fid.write(' LSYNCTIME=           900, ! All processes are synchronized to this time interval (s)\n')
    fid.write(
        ' CTL=               5.000, ! CTL>1, ABL time step = (Lagrangian timescale (TL))/CTL, uses LSYNCTIME if CTL<0\n')
    fid.write(' IFINE=                 10, ! Reduction for time step in vertical transport, used only if CTL>1\n');
    fid.write(
        ' IOUT=                  1, ! Output type: [1]mass 2]pptv 3]1&2 4]plume 5]1&4, +8 for NetCDF output     \n');
    fid.write(' IPOUT=                 1, ! Particle position output: 0]no 1]every output 2]only at end   \n');
    fid.write(
        ' LSUBGRID=              0, ! Increase of ABL heights due to sub-grid scale orographic variations;[0]off 1]on \n');
    fid.write(' LCONVECTION=          	1, ! Switch for convection parameterization;0]off [1]on    \n');
    fid.write(' LTURBULENCE=           1, ! Switch for turbulence parameterisation;0]off [1]on\n')
    fid.write(
        ' LTURBULENCE_MESO=      0, ! Switch for mesoscale turbulence parameterisation;0]off (recommended) [1]on\n')
    fid.write(
        ' LAGESPECTRA=          	0, ! Switch for calculation of age spectra (needs AGECLASSES);[0]off 1]on  \n');
    fid.write(
        ' IPIN=                  0, ! Warm start from particle dump (needs previous partposit_end file); [0]no 1]yes  \n');
    fid.write(
        ' IOUTPUTFOREACHRELEASE= 1, ! Separate output fields for each location in the RELEASE file; [0]no 1]yes \n');
    fid.write(' IFLUX=                 0, ! Output of mass fluxes through output grid box boundaries\n');
    fid.write(
        ' MDOMAINFILL=          	0, ! Switch for domain-filling, if limited-area particles generated at boundary\n');
    fid.write(' IND_SOURCE=            1, ! Unit to be used at the source   ;  [1]mass 2]mass mixing ratio \n');
    fid.write(
        ' IND_RECEPTOR=          1, ! Unit to be used at the receptor; [1]mass 2]mass mixing ratio 3]wet depo. 4]dry depo.\n');
    fid.write(' MQUASILAG=             0, ! Quasi-Lagrangian mode to track individual numbered particles \n');
    fid.write(' NESTED_OUTPUT=         0, ! Output also for a nested domain \n');
    fid.write(' LNETCDFOUT=            1, ! Gridded netcdf output: [0]no [1]yes\n')
    fid.write(
        ' LINIT_COND=            0, ! Output sensitivity to initial conditions (bkw mode only) [0]off 1]conc 2]mmr \n');
    fid.write(' SFC_ONLY=              0, ! Output only for the lowest model layer, used w/ LINIT_COND=1 or 2\n');
    fid.write(
        ' CBLFLAG=               1, ! Skewed, not Gaussian turbulence in the convective ABL, need large CTL and IFINE\n');
    fid.write(' OHFIELDS_PATH= "/people/corn062/flexpart11/options/oh_fields", ! Default path for OH file\n');

    fid.write(
        ' NXSHIFT=               0, ! Shift of the global meteorological data. Default 359 for ECMWF and 0 for GFS if not given\n')
    fid.write(
        ' MAXTHREADGRID=         4, ! Set maximum number of threads for doing grid computations. Recommended to set this no higher than 16. High numbers create more overhead and a larger memory footprint, 1=no parallelisation on grid.\n')
    fid.write(
        ' MAXFILESIZE=       10000, ! Maximum output of each partoutput NetCDF-4 file in Mb before a new one is created\n')
    fid.write(' LOGVERTINTERP=         0, ! Flag to set all vertical interpolation to logarithmic instead of linear\n')
    fid.write(' LCMOUTPUT=             0, ! Switch for the Linear Chemistry Module; [0] off [1] on\n')
    fid.write('/\n');
    fid.close();
    return


def writePathname(releaseFolder, outputFolder, metFolder, folder):
    import os
    if releaseFolder[-1] != '/':
        releaseFolder = releaseFolder + '/'
    if outputFolder[-1] != '/':
        outputFolder = outputFolder + '/'
    if metFolder[-1] != '/':
        metFolder = metFolder + '/'
    fid = open(os.path.join(folder, 'pathnames'), 'w')
    fid.write(releaseFolder + '\n')  # folder with the COMMAND, RELEASE and other run associated files
    fid.write(outputFolder + '\n')  # folder with the COMMAND, RELEASE and other run associated files
    fid.write(metFolder + '\n');  # folder with metData
    fid.write(metFolder + 'AVAILABLE\n')  # location of available file
    fid.write('============================================')
    fid.close()
    return


def writeRELEASES(gridCorner1, gridCorner2, altitude, releaseDateStart, releaseDateEnd, trajName, numParticles,
                  writeFolder, releaseID):
    # Script to automate the writing of the RELEASES file for FLEXPART
    # backtrajectory calculations
    # gridCorner1 is a Nx2 matrix populated with pairs of the lower left corner of the
    # desired release box
    # gridCorner2 is a Nx2 matrix populated with pairs of the upper right corner of the
    # desired release box
    # Altitude is a Nx1 or Nx2 matrix populated with the release height or the upper and lower bounds of the
    # desired release box--should be in meters above sea level
    # Release date a Nx1 matrix of matlab date numbers corresponding to the
    # initialization of the backtrajectory
    # trajName is a Nx1 cell array corresponding to the desired trajectory
    # names
    # Each release is currently set to be one hour and 500 particles in a
    # domain filled box. These can be changed in the first few lines of code if
    # desired.
    # print header
    import os
    import matplotlib.dates as dt
    import numpy as np

    fid = open(os.path.join(writeFolder, 'RELEASES'), 'w')
    fid.write('*************************************************************************\n')
    fid.write('*                                                                       *\n')
    fid.write('*                                                                       *\n')
    fid.write('*                                                                       *\n')
    fid.write('*   Input file for the Lagrangian particle dispersion model FLEXPART    *\n')
    fid.write('*                        Please select your options                     *\n')
    fid.write('*                                                                       *\n')
    fid.write('*                                                                       *\n')
    fid.write('*                                                                       *\n')
    fid.write('*************************************************************************\n')
    fid.write('&RELEASES_CTRL\n')
    fid.write(' NSPEC      =          %d, ! Total number of species\n' % int(len(releaseID)))
    fid.write(' SPECNUM_REL=          %d, ! Species numbers in directory SPECIES\n' % np.array(releaseID, dtype=int))
    fid.write('/\n')
    fid.write('&RELEASE\n')
    format1 = '%Y%m%d %H%M%S'
    format2 = '%Y%m%d_%H%M%S'
    for iidx, i in enumerate(gridCorner1):
        tmpRDS = releaseDateStart[iidx]
        if type(tmpRDS) != dt.datetime.datetime:
            tmpRDS = dt.num2date(tmpRDS)
        tmpRDE = releaseDateEnd[iidx]
        if type(tmpRDE) != dt.datetime.datetime:
            tmpRDE = dt.num2date(tmpRDE)
        startDateStr = dt.datetime.datetime.strftime(tmpRDS, format=format1)
        endDateStr = dt.datetime.datetime.strftime(tmpRDE, format=format1)
        trajNameFinal = trajName + '_' + dt.datetime.datetime.strftime(tmpRDE, format=format2) + '_' + str(
            iidx + 1).zfill(3)
        fid.write(' IDATE1  =       ' + startDateStr[
                                        0:8] + ', ! Release start date, YYYYMMDD: YYYY=year, MM=month, DD=day\n')  # start time of release
        fid.write(' ITIME1  =        ' + startDateStr[
                                         8:] + ', ! Release start time in UTC HHMISS: HH hours, MI=minutes, SS=seconds\n')  # end time of release
        fid.write(' IDATE2  =       ' + endDateStr[
                                        0:8] + ', ! Release start date, YYYYMMDD: YYYY=year, MM=month, DD=day\n')  # start time of release
        fid.write(' ITIME2  =        ' + endDateStr[
                                         8:] + ', ! Release start time in UTC HHMISS: HH hours, MI=minutes, SS=seconds\n')  # end time of release
        fid.write(' LON1    =         ' + str(
            gridCorner1[iidx][0]) + ', ! Left longitude of release box -180 < LON1 <180 \n')  # Lower left latitude
        fid.write(' LON2    =         ' + str(
            gridCorner2[iidx][0]) + ', ! Right longitude of release box, same as LON1 \n')  # Lower left latitude
        fid.write(' LAT1    =           ' + str(
            gridCorner1[iidx][1]) + ', ! Lower latitude of release box, -90 < LAT1 < 90 \n')  # Lower left latitude
        fid.write(' LAT2    =           ' + str(
            gridCorner2[iidx][1]) + ', ! Upper latitude of release box same format as LAT1 \n')  # Lower left latitude
        fid.write(' Z1      =             ' + str(
            altitude[iidx][0]) + ', ! Lower height of release box meters/hPa above reference level\n')
        fid.write(' Z2      =             ' + str(
            altitude[iidx][1]) + ', ! Upper height of release box meters/hPa above reference level\n')
        fid.write(
            ' ZKIND   =              1, ! Reference level 1=above ground, 2=above sea level, 3 for pressure in hPa\n')
        fid.write(' MASS    =           1000, ! Total mass emitted, only relevant for fwd simulations\n')
        fid.write(' PARTS   =          ' + str(numParticles) + ', ! Total number of particles to be released\n')
        fid.write(' COMMENT =    "' + trajNameFinal + '", ! Comment, written in the outputfile\n')
        fid.write('/\n')
    fid.close  # close file
    return


def writeRELEASES(gridCorner1, gridCorner2, altitude, releaseDateStart, releaseDateEnd, trajName, numParticles,
                  writeFolder, releaseID, releaseNames):
    # Script to automate the writing of the RELEASES file for FLEXPART
    # backtrajectory calculations
    # gridCorner1 is a Nx2 matrix populated with pairs of the lower left corner of the
    # desired release box
    # gridCorner2 is a Nx2 matrix populated with pairs of the upper right corner of the
    # desired release box
    # Altitude is a Nx1 or Nx2 matrix populated with the release height or the upper and lower bounds of the
    # desired release box--should be in meters above sea level
    # Release date a Nx1 matrix of matlab date numbers corresponding to the
    # initialization of the backtrajectory
    # trajName is a Nx1 cell array corresponding to the desired trajectory
    # names
    # Each release is currently set to be one hour and 500 particles in a
    # domain filled box. These can be changed in the first few lines of code if
    # desired.
    # print header
    import os
    import matplotlib.dates as dt
    import numpy as np

    fid = open(os.path.join(writeFolder, 'RELEASES'), 'w')
    fid.write('*************************************************************************\n')
    fid.write('*                                                                       *\n')
    fid.write('*                                                                       *\n')
    fid.write('*                                                                       *\n')
    fid.write('*   Input file for the Lagrangian particle dispersion model FLEXPART    *\n')
    fid.write('*                        Please select your options                     *\n')
    fid.write('*                                                                       *\n')
    fid.write('*                                                                       *\n')
    fid.write('*                                                                       *\n')
    fid.write('*************************************************************************\n')
    fid.write('&RELEASES_CTRL\n')
    fid.write(' NSPEC      =          %d, ! Total number of species\n' % len(releaseID))
    fid.write(' SPECNUM_REL=          %s, ! Species numbers in directory SPECIES\n' % str(releaseID)[1:-1])
    fid.write('/\n')
    fid.write('&RELEASE\n')
    format1 = '%Y%m%d %H%M%S'
    format2 = '%Y%m%d_%H%M%S'
    for iidx, i in enumerate(gridCorner1):
        tmpRDS = releaseDateStart[iidx]
        if type(tmpRDS) != dt.datetime.datetime:
            tmpRDS = dt.num2date(tmpRDS)
        tmpRDE = releaseDateEnd[iidx]
        if type(tmpRDE) != dt.datetime.datetime:
            tmpRDE = dt.num2date(tmpRDE)
        startDateStr = dt.datetime.datetime.strftime(tmpRDS, format=format1)
        endDateStr = dt.datetime.datetime.strftime(tmpRDE, format=format1)
        trajNameFinal = trajName + '_' + dt.datetime.datetime.strftime(tmpRDE, format=format2) + '_' + str(
            iidx + 1).zfill(3)
        fid.write(' IDATE1  =       ' + startDateStr[
                                        0:8] + ', ! Release start date, YYYYMMDD: YYYY=year, MM=month, DD=day\n')  # start time of release
        fid.write(' ITIME1  =        ' + startDateStr[
                                         8:] + ', ! Release start time in UTC HHMISS: HH hours, MI=minutes, SS=seconds\n')  # end time of release
        fid.write(' IDATE2  =       ' + endDateStr[
                                        0:8] + ', ! Release start date, YYYYMMDD: YYYY=year, MM=month, DD=day\n')  # start time of release
        fid.write(' ITIME2  =        ' + endDateStr[
                                         8:] + ', ! Release start time in UTC HHMISS: HH hours, MI=minutes, SS=seconds\n')  # end time of release
        fid.write(' LON1    =      ' + str(
            gridCorner1[iidx][0]) + ', ! Left longitude of release box -180 < LON1 <180 \n')  # Lower left latitude
        fid.write(' LON2    =      ' + str(
            gridCorner2[iidx][0]) + ', ! Right longitude of release box, same as LON1 \n')  # Lower left latitude
        fid.write(' LAT1    =       ' + str(
            gridCorner1[iidx][1]) + ', ! Lower latitude of release box, -90 < LAT1 < 90 \n')  # Lower left latitude
        fid.write(' LAT2    =       ' + str(
            gridCorner2[iidx][1]) + ', ! Upper latitude of release box same format as LAT1 \n')  # Lower left latitude
        fid.write(' Z1      =             ' + str(
            altitude[iidx][0]) + ', ! Lower height of release box meters/hPa above reference level\n')
        fid.write(' Z2      =             ' + str(
            altitude[iidx][1]) + ', ! Upper height of release box meters/hPa above reference level\n')
        fid.write(
            ' ZKIND   =              1, ! Reference level 1=above ground, 2=above sea level, 3 for pressure in hPa\n')
        fid.write(' MASS    =         100000, ! Total mass emitted, only relevant for fwd simulations\n')
        fid.write(' PARTS   =          ' + str(numParticles) + ', ! Total number of particles to be released\n')
        fid.write(' COMMENT =    "' + releaseNames[iidx] + '", ! Comment, written in the outputfile\n\n')
    fid.write('/\n')
    fid.close  # close file
    return
