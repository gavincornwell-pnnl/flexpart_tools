def writeFLEXWRF_input(StartTime,EndTime,NumPart,Longitude,Latitude,NewInputFile,OutputPath,MetPath,Available,TemplateFile='./input.hrrr.template'):
    '''
    Function to write a flexpart-wrf file using a template file. The key parameters, are:
    StartTime: as python datetiem
    StopTime:  ""
    Numpart: number of particles to release
    Longitude/Latitude
    MetPath: folder where met files arelocated
    Available: file where the AVAILABLE file (telling FLEXWRF what met files correpsond to what times
    templatefile: file to use as a template, to capture most of the standard variable inputs
    '''

    import os
    import datetime as dt

    ### read in file
    fid = open(TemplateFile)
    lines = fid.readlines()
#    print(len(lines))
    fid.close()
    
    ### write paths
    if OutputPath[-1] != '/':
        OutputPath = OutputPath+'/'
    if MetPath[-1] != '/':
        MetPath = MetPath+'/'
    lines[1] = OutputPath + '\n'
    lines[2] = MetPath + '\n'
    lines[3] = Available + '\n'

    ### write start and end time
    StartTimeFix = StartTime - dt.timedelta(days=5) # go back five days
    lines[7] = '    %s    YYYYMMDD HHMISS   starting date of simulation\n' % StartTimeFix.strftime('%Y%m%d %H%M%S')
    lines[8] = '    %s    YYYYMMDD HHMISS   end date of simulation\n' % EndTime.strftime('%Y%m%d %H%M%S')
    
    #### release information
    lines[76] = '   %s   YYYYMMDD HHMISS   beginning date of simulation\n' % StartTime.strftime('%Y%m%d %H%M%S')
    lines[77] = '   %s   YYYYMMDD HHMISS   ending date of simulation\n' % EndTime.strftime('%Y%m%d %H%M%S')
    lines[78] = '   %s        XPOINT1 (real)  longitude [deg] of lower left corner\n' % str(Longitude)
    lines[79] = '   %s        XPOINT1 (real)  latitude [deg] of lower left corner\n' % str(Latitude)
    lines[80] = '   %s        XPOINT2 (real)  longitude [deg] of upper right corner\n' % str(Longitude)
    lines[81] = '   %s        YPOINT2 (real)  latitude [deg] of upper right corner\n' % str(Latitude)
    lines[85] = '  %s           NPART (int)    total number of particles to be released\n' % str(NumPart) 
    simname = os.path.basename(NewInputFile)
    lines[87] = '   %s       NAME OF RELEASE LOCATION^M\n' % simname 

    #### write out files
    fid = open(NewInputFile,'w')
    fid.writelines(lines)
    fid.close()
