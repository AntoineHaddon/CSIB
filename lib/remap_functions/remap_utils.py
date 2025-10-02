"""
These functions are used for remapping outputs from CanESM or other CMOR-like files.

Written by: Jonathan Izett (2024/2025)
"""

from argparse import Namespace
import glob
import math
from netCDF4 import Dataset
import numpy as np
import os
import re
import shlex
import subprocess
import sys
from typing import List, Optional, Tuple, Union
import xarray as xr

#####################
# GENERAL PROCESSING
#####################
def calc_nemo_chunk_dates(args: Namespace) -> None:
    """
    Compute the start and end date of the current chunk
    Copied and modified from: CanNEMO/lib/maestro/bin/compute_iteration_info.py

    Input:
        args - args passed to remap_canesm.py
    Output:
        Writes date info. to nemo_counter_info.cfg
    """
    nf = int(args.nemo_freq_months)

    if int(args.run_start_month) not in range(1,13):
        raise ValueError(f'Start month must be 1-12 ({args.run_start_month})')
 
    sm = int(args.run_start_month)
    sy = int(args.run_start_year)
    
    # get loop
    if args.loop == 0:
        # from directory name
        ll=int(os.path.split(os.getcwd())[-2].split('+')[-1])-1
    else:
        # passed as argument
        ll = int(args.loop)-1

    cl_start_nmonth = sm + ll*nf
    cl_end_nmonth = sm + (ll+1)*nf - 1

    # dictionary mapper to get days in month
    days_in_month = {1 : 31,
                 2 : 28,
                 3 : 31,
                 4 : 30,
                 5 : 31,
                 6 : 30,
                 7 : 31,
                 8 : 31,
                 9 : 30,
                 10 : 31,
                 11 : 30,
                 12 : 31
                 } 

    # Get the calendar start/end month for each loop segment
    cl_start_cal_month = (cl_start_nmonth)%12 if cl_start_nmonth%12 > 0 else 12
    cl_end_cal_month =  cl_end_nmonth%12 if cl_end_nmonth%12 > 0 else 12
    # Get the calendar start/end year for each loop segment
    cl_start_cal_year = int((sm + ll*nf-1)/12 + sy) 
    cl_end_cal_year = int(math.ceil((sm + (ll+1)*nf -1)/12.0)) + sy -1

    cl_start_cal_day = 0o1 
    cl_end_cal_day = days_in_month[cl_end_cal_month]
    
    # List of all years in this chunk
    chunk_years = range(cl_start_cal_year, cl_end_cal_year+1,1)
    chunk_years_str = " ".join('%04d' % year for year in chunk_years)

    if args.outfile is None:
        cfgFile='nemo_counter_info.cfg'
    else:
        cfgFile=os.path.join(os.path.dirname(os.path.abspath(args.outfile)),'nemo_counter_info.cfg')

    with open(cfgFile, 'w') as ff: 
        ff.write('NEMO_CHUNK_START_DAY=%s\n' % (cl_start_cal_day))
        ff.write('NEMO_CHUNK_START_MONTH=%s\n' % (cl_start_cal_month))
        ff.write('NEMO_CHUNK_START_YEAR=%s\n' % (cl_start_cal_year))
        ff.write('NEMO_CHUNK_END_DAY=%s\n' % (cl_end_cal_day))
        ff.write('NEMO_CHUNK_END_MONTH=%s\n' % (cl_end_cal_month))
        ff.write('NEMO_CHUNK_END_YEAR=%s\n' % (cl_end_cal_year))
        ff.write(f'NEMO_CHUNK_START_DATE={cl_start_cal_year:04}-{cl_start_cal_month:02}-{cl_start_cal_day:02}\n')
        ff.write(f'NEMO_CHUNK_END_DATE={cl_end_cal_year:04}-{cl_end_cal_month:02}-{cl_end_cal_day:02}\n')
        ff.write("NEMO_CHUNK_YEARS='%s'\n" % (chunk_years_str))
        ff.write("NEMO_MODEL_LOOP='%s'\n" % (ll+1))

    return

def parseArgYears(args: Namespace) -> List[int]:
    """
    Parse the -y argument to determine if single year, multiple years, or range of years.

    Input:
        args - arguments passed to remap_canesm.py (uses args.years)
    Output:
        years - returns list of all years desired

    """
    if len(args.years)==2 and (max(args.years)-min(args.years)) > 1:
        years=range(min(args.years),max(args.years)+1)
    elif len(np.unique(args.years))==1:
        years=[np.unique(args.years)[0]]
    else:
        years=np.array(args.years)
    return years

def matchFileYear(year: int, pathPattern: str, file0: Optional[str] = "") -> Tuple[Optional[List[str]], Optional[str], bool]:
    """
    Match the desired year to the file.
    
    Necessary because not all files will have the same year format.
    Inputs:
        year - desired year (int)
        pathPattern - path/pattern to search for files using glob (str)
        file0 - name of original file checked (str)
    Outputs:
        file - file name that matches desired year (str)
        fdates - datestring of original file (str)
        sameFile - logical; True if file found is same as file0
    """
    # initialize return values
    foundFile=False; sameFile=False; file=[]

    # get list of all files that match desired pattern
    flist=np.array(sorted(glob.glob(pathPattern)))

    # find file with either this year or a previous year in the name
    # not the most elegant, but it is pretty quick!
    # allows multiple years in a single file, potentially speeding up the file
    # generation for more than one desired year
    for fle in flist:
        # look for year in the file name
        if f'_{year:04}' in fle or f'-{year:04}' in fle:
            file.append(fle)
    if len(file) > 0:
        foundFile=True
        sameFile=file[0]==file0
        fdates=file[0].split('_')[-1].split('-')[0]+'-'+file[-1].split('_')[-1].split('-')[1]
    else:
        for tryYear in range(year-1,-1,-1):
            # check for a file that starts before the desired year and ends on or after the desired year
            for fle in flist:
                if f'_{tryYear:04}' in fle:
                    # check that a later year is also in the file name (or year at the end)
                    for tryYear2 in range(year,year+201):
                        if f'-{tryYear2:04}' in fle:
                            file=[fle]
                            fdates=os.path.basename(fle).split('_')[-1]
                            foundFile=True
                            break
                    if foundFile:
                        break
    if foundFile:
        sameFile=file0==file
        return file,fdates,sameFile
    else:
        # let user know no files were found and give and indication of available years
        print(f"No files found (Year: {year})")
        if len(flist) > 0:
            fyears=[]
            for fle in flist:
                # get dates from file name
                fyears.append([int(yy[0:-2]) for yy in re.search("([0-9]{6}\-[0-9]{6})", fle)[0].split('-')])
            print(fyears)
            print(f"  The closest year appears to be {fyears[np.argmin(np.array(fyears)-year)]}")
            print(f"  Range: {np.nanmin(fyears)} to {np.nanmax(fyears)}")
        else:
            print('No files available.')
        return None,None,False

def checkMeshFile(meshFile: str) -> Tuple[float, float]:
    """
    Check/modify mesh file to ensure it is in proper format to use cdo

    Input:
        meshFile - path to netCDF mesh file (str)
    Output:
        Writes a temporary file: grd.tmp.nc with cdo-compliant coordinates
        bN - northernmost point of grid
        bS - southernmost point of grid
    """
    with xr.open_dataset(meshFile) as ds0:   
        try:
            x=len(ds0['x']); y=len(ds0['y'])
        except:
            x=np.shape(ds0['nav_lat'].values)[1]; y=len(ds0['nav_lat'])
        # need some 'data'
        nav_ones=np.ones((1,y,x))
        # convert to dataset
        msh=xr.Dataset.from_dict(
            {'nav_lon':{'dims':('y','x'),'data':ds0['nav_lon'].values,'attrs':{'_CoordinateAxisType':'Lon','units':'degrees_east'}},
            'nav_lat':{'dims':('y','x'),'data':ds0['nav_lat'].values,'attrs':{'_CoordinateAxisType':'Lat','units':'degrees_north'}},
            'nav_ones':{'dims':('time_counter','y','x'),'data':nav_ones,'attrs':{'coordinates':'nav_lat nav_lon'}},
            'time_counter':{'dims':('time_counter'),'data':[0.]}})
        # write to temporary netCDF file
        msh.to_netcdf(f'grd.tmp.nc')
    with xr.open_dataset('grd.tmp.nc') as ds0:
        # latitude range of grid
        bN=np.nanmax(ds0['nav_lat'])
        bS=np.nanmin(ds0['nav_lat'])
    return bN,bS

def getZ(meshFile: str, outFile: str) -> None:
    """
    Get a file with just z coordinates.

    Inputs:
        meshFile - path to netCDF file containing grid coordiantes (str)
        outFile - path/name of desired output file (str)
    Output:
        Writes outFile with only vertical coordinate variables
    """
    if not os.path.isfile(f'{outFile}.onlyz.tmp.nc'):
        with xr.open_dataset(meshFile) as mF:
            if 'e3t_0' in mF.keys():
                err=subproc(f"ncks -h -v e3t_0 {meshFile} -O {outFile}.onlyz.tmp.nc")
            elif 'tmask' in mF.keys():
                err=subproc(f"ncks -h -v tmask {meshFile} -O {outFile}.onlyz.tmp.nc")
            elif 'gdept_0' in mF.keys():
                err=subproc(f"ncks -h -v gdept_0 {meshFile} -O {outFile}.onlyz.tmp.nc")
            elif 'votemper' in mF.keys():
                err=subproc(f"ncks -h -v votemper {meshFile} -O {outFile}.onlyz.tmp.nc")
            else:
                sys.exit('Error: cannot find depth coordinate.')
    return

def intZ(inFile: str, outFile: str, zFile: str) -> None:
    """
    Interpolate in z-coordinates.

    Inputs:
        inFile - original path/name of file to be interpolated (str)
        outFile - path/name of desired output file (str)
        zFile - file containing z coordiantes (e.g., output from getZ) (str)
    Output:
        Writes outFile interpolated to zFile coordinates
    """
    
    # have to use nested CDO command, so write to shell script
    cmd=f"cdo -w --no_history intlevelx$(cdo -s showlevel {zFile} | tr ' ' ',') {inFile} {outFile}"
    with open("intz.tmp.sh","w") as shFile:
        shFile.write(f"#!/bin/bash\n\nsource ~/.profile\n\n{cmd}")
    err=subproc("chmod +x intz.tmp.sh")
    err=subproc("./intz.tmp.sh")          # run temporary shell script
    err=subproc("rm -f intz.tmp.sh")    # remove temporary shell script
    
    return

def cellAreas(meshFile: str) -> np.ndarray:
    """
    Calculate grid cell area from a given mesh

    Input:
        meshFile - path to file containing e1t and e2t variables
                   (grid cell dimensions)
    Output:
        gridArea - area of grid cells (float)
    """
    with xr.open_dataset(meshFile) as mF:
        e1t = mF['e1t'].values.squeeze() # y,x
        e2t = mF['e2t'].values.squeeze() # y,x

        # calulate grid area
        gridArea = e1t*e2t

    return gridArea

def splitCMD(cmd: str) -> List[str]:
    """
    Split a command, preserving values in quotations.

    Input:
        cmd - string of desired command
    Output:
        Returns string split into a list, preserving quotation marks
    """
    if "\"" in cmd or "(" in cmd or "'" in cmd:
        # remove leading/trailing spaces
        cmd=cmd.strip(" ")
        # split, preserving brackets and quotation marks
        nb_brackets=0; nb_quotes=0 ; nb_quote=0
        l=[0]
        for iC,cC in enumerate(cmd):
            if cC=="(":
                nb_brackets+=1
            elif cC==")":
                nb_brackets-=1
            elif cC=="\"" and nb_quotes==0:
                nb_quotes+=1
            elif cC=="\"" and nb_quotes==1:
                nb_quotes-=1
            elif cC=="'" and nb_quote==0:
                nb_quote+=1
            elif cC=="'" and nb_quote==1:
                nb_quote-=1
            elif cC==" " and nb_brackets==0 and nb_quotes==0 and nb_quote==0:
                l.append(iC)
        l.append(len(cmd))
        qbsplit=[cmd[i:j].strip(" ") for i,j in zip(l,l[1:])]
        return qbsplit
    else:
        return cmd.split(" ")

def subproc(cmd: Union[str, List[str]]) -> bytes:
    """
    Run a subprocess command in shell and check if successful.

    Input:
        cmd - list of strings to pass as command to subprocess
              If not a list, splits by spaces; unless it is the
              path to a shell script, in which case it runs the script.
              All options avoid invoking shell=True
    Runs:
        sbpr=subprocess.run(cmd.split(' '),capture_output=True)
    """

    if type(cmd) is list or (".sh" in cmd and " " not in cmd):
        # list or shell script passed
        sbpr=subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        # split by spaces, preserving quotation marks
        sbpr=subprocess.run(splitCMD(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if sbpr.returncode > 0:
        raise RuntimeError(f"\n{cmd}\n\n{sbpr.stdout.decode()}\n{sbpr.stderr.decode()}")

    return sbpr.stderr

def subscript(cmd: str) -> bytes:
    """
    Create a temporary shell script and run it.

    Input:
        cmd - string command 
    """

    with open("tmp.sh","w") as shFile:
        shFile.write(f"#!/bin/bash\n\nsource ~/.profile\n\n{cmd}")
    err=subproc("chmod +x tmp.sh")
    err=subproc("./tmp.sh")          # run temporary shell script
    err=subproc("rm -f tmp.sh")    # remove temporary shell script

    return

def changeCoords(inFile: str, outFile: str) -> None:
    """
    Use ncatted to change coordinate name and ensure CDO-compliant.

    Inputs:
        inFile: input file name
        outFile: file name for output
    Output:
        Writes to output file name
    """

    err=subproc(f"ncatted -h -a _CoordinateAxisType,nav_lon,o,c,Lon -a units,nav_lon,o,c,degrees_east -a _CoordinateAxisType,nav_lat,o,c,Lat -a units,nav_lat,o,c,degrees_north {inFile} -O {outFile}")

    return

def addCoords(varName: str, fName: str) -> None:
    """
    Add nav_lat and nav_lon coordinates to netCDF variable.

    Inputs:
        varName - variable name to add coordinates
        fName - path to file to add coordinate info (file is change in place)
    """
    
    err=subproc(f"ncatted -h -a coordinates,{varName},o,c,\"nav_lat nav_lon\" {fName} -O {fName}")

    return

def selYear(year1: int, year2: int, inFile: str, outFile: str) -> None:
    """
    Use cdo selyear to getonly desired years in file.

    Inputs:
        year1/year2 - desired years (start and end of range)
        inFile - path to input file
        outFile - path to output file
    Output:
        Writes new file at outFile with desired year.
    """
    
    err=subproc(f"cdo --no_history selyear,{year1}/{year2} {inFile} {outFile}")

    return

def selBox(lon1: float, lon2: float, lat1: float, lat2: float, inFile: str, outFile: str) -> None:
    """
    Use cdo sellonlatbox to restrict the file domain.

    Inputs:
        lon1 - Westernmost longitude
        lon2 - Easternmost longitude
        lat1 - Southernmost latitude
        lat2 - Northernmost latitude
        inFile - path to input file
        outFile - path to output file
    Output:
        Writes new file at outFile with restricted domain
    """
    
    err=subproc(f"cdo --no_history sellonlatbox,{lon1},{lon2},{lat1},{lat2} {inFile} {outFile}")

    return

def fillMiss2(inFile: str, outFile: str) -> None:
    """
    Use cdo fillmiss2,2 to fill any gaps in file.

    Inputs:
        inFile - path to input file
        outFile - path for writing new filled file
    Output:
        writes filled file to outFile
    """

    err=subproc(f"cdo -w --no_history fillmiss2,2 {inFile} {outFile}")

    return

def remap(grdFile: str, inFile: str, outFile: str) -> None:
    """
    Use cdo remapdis to remap from input file to a new grid.

    Inputs:
        grdFile - path grid file to which you want to remap
        inFile - path to original file on original grid
        outFile - path for writing new file remapped to grdFile
    Output:
        Writes remapped field to outFile.
    """

    err=subproc(f"cdo -w --no_history remapdis,{grdFile} {inFile} {outFile}")

    return

def cleanTmp(pattern: str) -> None:
    """
    Clean all files matching the desired pattern.
    """
    rmFiles="rm -f"
    for fF in glob.glob(pattern):
        rmFiles+=f" {fF}"

    err=subproc(rmFiles)

    return

def his2cmor(var4cmor: str, fYear: int, fOut: str, args: Namespace) -> str:
    """
    Naive conversion from history files to a "CMOR-like" file (i.e., changes variable names, coordinates, etc.)

    Inputs:
        var4cmor - variable name (str)
        fYear - desired file year (int)
        fOut - output file name/path for new "CMOR-like" file (str)
        args - args passed to remap_canesm.py
    Output:
        Writes a new "CMOR-like" file to fOut
    """
    if args.parent_name.lower()=='glorys':
        ftype={'thetao':'grid_[tT]','so':'grid_[tT]','zos':'grid_[tT]','uo':'grid_[tT]','vo':'grid_[tT]'}
    else:
        ftype={'thetao':'grid_[tT]','so':'grid_[tT]','zos':'grid_[tT]','uo':'grid_[uU]','vo':'grid_[vV]'}
    lev={'grid_[tT]':'deptht','grid_[uU]':'depthu','grid_[vV]':'depthv'}

    # find correct file
    pPath=args.parent_path.split('nc_output')[0]
    varPath=os.path.join(pPath,f"mc_{args.parent_name.replace('CanESM5-NEMO4-','')}_{fYear}_*{ftype[var4cmor]}.nc*")
    
    # find all files that match format
    flist=np.array(sorted(glob.glob(varPath)))
    if len(flist)==0 and 'glorys' in args.parent_name.lower():
      varPath=os.path.join(pPath,f"cmems_mod_glo_*{fYear}-01-01-{fYear}-12-01.nc")
      print(varPath)
      flist=np.array(sorted(glob.glob(varPath)))

    if len(flist) == 0:
        cmfile=None
    else:
        # get highest value of file version
        hisFile=flist[-1]
        cmfile=f"{fOut.replace('VAR',var4cmor)}.{var4cmor}.2cmor.tmp.nc"
        
        # slice file to get desired variable
        if 'glorys' in args.parent_name.lower() and args.type==1:
            # restrict latitude further as files are BIG??
            if args.type==1:
              err=subproc(f"cdo -w --no_history sellonlatbox,-180,180,30,45 {hisFile} {cmfile}")
              err=subproc(f"ncks -v {var4cmor} {cmfile} -O {cmfile}")
        else:
            err=subproc(f"ncks -v {var4cmor} {hisFile} -O {cmfile}")
        # rename coordinates (if applicable)
        err=subproc(f"ncrename -h -d .x,i -d .y,j -d .longitude,i -d .latitude,j -v .nav_lon,longitude -v .nav_lat,latitude {cmfile} -O {cmfile}")
        # rename vertical coordinate (if applicable)
        # either a grid-specific name, or just generic 'depth'
        err=subproc(f"ncrename -h -d .{lev[ftype[var4cmor]]},lev {cmfile} -O {cmfile}")
        err=subproc(f"ncrename -h -d .depth,lev {cmfile} -O {cmfile}")
        try:
          err=subscript(f"ncrename -h -a coordinates,{var4cmor},o,c,'latitude longitude' {cmfile} -O {cmfile}")
        except:
          err=subscript(f"ncatted -h -a coordinates,{var4cmor},o,c,'latitude longitude' {cmfile} -O {cmfile}")
        # add gridded lat/lon if necessary
        with xr.open_dataset(cmfile) as cmtmp:
            shp=np.shape(cmtmp.latitude.values)
            try:
                shpi=shp[1]
                shpj=shp[0]
                reLat=False
            except:
                ln=cmtmp.longitude.values
                lt=cmtmp.latitude.values
                LN,LT=np.meshgrid(ln,lt)
                shpi=len(ln)
                shpj=len(lt)
                reLat=True
        #add dimension variables to file         
        with xr.open_dataset(cmfile) as cmtmp:
            if var4cmor=='zos':
                ij=xr.Dataset.from_dict(
                    {'i':{'dims':('i'),'data':range(shpi),'attrs':{'units':'1','long_name':'first spatial index for variables stored on an unstructured grid'}},
                    'j':{'dims':('j'),'data':range(shpj),'attrs':{'units':'1','long_name':'second spatial index for variables stored on an unstructured grid'}}})
            else:
                with xr.open_dataset(hisFile) as hstmp:
                    try:
                        depth=hstmp.variables[lev[ftype[var4cmor]]]
                    except:
                        try:
                            depth=hstmp.variables['lev']
                        except:
                            depth=hstmp.variables['depth']
                    ij=xr.Dataset.from_dict(
                        {'i':{'dims':('i'),'data':range(shpi),'attrs':{'units':'1','long_name':'first spatial index for variables stored on an unstructured grid'}},
                        'j':{'dims':('j'),'data':range(shpj),'attrs':{'units':'1','long_name':'second spatial index for variables stored on an unstructured grid'}},
                        'lev':{'dims':('lev'),'data':depth,'attrs':{'long_name':"Vertical T levels",'units':'m','positive':"down",'axis':"Z"}}})
            cmtmp2=xr.merge([cmtmp,ij])
            cmtmp2.to_netcdf(f'{cmfile}.ij.tmp.nc')
            if 'zos' not in cmfile:
              err=subproc(f"ncatted -h -a positive,lev,o,c,down {cmfile}.ij.tmp.nc -O {cmfile}.ij.tmp.nc")
              err=subproc(f"ncatted -h -a axis,lev,o,c,Z {cmfile}.ij.tmp.nc -O {cmfile}.ij.tmp.nc")
        if reLat:
            # want lat/lon to be gridded, not vectors
            # remove lon & lat from file
            err=subproc(f"ncks -v {var4cmor} {cmfile}.ij.tmp.nc {cmfile}.nolonlat.tmp.nc")
            lnlt=xr.Dataset.from_dict(
                    {'longitude':{'dims':('j','i'),'data':LN,'attrs':{'_CoordinateAxisType':'Lon','units':'degrees_east'}},
                     'latitude':{'dims':('j','i'),'data':LT,'attrs':{'_CoordinateAxisType':'Lat','units':'degrees_north'}}})
            lnlt.to_netcdf(f'{cmfile}.lonlat.tmp.nc')
            err=subproc(f"cdo merge {cmfile}.nolonlat.tmp.nc {cmfile}.lonlat.tmp.nc {cmfile}.gridded.tmp.nc")
            err=subscript(f"ncatted -h -a coordinates,{var4cmor},o,c,'latitude longitude' {cmfile}.gridded.tmp.nc -O {cmfile}.gridded.tmp.nc")
            cmfile=[f'{cmfile}.gridded.tmp.nc']
        else:
            cmfile=[f'{cmfile}.ij.tmp.nc']

    return cmfile

def fixVelocities(vFile: str, outVFile: str) -> str:
    """
    Flip meridional velocities in the "northern" portion of a domain (if applicable)
    
    Inputs:
        vFile - file name/path containing velocities (str)
        outVFile - file name/path with flipped velocities (str)
    Output:
        Writes flipped velocities to outVFile

    NOTE: This is only relevant for, e.g., the CREG domain, which has an inverted Pacific Ocean
    """
    with xr.open_dataset(vFile) as vds:
        if 'vo' in vds.keys() or 'uo' in vds.keys():
            # history or CMORized file
            vV='vo'; uU='uo'
            if 'nav_lon' in vds.keys():
                # history file
                lt='nav_lat'; ydim='y'
            else:
                lt='latitude'; ydim='j'
            velFile=True
        elif 'vn' in vds.keys():
            # restart file
            vV='vn'; uU='un'
            lt='nav_lat'; ydim='y'
            velFile=True
        else:
            velFile=False
            outVFile=vFile
            print(f'Velocity not found in {vFile}.')

        if velFile:
            # check if sorted order is identical to original order
            vlat=vds[lt].values
            vInFile=False ; uInFile=False
            try:
                # may only have v in file
                vlat[np.isnan(vds[vV].values[0,0,:,:])]=np.nan
                vInFile=True
            except:
                # may only have u in file
                vlat[np.isnan(vds[uU].values[0,0,:,:])]=np.nan
                uInFile=True
            iX=np.argmax(np.nanmax(vlat,axis=0)) ; vlat=vlat[:,iX]
            vy=vds[ydim].values
            if np.nansum(np.abs(vy-np.argsort(vlat)))!=0 and np.nanmax(vlat) > 70:
                # there is a transition over the pole if latitudes are not monotonic
                # i.e., sorting doesn't match index order
                # restrict to only search near pole
                i90=vlat>80
                vy=vy[i90]
                vlat=vlat[i90]
                
                # find transition by looping over latitudes from top
                iTransY=None
                for iY,VY in enumerate(vy):
                    # check if latitude is larger than previous
                    # if not, then assume transition point and break
                    if vlat[iY] < vlat[iY-1]:
                        iTransY=VY
                        break
                if iTransY is not None:
                    # create mask with negative signs in "northern" portion
                    if vInFile:
                        nMask=np.ones(np.shape(vds[vV].values))
                    else:
                        nMask=np.ones(np.shape(vds[uU].values))
                    nMask[:,:,iTransY::,:]=-1
                    # replace velocities with those multiplied by negative sign
                    if vInFile:
                        v2=nMask*vds[vV].values
                    try:
                        # u not always with v in file
                        u2=nMask*vds[uU].values
                        uInFile=True
                    except:
                        pass
                    # save velocities to new file
                    vds2=vds.copy()
                    if vInFile:
                        vds2[vV][:]=v2
                    if uInFile:
                        vds2[uU][:]=u2
                    if vV == 'vn':
                        # also adjust vb and 'ub'
                        if vInFile:
                            vds2['vb'][:]=nMask*vds['vb'].values
                        if uInFile:
                            vds2['ub'][:]=nMask*vds['ub'].values
                    vds2.to_netcdf(outVFile)
                    print(f'Saved {outVFile}')
                else:
                    outVFile=vFile

    return outVFile
