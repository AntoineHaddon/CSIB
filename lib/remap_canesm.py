"""
Functions to remap CMORIZED outputs from CanESM runs onto a different domain
(including CanTODS)

Written by J. G. Izett (2024)
"""

#TODO: Add rivers
#TODO: hot start

import argparse
import glob
import math
from netCDF4 import Dataset
import numpy as np
import os
import subprocess
import sys
import time
import xarray as xr

# argument parser
# allows to be calculated at runtime, or a posteriori
parser=argparse.ArgumentParser(description='Calculate CanTODS diagnostics.')

# add arguments
parser.add_argument('-y','--years',help='Years to process. It two years passed, processes years in range(year1,year2+1). If not two years, processes as a list.',action='append',required=True,type=int)
parser.add_argument('-P','--parent_path',help='Full path to files containing model output from parent run, ending before the ensemble identifier.\ne.g., /fs/site5/eccc/crd/ccrn/model_output/CMIP6/final/CMIP6/CMIP/CCCma/CanESM5',required=True)
parser.add_argument('-p','--parent_name',help='Name of parent run to help identify output files. Required for finding files.\ne.g., CanESM5',required=True)
parser.add_argument('-e','--parent_ensemble',help='Ensemble identifier for parent run.',required=True)
parser.add_argument('-x','--parent_experiment',help='Experiment (e.g., piControl) for parent run.',required=True)
parser.add_argument('-g','--parent_grid',help='Parent meshfile. Only needed if remapping rivers.',default=None)
parser.add_argument('-o','--outfile',help='Prefix for output file. If None, default name is created based on type of file being produced.',default=None)
parser.add_argument('-m','--meshfile',help='Meshfile for remapping.',default='CREG_NEMO4_domain_meshfile_forRemapping.nc')
parser.add_argument('-t','--type',help='Type of files to produce. Either:\n  -1 - nemo configuration info\n  0 - Initial conditions (default)\n  1 - Boundary conditions\n  2 - Atmospheric forcing (does not remap; only gets correct file time)',type=int,default=0)
parser.add_argument('-X','--Xdeg',help='Slice to contain X degrees either side of the boundary. Should be larger than the parent grid resolution to guarantee border. Default = 1.25.',default=1.25,type=float)
parser.add_argument('-B','--bdyFile',help='Map to a specific boundary coordinate file (if type=1), e.g., to the Med. Otherwise maps to presumed north/south boudary.',default=None)
parser.add_argument('-F','--forcing',help='Type of forcing CanESM (default) or OMIP',default='CanESM')
parser.add_argument('-i','--ic_ind',help='Index in file of desired time for initial condition. Default: 0',default=0)
parser.add_argument('-H','--hot_start',help='If remapping IC, generate hotstart file from NEMO restart. Default: 0',default=0)
parser.add_argument('-V','--flip_vel',help='Flip velocities in the "north" of the domain if regional fold',default=0)
parser.add_argument('-R','--run_start_year',help='Run start year if type==-1',default=0)
parser.add_argument('-r','--run_start_month',help='Run start month if type==-1',default=0)
parser.add_argument('-l','--loop',help='Sequencer loop if type==-1',default=0)
parser.add_argument('-f','--nemo_freq_months',help='If type==-1',default=0)
parser.add_argument('-M','--mor',help='CMOR frequency/directory (e.g., 3hr or Amon) for atmospheric forcing.',default='Amon')
parser.add_argument('-v','--rvr',help='River file for remapping/scaling of frehswater inputs (if type==3)',default=None)
parser.add_argument('-A','--iaf_year_offset',help='Year offset if wanting to use cyclical forcing.',default=None)
parser.add_argument('-a','--iaf_loop_year',help='Reference year for loop if offset for cyclical forcing.',default=None)
parser.add_argument('-c','--his2cmor',help='If 1, look for history file rather than CMORized file. Default: 0',default=0)

#-----------#
# FUNCTIONS #
#-----------#

def matchFileYear(year,flist,file0=''):
    """
    Match the desired year to the file. Necessary because not all files will have
    the same year format.
    """
    foundFile=False; sameFile=False; file=[]
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
        return None,None,False

def checkMeshFile(meshfile):
    # ensure mesh file is in proper format to be able to remap using cdo
    with xr.open_dataset(meshfile) as ds0:   
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
    return

def getZ(meshFile,outFile):
    """
    Get a file with just z coordinates.
    """
    if not os.path.isfile(f'{outFile}.onlyz.tmp.nc'):
        with xr.open_dataset(args.meshfile) as mF:
            if 'e3t_0' in mF.keys():
                subprocess.run(f'ncks -h -v e3t_0 {args.meshfile} -O {outFile}.onlyz.tmp.nc',shell=True)
            elif 'tmask' in mF.keys():
                subprocess.run(f'ncks -h -v tmask {args.meshfile} -O {outFile}.onlyz.tmp.nc',shell=True)
            elif 'gdept_0' in mF.keys():
                subprocess.run(f'ncks -h -v gdept_0 {args.meshfile} -O {outFile}.onlyz.tmp.nc',shell=True)
            elif 'votemper' in mF.keys():
                subprocess.run(f'ncks -h -v votemper {args.meshfile} -O {outFile}.onlyz.tmp.nc',shell=True)
            else:
                sys.exit('Error: cannot find depth coordinate.')
    return

def cellAreas(meshFile):
    with xr.open_dataset(meshFile) as mF:
        e1t = mF['e1t'].values.squeeze() # y,x
        e2t = mF['e2t'].values.squeeze() # y,x

        # calulate grid area
        gridArea = e1t*e2t

    return gridArea

def rs_remap(args):
    """
    Remap restart files for hot start.
    """

    # loop through each variable
    print(f'Finding and processing data from {args.parent_name}_{args.parent_experiment}_{args.parent_ensemble}.')

    # get mesh file in correct format and find lat/lon bounds
    checkMeshFile(args.meshfile)
    with xr.open_dataset('grd.tmp.nc') as ds0:
        # latitude range of grid
        bN=np.nanmax(ds0['nav_lat'])
        bS=np.nanmin(ds0['nav_lat'])

    if args.outfile is None:
        outFile='TYPE'
    else:
        outFile=args.outfile

    # find directory containing restart files
    filePath=np.array(sorted(glob.glob(os.path.join(args.parent_path,f'mc_{args.parent_name}_{min(args.years)}_m*_nemors.*'))))
    if len(filePath) > 1:
        # do not want tiled files
        iP = np.array(['tiled' not in fP for fP in filePath])
        filePath=filePath[iP]
    if len(filePath) > 1:
        # only keep the highest suffix
        sfxs=np.array([int(fP[-3::]) for fP in filePath])
        iP = np.array([f'nemors.{max(sfxs):03}' in fP for fP in filePath])
        filePath=filePath[iP]
    if len(filePath) > 1:
        # only keep the correct month
        iP = np.array([f'm{args.ic_ind+1:02}_nemors' in fP for fP in filePath])
        filePath=filePath[iP]
    if len(filePath) != 1:
        sys.exit('No restart files found!')

    # loop through different restart file types
    for tT in ['restart','restart_ice','restart_trc']:
        # ensure coordinates are correctly assigned to variable metadata (for remapping)
        subprocess.run(f"ncatted -h -a _CoordinateAxisType,nav_lon,o,c,Lon -a units,nav_lon,o,c,degrees_east -a _CoordinateAxisType,nav_lat,o,c,Lat -a units,nav_lat,o,c,degrees_north {os.path.join(filePath[0],f'{args.parent_name}*_{tT}.nc')} -O {outFile.replace('TYPE',tT)}.atts.tmp.nc",shell=True)
        # loop through variables in file and assign coordinates as appropriate
        geoVars=[]; scalars=[]
        with xr.open_dataset(f"{outFile.replace('TYPE',tT)}.atts.tmp.nc") as attFle:
            for vV in attFle.keys():
                if vV not in ['nav_lon','nav_lat'] and len(np.shape(attFle.variables[vV]))>2:
                    # keep track of variables with geographic coordinates to modify metadata
                    geoVars.append(vV)
                elif len(np.shape(attFle.variables[vV]))==0:
                    # keep track of scalars to add back to file at end
                    scalars.append(vV)
        for gV in geoVars:
            subprocess.run(f"ncatted -h -a coordinates,{gV},o,c,\"nav_lat nav_lon\" {outFile.replace('TYPE',tT)}.atts.tmp.nc -O {outFile.replace('TYPE',tT)}.atts.tmp.nc",shell=True)

        # subsample srcFile to be within +/- X degrees of southernmost point to speed things up
        subprocess.run(f"cdo --no_history sellonlatbox,-180,180,{np.nanmax([-90,bS-args.Xdeg])},{np.nanmin([90,bN+args.Xdeg])} {outFile.replace('TYPE',tT)}.atts.tmp.nc {outFile.replace('TYPE',tT)}.sliced.tmp.nc",shell=True)
        
        # fill any gaps in the sliced srcFile (two iterations)
        subprocess.run(f"cdo --no_history fillmiss2,2 {outFile.replace('TYPE',tT)}.sliced.tmp.nc {outFile.replace('TYPE',tT)}.filled.tmp.nc",shell=True)

        # remap to child grid
        subprocess.run(f"cdo --no_history remapdis,grd.tmp.nc {outFile.replace('TYPE',tT)}.filled.tmp.nc {outFile.replace('TYPE',tT)}.remapped.tmp.nc",shell=True)

        # interpolate vertically if appropriate
        if args.outfile is None and tT != 'restart':
            tSuff='_in'
        else:
            tSuff=''
        # Interpolate vertical levels...doesn't seem to work properly?! For now, leave out since grids all eORCA
        # with same vertical resolution.
        # if tT == 'restart_ice':
        #     # No vertical interpolation for ice files
        #     subprocess.run(f"mv {outFile.replace('TYPE',tT)}.remapped.tmp.nc {outFile.replace('TYPE',tT)}.tmp.zmap.nc",shell=True)
        # else:
        #     # interpolate vertical levels
        #     getZ(args.meshfile,outFile)
        #     subprocess.run(f"cdo --no_history intlevelx$(cdo -s showlevel {outFile}.onlyz.tmp.nc | tr ' ' ',') {outFile.replace('TYPE',tT)}.remapped.tmp.nc {outFile.replace('TYPE',tT)}.tmp.zmap.nc",shell=True)

        # add scalars back to file
        subprocess.run(f"cp {outFile.replace('TYPE',tT)}.remapped.tmp.nc {outFile.replace('TYPE',tT)}{tSuff}.nc",shell=True)
        for sS in scalars:
            subprocess.run(f"ncks -h -v {sS} {outFile.replace('TYPE',tT)}.atts.tmp.nc -A {outFile.replace('TYPE',tT)}{tSuff}.nc",shell=True)

        if tT == 'restart' and int(args.flip_vel)==1:
            fixVelocities(f"{outFile.replace('TYPE',tT)}{tSuff}.nc",f"{outFile.replace('TYPE',tT)}{tSuff}.nc.velflip.nc")
            subprocess.run(f"rm -f {outFile.replace('TYPE',tT)}{tSuff}.nc ; mv {outFile.replace('TYPE',tT)}{tSuff}.nc.velflip.nc {outFile.replace('TYPE',tT)}{tSuff}.nc",shell=True)

        # remove intermediate variable files
        subprocess.run(f"rm -f {outFile.replace('TYPE',tT)}*.tmp.nc*",shell=True)
    # remove intermediate files
    subprocess.run(f'rm -f grd.tmp.nc {outFile}*.tmp.nc',shell=True)

def ic_remap(args):
    """
    Remap initial condition for cold start.
    """

    # get month index as integer
    args.ic_ind=int(args.ic_ind)

    # loop through each variable
    print(f'Finding and processing data from {args.parent_name}_{args.parent_experiment}_{args.parent_ensemble}.')

    # get mesh file in correct format and find lat/lon bounds
    checkMeshFile(args.meshfile)
    with xr.open_dataset('grd.tmp.nc') as ds0:
        # latitude range of grid
        bN=np.nanmax(ds0['nav_lat'])
        bS=np.nanmin(ds0['nav_lat'])

    if args.outfile is None:
        outFile='data_1m_VAR_nomask.nc'
    else:
        outFile=args.outfile
    vRep={'thetao':'potential_temperature','so':'salinity'}

    for iV,vV in enumerate(['thetao','so']):
        print(f'{vV}')
        startYear=min(args.years)
        if int(args.his2cmor)==1:
            file=his2cmor(vV,startYear,outFile,args)
        else:
            varPath=os.path.join(args.parent_path,args.parent_experiment,args.parent_ensemble,f'Omon/{vV}/gn/v20190429/')
            # find all files that match format in Omon
            flist=np.array(sorted(glob.glob(os.path.join(varPath,f'{vV}_Omon_{args.parent_name}_{args.parent_experiment}_{args.parent_ensemble}_gn_*.nc'))))

            # identify file based on desired year
            file,fileDates,_=matchFileYear(startYear,flist)
        if file is None:
            print(f"No files found.\n{vV} {startYear}")
        else:
            # extract desired year from file
            if len(file) > 1:
                # concatenate multiple files
                fstr=''
                for fle in file:
                    fstr+=f'{fle} '
                subprocess.run(f'ncrcat -h {fstr} -O {outFile}_y{fdates}_{vV}.nc',shell=True)
                subprocess.run(f'cdo --no_history selyear,{startYear}/{startYear} {outFile}_y{fdates}_{vV}.nc {outFile}_y{startYear}_{vV}.nc',shell=True)
                subprocess.run(f'rm -f {outFile}_y{fdates}_{vV}.nc',shell=True)
            else:
                subprocess.run(f'cdo --no_history selyear,{startYear}/{startYear} {file[0]} {outFile}_y{startYear}_{vV}.nc',shell=True)
            # get desired index/month
            subprocess.run(f'ncks -h -d time,{int(args.ic_ind)},{int(args.ic_ind)} {outFile}_y{startYear}_{vV}.nc -O {outFile}_y{startYear}{args.ic_ind+1:02}_{vV}.nc',shell=True)

            # subsample srcFile to be within +/- X degrees of southernmost point to speed things up
            subprocess.run(f'cdo --no_history sellonlatbox,-180,180,{bS-args.Xdeg},90 {outFile}_y{startYear}{args.ic_ind+1:02}_{vV}.nc {outFile}_y{startYear}{args.ic_ind+1:02}_{vV}.sliced.tmp.nc',shell=True)
            
            # fill any gaps in the sliced srcFile (two iterations)
            subprocess.run(f'cdo --no_history fillmiss2,2 {outFile}_y{startYear}{args.ic_ind+1:02}_{vV}.sliced.tmp.nc {outFile}_y{startYear}{args.ic_ind+1:02}_{vV}.filled.tmp.nc',shell=True)

            # remap to child grid
            subprocess.run(f"cdo --no_history remapdis,grd.tmp.nc {outFile}_y{startYear}{args.ic_ind+1:02}_{vV}.filled.tmp.nc {outFile}_y{startYear}{args.ic_ind+1:02}_{vV}.remapped.tmp.nc",shell=True)

            # interpolate vertically (if not SSH or sea ice)
            if vV in ['zos','siconc','sithick']:
                # interpolate filled src file to boundary points
                subprocess.run(f"mv {outFile}_y{startYear}{args.ic_ind+1:02}_{vV}.remapped.tmp.nc {outFile.replace('VAR',vRep[vV]).replace('yYYYY',f'y{startYear}')}.nc",shell=True)
            else:
                getZ(args.meshfile,outFile)
                # if not os.path.isfile(f'{outFile}.onlyz.tmp.nc'):
                #     with xr.open_dataset(args.meshfile) as mF:
                #         if 'e3t_0' in mF.keys():
                #             subprocess.run(f'ncks -h -v e3t_0 {args.meshfile} -O {outFile}.onlyz.tmp.nc',shell=True)
                #         elif 'tmask' in mF.keys():
                #             subprocess.run(f'ncks -h -v tmask {args.meshfile} -O {outFile}.onlyz.tmp.nc',shell=True)
                #         else:
                #             subprocess.run(f'ncks -h -v gdept_0 {args.meshfile} -O {outFile}.onlyz.tmp.nc',shell=True)
                subprocess.run(f"cdo --no_history intlevelx$(cdo -s showlevel {outFile}.onlyz.tmp.nc | tr ' ' ',') {outFile}_y{startYear}{args.ic_ind+1:02}_{vV}.remapped.tmp.nc {outFile.replace('VAR',vRep[vV]).replace('yYYYY',f'y{startYear}')}",shell=True)

            # remove intermediate variable files
            subprocess.run(f"rm -f {outFile.replace('VAR','*')}*.tmp.nc*",shell=True)
            subprocess.run(f'rm -f {outFile}*{vV}.nc',shell=True)
    # remove intermediate files
    subprocess.run(f'rm -f grd.tmp.nc',shell=True)

def bdy_slc(args):
    """
    Take global field and simply use a weight file. Only extract desired time and interpolate vertically.
    """
    print(f'Finding and processing boundaries from {args.parent_name}_{args.parent_experiment}_{args.parent_ensemble}.')

    if len(args.years) == 2 and (max(args.years)-min(args.years)) > 1:
        years=range(min(args.years),max(args.years)+1)
    else:
        years=np.array(args.years)

    # get output file name
    if args.outfile is None:
        # BDY gets replaced with north/south and YYYY gets replaced with the year
        outFile='obc_cantods025_yYYYY.nc'
    else:
        outFile=args.outfile

    # loop through each variable and interpolate to the regional boundaries
    vars2interp=['thetao','so','uo','vo','zos']
    for iV,vV in enumerate(vars2interp):
        varPath=os.path.join(args.parent_path,args.parent_experiment,args.parent_ensemble,f'Omon/{vV}/gn/v20190429/')

        # find all files that match format
        flist=np.array(sorted(glob.glob(os.path.join(varPath,f'{vV}_Omon_{args.parent_name}_{args.parent_experiment}_{args.parent_ensemble}_gn_*.nc'))))
        
        # identify and process file based on desired year
        # slow on first pass, but then much quicker on subsequent years if reading from same file
        file0=''
        for year in years:
            print(f'{vV} - {year}')
            file,fdates,sameFile=matchFileYear(year,flist,file0=file0)
            if (file is not None):
                # concatenate files if necessary
                if len(file) > 1:
                    fstr=''
                    for fle in file:
                        fstr+=f'{fle} '
                    subprocess.run(f'ncrcat -h {fstr} -O {outFile}.{vV}.concat.tmp.nc',shell=True)
                else:
                    subprocess.run(f'ln -s {file[0]} {outFile}.{vV}.concat.tmp.nc',shell=True)
                
                # subsample file to only include desired years
                subprocess.run(f"cdo --no_history selyear,{min(years)}/{max(years)} {outFile}.{vV}.concat.tmp.nc {outFile}.{vV}.sliced2.tmp.nc",shell=True)
                
                # fill any gaps in the sliced srcFile (two iterations)
                subprocess.run(f'cdo --no_history fillmiss2,2 {outFile}.{vV}.sliced2.tmp.nc {outFile}.{vV}.filled.tmp.nc',shell=True)

                # interpolate vertically (if not SSH)
                if vV == 'zos':
                    # interpolate filled src file to boundary points
                    subprocess.run(f"mv {outFile}.{vV}.filled.tmp.nc {outFile.replace('YYYY',f'{year}')}.{vV}.nc",shell=True)
                else:
                    getZ(args.meshfile,outFile)
                    subprocess.run(f"cdo --no_history intlevelx$(cdo -s showlevel {outFile}.onlyz.tmp.nc | tr ' ' ',') {outFile}.{vV}.filled.tmp.nc {outFile.replace('YYYY',f'{year}')}.{vV}.nc",shell=True)

                # remove intermediate files
                subprocess.run(f"rm -f {outFile}.*.tmp.nc*",shell=True)    

    # remove all remaining intermediate files generated above
    subprocess.run(f"rm -f {outFile}*.tmp.nc*",shell=True)

    # now concatenate physical variables and rename currents
    print(f'\rConcatenating files')
    for year in years:
        flist=''; ccount=0
        for vV in vars2interp:
            if os.path.isfile(f"{outFile.replace('YYYY',f'{year}')}.{vV}.nc"):
                flist+=f"{outFile.replace('YYYY',f'{year}')}.{vV}.nc "
                ccount+=1
        if ccount > 0:
            subprocess.run(f"cdo merge {flist} {outFile.replace('YYYY',f'{year}')}",shell=True)
            # remove individual variable files
            subprocess.run(f'rm -f {flist}',shell=True)
        
    # rename variables, dimensions, etc. for NEMO
    dimNames={'time':'t','lev':'z','i':'x','j':'y'}
    varNames={'longitude':'nav_lon','latitude':'nav_lat','i':'x','j':'y','thetao':'votemper',
            'so':'vosaline','uo':'vozocrtx','vo':'vomecrty','zos':'sossheig'}
    fcount=0
    for fF in sorted(glob.glob(f"{outFile.replace('YYYY','*')}")):
        # rename depth variable (shouldn't matter since all the same depth...)
        # and add attributes
        subprocess.run(f'cdo --no_history setattribute,thetao@grid=T,so@grid=T,zos@grid=T,uo@grid=U,vo@grid=V {fF} {fF}2',shell=True)
        subprocess.run(f'ncrename -h -v .lev,deptht {fF}2 -O {fF}',shell=True)
        subprocess.run(f'rm -f {fF}2',shell=True)

        # rename other variables and dimensions
        for iD,dD in enumerate(dimNames.keys()):
            subprocess.run(f'ncrename -h -d .{dD},{dimNames[dD]} {fF} -O {fF}',shell=True)
        for iV,vV in enumerate(varNames.keys()):
            subprocess.run(f'ncrename -h -v .{vV},{varNames[vV]} {fF} -O {fF}',shell=True)
        
        # count files
        fcount+=1

    # remove sliced meshfiles and intermediate files
    subprocess.run(f"rm -f {outFile.replace('YYYY','*')}.*.nc",shell=True)

    return fcount

def bdy_remap(args):
    """
    Remap to boundary.
    """
    print(f'Finding and processing boundaries from {args.parent_name}_{args.parent_experiment}_{args.parent_ensemble}.')

    if len(args.years) == 2 and (max(args.years)-min(args.years)) > 1:
        years=range(min(args.years),max(args.years)+1)
    else:
        years=np.array(args.years)

    # get output file name
    if args.outfile is None:
        # BDY gets replaced with north/south and YYYY gets replaced with the year
        outFile='obc_BDY_cantods025_yYYYY.nc'
    else:
        outFile=args.outfile

    # get boundaries from meshfile
    bN=-90.; bS=90.
    checkMeshFile(args.meshfile)
    # extract boundaries
    if args.bdyFile is not None:
        # boundary from defined file
        # identify boundary with file name
        bbase=os.path.basename(args.bdyFile).replace('.nc','')
        # get coordinates, include buffer zone, write to file
        with xr.open_dataset(args.bdyFile) as bdyFile:
            # write to temporary file
            msh=xr.Dataset.from_dict(
                {'nav_lon':{'dims':('y','x'),'data':bdyFile.nav_lon.values,'attrs':{'_CoordinateAxisType':'Lon','units':'degrees_east'}},
                'nav_lat':{'dims':('y','x'),'data':bdyFile.nav_lat.values,'attrs':{'_CoordinateAxisType':'Lat','units':'degrees_north'}},
                'nav_ones':{'dims':('time_counter','y','x'),'data':np.ones((1,np.shape(bdyFile.nav_lat.values)[0],np.shape(bdyFile.nav_lat.values)[1])),'attrs':{'coordinates':'nav_lat nav_lon'}},
                'time_counter':{'dims':('time_counter'),'data':[0.]}})
            # write to temporary netCDF file
            msh.to_netcdf(f'bdy.{bbase}.nc')
            blist=[bbase]
            # boundary max/min latitude
            bN=np.nanmax([bN,np.nanmax(bdyFile.nav_lat)])
            bS=np.nanmin([bS,np.nanmin(bdyFile.nav_lat)])
    else:
        # assumed boundary from grid
        blist=['north','south']
        for bdy in blist:
            # sub-sample meshfile to only include the boundary and 10 rows around it
            if bdy=='south':
                subprocess.run(f'ncks -h -d y,1,10 grd.tmp.nc -O bdy.south.nc',shell=True)
            elif bdy=='north':
                subprocess.run(f'ncks -h -d y,-11,-2 grd.tmp.nc -O bdy.north.nc',shell=True)

            # get the maximum latitude extent of the boundary to slice file (quicker processing)
            with xr.open_dataset(f'bdy.{bdy}.nc') as ds0:
                # boundary max/min latitude and longitude
                bN=np.nanmax([bN,np.nanmax(ds0['nav_lat'].values)])
                bS=np.nanmin([bS,np.nanmin(ds0['nav_lat'].values)])
        subprocess.run(f'rm -f grd.tmp.nc',shell=True)

    # loop through each variable and interpolate to the regional boundaries
    vars2interp=['thetao','so','uo','vo','zos']
    for iV,vV in enumerate(vars2interp):
         
        if int(args.his2cmor) != 1:
            varPath=os.path.join(args.parent_path,args.parent_experiment,args.parent_ensemble,f'Omon/{vV}/gn/v20190429/')

            # find all files that match format
            flist=np.array(sorted(glob.glob(os.path.join(varPath,f'{vV}_Omon_{args.parent_name}_{args.parent_experiment}_{args.parent_ensemble}_gn_*.nc'))))
        
        # identify and process file based on desired year
        # slow on first pass, but then much quicker on subsequent years if reading from same file
        file0=''
        if (args.iaf_year_offset is not None) and (args.iaf_loop_year is not None):
            # get max and min years in forcing range
            mnY=np.inf; mxY=-np.inf
            for year in years:
                fy=year + int(args.iaf_year_offset)
                yd=int(args.iaf_loop_year)-int(args.iaf_year_offset)
                yr=fy-int((year-1)/yd)*yd
                mnY=np.nanmin([mnY,yr])
                mxY=np.nanmax([mxY,yr])
        else:
            mnY=min(years); mxY=max(years)
        mnY=int(mnY)
        mxY=int(mxY)
        for year in years:
            print(f'{vV} - {year}')
            if (args.iaf_year_offset is not None) and (args.iaf_loop_year is not None):
                fy=year + int(args.iaf_year_offset)
                yd=int(args.iaf_loop_year)-int(args.iaf_year_offset)
                yr=fy-int((year-1)/yd)*yd
                if int(args.his2cmor)==1:
                    # if args.parent_name.lower()=='glorys':
                    #     file=glorys2cmor(vV,yr,outFile.replace('YYYY',f'{yr:04}'),args)
                    # else:
                    file=his2cmor(vV,yr,outFile.replace('YYYY',f'{yr:04}'),args)
                    fdates=f'{yr}'
                    sameFile=False
                else:
                    file,fdates,sameFile=matchFileYear(yr,flist,file0=file0)
            else:
                # find matching file
                if int(args.his2cmor)==1:
                    # if args.parent_name.lower()=='glorys':
                    #     file=glorys2cmor(vV,year,outFile.replace('YYYY',f'{year:04}'),args)
                    # else:
                    file=his2cmor(vV,year,outFile.replace('YYYY',f'{year:04}'),args)
                    fdates=f'{year}'
                    sameFile=False
                else:
                    file,fdates,sameFile=matchFileYear(year,flist,file0=file0)
                yr=int(year)

            if (file is not None) and (not sameFile):
                file0=file
                            
                # subsample srcFile to be within +/- X degrees of the boundary
                if len(file) > 1:
                    fstr=''
                    for fle in file:
                        fstr+=f'{fle} '
                    subprocess.run(f'ncrcat -h {fstr} -O {outFile}.{vV}.concat.tmp.nc',shell=True)
                else:
                    subprocess.run(f'ln -s {file[0]} {outFile}.{vV}.concat.tmp.nc',shell=True)
                with xr.open_dataset(file[0]) as src:
                    # find indices of region north and south of the boundary
                    try:
                        ilat=np.where(np.sum(np.logical_and(src['latitude']>=bS-args.Xdeg,src['latitude']<=bN+args.Xdeg),axis=1))[0]
                        latLen=len(src['latitude'].isel(i=0))
                    except:
                        ilat=np.where(np.logical_and(src['latitude']>=bS-args.Xdeg,src['latitude']<=bN+args.Xdeg))[0]
                        latLen=len(src['latitude'])
                # add some extra latitude points if min and max the same, or only one point
                if (np.nanmax(ilat)-np.nanmin(ilat)) <= 1:
                    ilat=[np.nanmax([0,np.nanmin(ilat)-1]),np.nanmin([np.nanmax(ilat)+1,latLen])]
                subprocess.run(f'ncks -h -d j,{np.nanmin(ilat)}.,{np.nanmax(ilat)}. {outFile}.{vV}.concat.tmp.nc -O {outFile}.{vV}.sliced.tmp.nc',shell=True)

                # subsample file to only include desired years
                subprocess.run(f"cdo --no_history selyear,{mnY}/{mxY} {outFile}.{vV}.sliced.tmp.nc {outFile}.{vV}.sliced2.tmp.nc",shell=True)
                
                # fill any gaps in the sliced srcFile (two iterations)
                subprocess.run(f'cdo --no_history fillmiss2,2 {outFile}.{vV}.sliced2.tmp.nc {outFile}.{vV}.filled.tmp.nc',shell=True)

                # interpolate vertically (if not SSH)
                if vV == 'zos':
                    # interpolate filled src file to boundary points
                    subprocess.run(f'mv {outFile}.{vV}.filled.tmp.nc {outFile}.{vV}.z.tmp.nc',shell=True)
                else:
                    getZ(args.meshfile,outFile)
                    subprocess.run(f"cdo --no_history intlevelx$(cdo -s showlevel {outFile}.onlyz.tmp.nc | tr ' ' ',') {outFile}.{vV}.filled.tmp.nc {outFile}.{vV}.z.tmp.nc",shell=True)

                # remap to meshfiles to boundaries
                for bdy in blist:
                    subprocess.run(f"cdo --no_history remapdis,bdy.{bdy}.nc {outFile}.{vV}.z.tmp.nc {outFile.replace('BDY',bdy)}_y{fdates}.{vV}.tmp.nc",shell=True)

                # remove intermediate variable files
                subprocess.run(f'rm -f {outFile}.{vV}*.tmp.nc*',shell=True)
            
            if file is not None:
                # slice files to match the desired date range
                for bdy in blist:
                    subprocess.run(f"cdo --no_history selyear,{yr}/{yr} {outFile.replace('BDY',bdy)}_y{fdates}.{vV}.tmp.nc {outFile.replace('BDY',bdy).replace('YYYY',f'{year}')}.{vV}.nc",shell=True)                    
            else:
                print('No files found.')
        
        # remove intermediate files
        subprocess.run(f"rm -f {outFile.replace('BDY','*')}_y*-*.{vV}*",shell=True)
    
    # remove all remaining intermediate files generated above
    subprocess.run(f"rm -f {outFile.replace('YYYY','*')}*.tmp.nc*",shell=True)

    # now concatenate physical variables and rename currents
    print(f'\rConcatenating files')
    for bdy in blist:
        outFileBdy=f"{outFile.replace('BDY',bdy)}"
        for year in years:
            flist=''; ccount=0
            for vV in vars2interp:
                if os.path.isfile(f"{outFileBdy.replace('YYYY',f'{year}')}.{vV}.nc"):
                    print(f"Adding: {outFileBdy.replace('YYYY',f'{year}')}.{vV}.nc")
                    flist+=f"{outFileBdy.replace('YYYY',f'{year}')}.{vV}.nc "
                    ccount+=1
                else:
                    print(f"{outFileBdy.replace('YYYY',f'{year}')}.{vV}.nc does not exist")
            if ccount > 0:
                print(f"cdo merge {flist} {outFileBdy.replace('YYYY',f'{year}')}")
                subprocess.run(f"cdo merge {flist} {outFileBdy.replace('YYYY',f'{year}')}",shell=True)
                # remove individual variable files
                subprocess.run(f'rm -f {flist}',shell=True)
            
        # rename variables, dimensions, etc. for NEMO
        dimNames={'time':'t','lev':'z','i':'x','j':'y'}
        varNames={'longitude':'nav_lon','latitude':'nav_lat','i':'x','j':'y','thetao':'votemper',
                'so':'vosaline','uo':'vozocrtx','vo':'vomecrty','zos':'sossheig'}
        fcount=0
        for fF in sorted(glob.glob(f"{outFileBdy.replace('YYYY','*')}")):
            # rename depth variable (shouldn't matter since all the same depth...)
            # and add attributes
            subprocess.run(f'cdo --no_history setattribute,thetao@grid=T,so@grid=T,zos@grid=T,uo@grid=U,vo@grid=V {fF} {fF}2',shell=True)
            subprocess.run(f'ncrename -h -v .lev,deptht {fF}2 -O {fF}',shell=True)
            subprocess.run(f'rm -f {fF}2',shell=True)

            # rename other variables and dimensions
            for iD,dD in enumerate(dimNames.keys()):
                subprocess.run(f'ncrename -h -d .{dD},{dimNames[dD]} {fF} -O {fF}',shell=True)
            for iV,vV in enumerate(varNames.keys()):
                subprocess.run(f'ncrename -h -v .{vV},{varNames[vV]} {fF} -O {fF}',shell=True)
            
            # count files
            fcount+=1

    # remove sliced meshfiles and intermediate files
    subprocess.run(f"rm -f {outFile.replace('BDY',bdy).replace('YYYY',f'{year}')}.{vV}.nc",shell=True)
    subprocess.run(f"rm -f {outFile.replace('BDY',bdy)}*.tmp.nc",shell=True)
    for bdy in blist:
        subprocess.run(f'rm -f bdy.{bdy}.nc',shell=True)

    return fcount

def frc_slice(args):
    """
    Slice atmospheric forcing into annual files (if not already)
    """
    # variables to process
    if args.forcing == 'OMIP':
        vars2process=['ncar_precip','ncar_rad','q_10','slp','t_10','u_10','v_10']
        vRep={'ncar_precip':'prec','ncar_rad':'rad','q_10':'humi','slp':'slp',
              't_10':'tair','u_10':'u','v_10':'v'}
    else:
        vars2process=['tas','uas','vas','huss','rsds','rlds','pr','prsn','psl','ps']
        vRep={'tas':'tair','uas':'u','vas':'v','huss':'humi','rsds':'qsr',
            'rlds':'qlw','pr':'prec','prsn':'snow','psl':'slp','ps':'slp'}

    if args.outfile is None:
        outFile=f'VAR_cantods025_yYYYY'
    else:
        outFile=args.outfile

    if len(args.years)==2 and (max(args.years)-min(args.years)) > 1:
        years=range(min(args.years),max(args.years)+1)
    else:
        years=np.array(args.years)

    if args.forcing=='OMIP':
        print(f'Finding OMIP forcing.')
    else:
        print(f'Finding forcing from {args.parent_name}_{args.parent_experiment}_{args.parent_ensemble}.')
    fexpect=0; fcount=0
    for iV,vV in enumerate(vars2process):
        if args.forcing=='OMIP':
            varPath='/space/hall6/sitestore/eccc/crd/ccrn/users/rdy001/forcing/corev2-ciaf'
            # find all files that match format
            flist=glob.glob(os.path.join(varPath,f'{vV}.1948-2009.23OCT2012.nc'))
        else:
            varPath=os.path.join(args.parent_path,args.parent_experiment,args.parent_ensemble,f'{args.mor}/{vV}/gn/v20190429/')
            # find all files that match format
            flist=np.array(sorted(glob.glob(os.path.join(varPath,f'{vV}_{args.mor}_{args.parent_name}_{args.parent_experiment}_{args.parent_ensemble}_gn_*.nc'))))

        # identify file based on desired year
        for year in years:
            fexpect+=1
            if (args.iaf_year_offset is not None) and (args.iaf_loop_year is not None):
                fy=year + int(args.iaf_year_offset)
                yd=int(args.iaf_loop_year)-int(args.iaf_year_offset)
                yr=fy-int((year-1)/yd)*yd
                file=matchFileYear(yr,flist)[0]
            else:
                # find matching file
                file=matchFileYear(year,flist)[0]
                yr=year
            print(f'{vV} - {year} ({yr})')

            if file is None:
                print(f'No {vV} files found for {year} ({yr}).')
            else:
                # get desired time from file
                if len(file) > 1:
                    fstr=''
                    for fle in file:
                        fstr+=f'{fle} '
                    subprocess.run(f'ncrcat -h {fstr} -O {outFile}.{vV}.concat.tmp.nc',shell=True)
                else:
                    subprocess.run(f'ln -s {file[0]} {outFile}.{vV}.concat.tmp.nc',shell=True)
                subprocess.run(f"cdo --no_history selyear,{yr}/{yr} {outFile}.{vV}.concat.tmp.nc {outFile.replace('VAR',vRep[vV]).replace('yYYYY',f'y{year:04}')}.nc",shell=True)
                subprocess.run(f"rm -f {outFile}.{vV}.concat.tmp.nc",shell=True)
                # fill missing points (e.g., in raw OMIP files)
                if args.forcing=='OMIP':
                    # need to rename lat/lon
                    subprocess.run(f"ncrename -h -v LON,nav_lon -v LAT,nav_lat {outFile.replace('VAR',vRep[vV]).replace('yYYYY',f'y{year:04}')}.nc -O {outFile.replace('VAR',vRep[vV]).replace('yYYYY',f'y{year:04}')}.nc",shell=True)
                    # subprocess.run(f"ncatted -h -a _CoordinateAxisType,nav_lon,o,c,Lon -a _CoordinateAxisType,nav_lat,o,c,Lat {outFile.replace('VAR',vRep[vV]).replace('yYYYY',f'y{year:04}')}.nc -O {outFile.replace('VAR',vRep[vV]).replace('yYYYY',f'y{year:04}')}.nc",shell=True)
                # subprocess.run(f"cdo --no_history fillmiss2 {outFile.replace('VAR',vRep[vV]).replace('yYYYY',f'y{year:04}')}.tmp.nc {outFile.replace('VAR',vRep[vV]).replace('yYYYY',f'y{year:04}')}.nc",shell=True)
                # subprocess.run(f"rm -f {outFile.replace('VAR',vRep[vV]).replace('yYYYY',f'y{year:04}')}.tmp*.nc",shell=True)
                # else:
                #     # subprocess.run(f"ncatted -a _FillValue,{vV},d,, {outFile.replace('VAR',vRep[vV]).replace('yYYYY',f'y{year:04}')}.nc -O {outFile.replace('VAR',vRep[vV]).replace('yYYYY',f'y{year:04}')}-2.nc",shell=True)
                #     # subprocess.run(f"ncatted -h -a coordinates,{vV},o,c,'lat lon' {outFile.replace('VAR',vRep[vV]).replace('yYYYY',f'y{year:04}')}.nc -O {outFile.replace('VAR',vRep[vV]).replace('yYYYY',f'y{year:04}')}-2.nc",shell=True)
                #     subprocess.run(f"cdo --no_history fillmiss2 {outFile.replace('VAR',vRep[vV]).replace('yYYYY',f'y{year:04}')}-2.nc {outFile.replace('VAR',vRep[vV]).replace('yYYYY',f'y{year:04}')}-3.nc",shell=True)
                #     subprocess.run(f"\mv {outFile.replace('VAR',vRep[vV]).replace('yYYYY',f'y{year:04}')}-3.nc {outFile.replace('VAR',vRep[vV]).replace('yYYYY',f'y{year:04}')}.nc",shell=True)
                #     subprocess.run(f"ncrename -h -d lon,x -d lat,y -v lon,nav_lon -v lat,nav_lat {outFile.replace('VAR',vRep[vV]).replace('yYYYY',f'y{year:04}')}.nc -O {outFile.replace('VAR',vRep[vV]).replace('yYYYY',f'y{year:04}')}.nc",shell=True)
                    # subprocess.run(f"rm -f {outFile.replace('VAR',vRep[vV]).replace('yYYYY',f'y{year:04}')}-*.nc",shell=True)
                # if OMIP precipitation, need to have snow and rain+snow variables
                    if vV == 'ncar_precip':
                        subprocess.run(f"ncap2 -O -s 'PRECIP=RAIN+SNOW' {outFile.replace('VAR',vRep[vV]).replace('yYYYY',f'y{year:04}')}.nc {outFile.replace('VAR',vRep[vV]).replace('yYYYY',f'y{year:04}')}.nc",shell=True)
                fcount+=1
    return fcount,fexpect

def rvr_remap(args):
    """
    (Very) simple river remapping (simply scales one file to match another).
    """
    #TODO: perform area-weighted sum of river discharges
    # calculate total disharge in kg/s
    # calculate ratio of total discharge
    # scale ratio of total discharge by ratio of water areas
    # scale discharge in second file by those values

    if args.outfile is None:
        outFile=f'rvr_cantods025'
    else:
        outFile=args.outfile
    
    if len(args.years)==2 and (max(args.years)-min(args.years)) > 1:
        years=range(min(args.years),max(args.years)+1)
    else:
        years=np.array(args.years)

    print(f'Finding river inputs forcing from {args.parent_name}_{args.parent_experiment}_{args.parent_ensemble}.')
    # find all files that match format
    rvrPath=os.path.join(args.parent_path,args.parent_experiment,args.parent_ensemble,f'Omon/friver/gn/v20190429/')
    # find all files that match format
    flist=np.array(sorted(glob.glob(os.path.join(rvrPath,f'friver_Omon_{args.parent_name}_{args.parent_experiment}_{args.parent_ensemble}_gn_*.nc'))))
    if len(flist)==0:
        print('No river files!')
    else:
        for year in years:
            if args.iaf_year_offset is not None and args.iaf_loop_year is not None:
                fy=year + int(args.iaf_year_offset)
                yd=int(args.iaf_loop_year)-int(args.iaf_year_offset)
                yr=fy-int((year-1)/yd)*yd
            else:
                yr=year
            file0=matchFileYear(yr,flist)[0]
            if file0 is None:
                print(f'No river file found for {year} ({yr}).')
            else:
                # concatenate multiple files
                if len(file0) > 1:
                    fstr=''
                    for fle in file0:
                        fstr+=f'{fle} '
                    subprocess.run(f'ncrcat -h {fstr} -O {outFile}.concat.tmp.nc',shell=True)
                else:
                    subprocess.run(f'ln -s {file0[0]} {outFile}.concat.tmp.nc',shell=True)
                # if only one file given, simply slice for correct dates
                if (args.rvr is None):
                    print(f'Slicing river file - {year} ({yr})')
                    subprocess.run(f"cdo --no_history selyear,{yr}/{yr} {outFile}.concat.tmp.nc {outFile.replace('yYYYY',f'y{year:04}')}.nc",shell=True)
                # otherwise, remap
                else:
                    print(f'Remapping {file0} - {year} ({yr})') 
                    # get river scaling and apply to new file
                    # sum up all values in file 0 and file 1
                    if 'flist1' not in locals():
                        flist1=np.array(sorted(glob.glob(args.rvr)))
                    if len(flist1)==0:
                        print(args.rvr)
                        print('No files for remapping rivers!')
                    else:
                        file1=matchFileYear(yr,flist1)[0]
                        if file1 is None:
                            print('No exact file found for remapping. Instead, taking first file found.')
                            file1=[flist1[0]]
                        if len(file1) > 1:
                            fstr=''
                            for fle in file1:
                                fstr+=f'{fle} '
                            subprocess.run(f'ncrcat -h {fstr} -O {outFile}.concat1.tmp.nc',shell=True)
                        else:
                            subprocess.run(f'ln -s {file1[0]} {outFile}.concat1.tmp.nc',shell=True)
                        # get time-sliced files
                        subprocess.run(f"cdo --no_history selyear,{yr}/{yr} {outFile}.concat.tmp.nc {outFile.replace('yYYYY',f'y{year:04}')}.tmp0.nc",shell=True)
                        # get single year from target file (may not match yr)
                        syr=os.path.basename(args.rvr).split('_')[-1].split('-')[0][0:4]
                        subprocess.run(f"cdo --no_history selyear,{syr}/{syr} {outFile}.concat1.tmp.nc {outFile.replace('yYYYY',f'y{year:04}')}.tmp1.nc",shell=True)
                        area0=cellAreas(args.parent_grid)
                        area1=cellAreas(args.meshfile)
                        with xr.open_dataset(f"{outFile.replace('yYYYY',f'y{year:04}')}.tmp0.nc") as rvr0:
                            with xr.open_dataset(f"{outFile.replace('yYYYY',f'y{year:04}')}.tmp1.nc") as rvr1:
                                sum0=np.nansum(rvr0.friver*area0)#; ar0=np.nansum(area0[rvr0.friver > 0])
                                sum1=np.nansum(rvr1.friver*area1)#; ar1=np.nansum(area1[rvr1.friver > 0])
                                # get ratio
                                ratio=sum0/sum1#; aratio=ar0/ar1
                                # check: does sum of rvr1*ratio == sum0??
                                diffCheck=100*np.abs(np.nansum(ratio*rvr1.friver*area1)-sum0)/sum0
                                if diffCheck > 0.1:
                                    sys.exit(f'Total fresh water does not agree!\n{sum0}\n{np.nansum(ratio*rvr1.friver*area1)}\n{diffCheck}%')
                                else:
                                    print(f'Total freshwater within 0.1% ({diffCheck}%)')
                        # apply ratio to files
                        subprocess.run(f"cdo expr,'friver={ratio}*friver' {outFile.replace('yYYYY',f'y{year:04}')}.tmp1.nc {outFile.replace('yYYYY',f'y{year:04}')}.nc",shell=True)
                    subprocess.run(f"rm -f {outFile.replace('yYYYY',f'*')}.tmp*.nc",shell=True)
                subprocess.run(f'rm -f {outFile}*concat*tmp.nc',shell=True)
                
                # add river mask to file
                with xr.open_dataset(f"{outFile.replace('yYYYY',f'y{year:04}')}.nc") as rvr:
                    rflag=np.squeeze(np.nansum(rvr.friver,axis=0))>0
                    rmask=np.full(np.shape(rflag),0.0)
                    rmask[rflag]=0.5
                rset=xr.Dataset.from_dict({'riv_mask':{'dims':('y','x'),'data':rmask}})
                
                # write to temporary netCDF file
                rset.to_netcdf(f"{outFile.replace('yYYYY',f'y{year:04}')}.rmask.nc")
                # # scale friver by 0
                # subprocess.run(f"cdo --no_history expr,\"friver=0.0*friver\" {outFile.replace('yYYYY',f'y{year:04}')}.nc {outFile.replace('yYYYY',f'y{year:04}')}.0scale.nc",shell=True)
                # add to original file
                subprocess.run(f"cdo --no_history merge {outFile.replace('yYYYY',f'y{year:04}')}.nc {outFile.replace('yYYYY',f'y{year:04}')}.rmask.nc {outFile.replace('yYYYY',f'y{year:04}')}.masked.nc",shell=True)
                # delete temporary file
                subprocess.run(f"rm -f {outFile.replace('yYYYY',f'y{year:04}')}.rmask.nc",shell=True)
    return

def his2cmor(var4cmor,fYear,fOut,args):
    """
    Naive conversion from history files to "CMORized" file to use with scripts above
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
              subprocess.run(f"cdo --no_history sellonlatbox,-180,180,30,45 {hisFile} {cmfile}",shell=True)
              subprocess.run(f"ncks -v {var4cmor} {cmfile} -O {cmfile}",shell=True)
        else:
            subprocess.run(f"ncks -v {var4cmor} {hisFile} -O {cmfile}",shell=True)
        # rename coordinates (if applicable)
        subprocess.run(f"ncrename -h -d .x,i -d .y,j -d .longitude,i -d .latitude,j -v .nav_lon,longitude -v .nav_lat,latitude {cmfile} -O {cmfile}",shell=True)
        # subprocess.run(f"ncrename -h -v .{lev[ftype[var4cmor]]},lev -d .{lev[ftype[var4cmor]]},lev {cmfile} -O {cmfile}",shell=True)
        # rename vertical coordinate (if applicable)
        # either a grid-specific name, or just generic 'depth'
        subprocess.run(f"ncrename -h -d .{lev[ftype[var4cmor]]},lev {cmfile} -O {cmfile}",shell=True)
        subprocess.run(f"ncrename -h -d .depth,lev {cmfile} -O {cmfile}",shell=True)
        # if var4cmor != 'zos':
        #     subprocess.run(f"ncks -h -x -v {lev[ftype[var4cmor]]} -d {lev[ftype[var4cmor]]}_bnds {cmfile} -O {cmfile}",shell=True)
        subprocess.run(f"ncatted -h -a coordinates,{var4cmor},o,c,'latitude longitude' {cmfile} -O {cmfile}",shell=True)
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
            subprocess.run(f"ncatted -h -a positive,lev,o,c,down {cmfile}.ij.tmp.nc -O {cmfile}.ij.tmp.nc",shell=True)
            subprocess.run(f"ncatted -h -a axis,lev,o,c,Z {cmfile}.ij.tmp.nc -O {cmfile}.ij.tmp.nc",shell=True)
        if reLat:
            # want lat/lon to be gridded, not vectors
            # remove lon & lat from file
            subprocess.run(f"ncks -v {var4cmor} {cmfile}.ij.tmp.nc {cmfile}.nolonlat.tmp.nc",shell=True)
            lnlt=xr.Dataset.from_dict(
                    {'longitude':{'dims':('j','i'),'data':LN,'attrs':{'_CoordinateAxisType':'Lon','units':'degrees_east'}},
                     'latitude':{'dims':('j','i'),'data':LT,'attrs':{'_CoordinateAxisType':'Lat','units':'degrees_north'}}})
            lnlt.to_netcdf(f'{cmfile}.lonlat.tmp.nc')
            # with xr.open_dataset(f"{cmfile}.nolonlat.tmp.nc") as cxlnlt:
            #     dstmp=xr.merge([cxlnlt,lnlt])
            # dstmp.to_netcdf(f'{cmfile}.lonlat.tmp.nc')
            subprocess.run(f"cdo merge {cmfile}.nolonlat.tmp.nc {cmfile}.lonlat.tmp.nc {cmfile}.gridded.tmp.nc",shell=True)
            subprocess.run(f"ncatted -h -a coordinates,{var4cmor},o,c,'latitude longitude' {cmfile}.gridded.tmp.nc -O {cmfile}.gridded.tmp.nc",shell=True)
            cmfile=[f'{cmfile}.gridded.tmp.nc']
        else:
            cmfile=[f'{cmfile}.ij.tmp.nc']

    return cmfile

# def glorys2cmor(var4cmor,fYear,fOut,args):
#     """
#     Naive conversion from GLORYS files to "CMORized" file to use with scripts above
#     """
#     ftype={'thetao':'grid_[tT]','so':'grid_[tT]','zos':'grid_[tT]','uo':'grid_[uU]','vo':'grid_[vV]'}
#     lev={'grid_[tT]':'deptht','grid_[uU]':'depthu','grid_[vV]':'depthv'}

#     # find correct file
#     pPath=args.parent_path.split('nc_output')[0]
#     varPath=os.path.join(pPath,f"mc_{args.parent_name.replace('CanESM5-NEMO4-','')}_{fYear}_*{ftype[var4cmor]}.nc*")
    
#     # find all files that match format
#     flist=np.array(sorted(glob.glob(varPath)))

#     if len(flist) == 0:
#         cmfile=None
#     else:
#         # get highest value of file version
#         hisFile=flist[-1]
#         cmfile=f"{fOut.replace('VAR',var4cmor)}.{var4cmor}.2cmor.tmp.nc"

#         # open glorys file for reading
#         with xr.open_dataset(hisFile) as hF:
#             # gridded coordinates
#             ln=hF.longitude.values
#             lt=hF.latitude.values
#             LN,LT=np.meshgrid(ln,lt)
#             shpi=len(ln)
#             shpj=len(lt)
#             # open new file for writing
#             with Dataset(cmfile,'w') as cF:
#                 # create dimensions
#                 cF.createDimension('time',None)
#                 if var4cmor != 'zos':
#                     cF.createDimension('lev',len(hF.depth.values))
#                 cF.createDimension('j',len(hF.latitude.values))
#                 cF.createDimension('i',len(hF.longitude.values))
#                 # create variables
#                 time=cF.createVariable('time','f8',('time',)); time.units="hours since 1950-01-01"; time.calendar="gregorian"
#                 if var4cmor!='zos':
#                     lev=cF.createVariable('lev','f4',('lev',))
#                     var=cF.createVariable(var4cmor,'f8',('time','lev','j','i')); var.coordinates="latitude longitude"
#                 else:
#                     var=cF.createVariable(var4cmor,'f8',('time','j','i')); var.coordinates="latitude longitude"
#                 lat=cF.createVariable('latitude','f4',('j','i')); lat._CoordinateAxisType='Lat'; lat.units='degrees_north'
#                 lon=cF.createVariable('longitude','f4',('j','i')); lon._CoordinateAxisType='Lon'; lon.units='degrees_east'
#                 ii=cF.createVariable('i','i8',('i',)); ii.units='1'; ii.long_name='first spatial index for variables stored on an unstructured grid'
#                 jj=cF.createVariable('j','i8',('j',)); jj.units='1'; jj.long_name='second spatial index for variables stored on an unstructured grid'
                
#                 # write variables
#                 time[:]=hF.time.values
#                 if var4cmor!='zos':
#                     lev[:]=hF.depth.values
#                 lat[:]=LT
#                 lon[:]=LN
#                 ii[:]=range(shpi)
#                 jj[:]=range(shpj)
#                 var[:]=hF[var4cmor].values
#             print(f'Wrote temporary GLORYS file {var4cmor}')
#         cmfile=[cmfile]
       
        # # slice file to get desired variable
        # if 'glorys' in cmfile:
        #     # restrict latitude further as files are BIG
        #     subprocess.run(f"cdo --no_history sellonlatbox,-180,180,30,45 {hisFile} {cmfile}",shell=True)
        #     subprocess.run(f"ncks -v {var4cmor} {cmfile} -O {cmfile}",shell=True)
        # else:
        #     subprocess.run(f"ncks -v {var4cmor} {hisFile} -O {cmfile}",shell=True)
        # # rename coordinates (if applicable)
        # subprocess.run(f"ncrename -h -d .x,i -d .y,j -d .longitude,i -d .latitude,j -v .nav_lon,longitude -v .nav_lat,latitude {cmfile} -O {cmfile}",shell=True)
        # # subprocess.run(f"ncrename -h -v .{lev[ftype[var4cmor]]},lev -d .{lev[ftype[var4cmor]]},lev {cmfile} -O {cmfile}",shell=True)
        # # rename vertical coordinate (if applicable)
        # # either a grid-specific name, or just generic 'depth'
        # subprocess.run(f"ncrename -h -d .{lev[ftype[var4cmor]]},lev {cmfile} -O {cmfile}",shell=True)
        # subprocess.run(f"ncrename -h -d .depth,lev {cmfile} -O {cmfile}",shell=True)
        # # if var4cmor != 'zos':
        # #     subprocess.run(f"ncks -h -x -v {lev[ftype[var4cmor]]} -d {lev[ftype[var4cmor]]}_bnds {cmfile} -O {cmfile}",shell=True)
        # subprocess.run(f"ncatted -h -a coordinates,{var4cmor},o,c,'latitude longitude' {cmfile} -O {cmfile}",shell=True)
        # # add gridded lat/lon if necessary
        # with xr.open_dataset(cmfile) as cmtmp:
        #     shp=np.shape(cmtmp.latitude.values)
        #     try:
        #         shpi=shp[1]
        #         shpj=shp[0]
        #         reLat=False
        #     except:
        #         ln=cmtmp.longitude.values
        #         lt=cmtmp.latitude.values
        #         LN,LT=np.meshgrid(ln,lt)
        #         shpi=len(ln)
        #         shpj=len(lt)
        #         reLat=True
        # #add dimension variables to file         
        # with xr.open_dataset(cmfile) as cmtmp:
        #     if var4cmor=='zos':
        #         ij=xr.Dataset.from_dict(
        #             {'i':{'dims':('i'),'data':range(shpi),'attrs':{'units':'1','long_name':'first spatial index for variables stored on an unstructured grid'}},
        #             'j':{'dims':('j'),'data':range(shpj),'attrs':{'units':'1','long_name':'second spatial index for variables stored on an unstructured grid'}}})
        #     else:
        #         with xr.open_dataset(hisFile) as hstmp:
        #             try:
        #                 depth=hstmp.variables[lev[ftype[var4cmor]]]
        #             except:
        #                 try:
        #                     depth=hstmp.variables['lev']
        #                 except:
        #                     depth=hstmp.variables['depth']
        #             ij=xr.Dataset.from_dict(
        #                 {'i':{'dims':('i'),'data':range(shpi),'attrs':{'units':'1','long_name':'first spatial index for variables stored on an unstructured grid'}},
        #                 'j':{'dims':('j'),'data':range(shpj),'attrs':{'units':'1','long_name':'second spatial index for variables stored on an unstructured grid'}},
        #                 'lev':{'dims':('lev'),'data':depth,'attrs':{'long_name':"Vertical T levels",'units':'m','positive':"down",'axis':"Z"}}})
        #     cmtmp2=xr.merge([cmtmp,ij])
        #     cmtmp2.to_netcdf(f'{cmfile}.ij.tmp.nc')
        #     subprocess.run(f"ncatted -h -a positive,lev,o,c,down {cmfile}.ij.tmp.nc -O {cmfile}.ij.tmp.nc",shell=True)
        #     subprocess.run(f"ncatted -h -a axis,lev,o,c,Z {cmfile}.ij.tmp.nc -O {cmfile}.ij.tmp.nc",shell=True)
        # if reLat:
        #     # want lat/lon to be gridded, not vectors
        #     # remove lon & lat from file
        #     subprocess.run(f"ncks -v {var4cmor} {cmfile}.ij.tmp.nc {cmfile}.nolonlat.tmp.nc",shell=True)
        #     lnlt=xr.Dataset.from_dict(
        #             {'longitude':{'dims':('j','i'),'data':LN,'attrs':{'_CoordinateAxisType':'Lon','units':'degrees_east'}},
        #              'latitude':{'dims':('j','i'),'data':LT,'attrs':{'_CoordinateAxisType':'Lat','units':'degrees_north'}}})
        #     lnlt.to_netcdf(f'{cmfile}.lonlat.tmp.nc')
        #     # with xr.open_dataset(f"{cmfile}.nolonlat.tmp.nc") as cxlnlt:
        #     #     dstmp=xr.merge([cxlnlt,lnlt])
        #     # dstmp.to_netcdf(f'{cmfile}.lonlat.tmp.nc')
        #     subprocess.run(f"cdo merge {cmfile}.nolonlat.tmp.nc {cmfile}.lonlat.tmp.nc {cmfile}.gridded.tmp.nc",shell=True)
        #     print("I AM HERE!!")
        #     subprocess.run(f"ncatted -h -a coordinates,{var4cmor},o,c,'latitude longitude' {cmfile}.gridded.tmp.nc -O {cmfile}.gridded.tmp.nc",shell=True)
        #     cmfile=[f'{cmfile}.gridded.tmp.nc']
        # else:
        #     cmfile=[f'{cmfile}.ij.tmp.nc']

    # return cmfile

def fixVelocities(vFile,outVFile):
    """
    Flip meridional velocities in the "northern" portion of a domain (if applicable)
    TODO: CanTODS CMORization will require correction step?!
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


def calc_nemo_chunk_dates(args):
    """
    Compute the start and end date of the current chunk
    Copied and modified from: CanNEMO/lib/maestro/bin/compute_iteration_info.py
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

#---------------#
# END FUNCTIONS #
#---------------#

# parse arguments
args=parser.parse_args()

start=time.time()

# get necessary info about run
if args.type==-1:
    calc_nemo_chunk_dates(args)
elif args.type==0:
    if args.hot_start==0:
        # restart from rest using history/CMORized files
        ic_remap(args)
    else:
        # hot start
        rs_remap(args)
    print(f'Done making IC files!\n{time.time()-start} s elapsed.')
elif args.type==1:
    fcount=bdy_remap(args)
    print(f'Done making the boundary files!\n{fcount} files created.\n{time.time()-start} s elapsed.')
elif args.type==2:
    fcount,fexpect=frc_slice(args)
    print(f'Done finding forcing files!\n({fcount}/{fexpect} found)\n{time.time()-start} s elapsed.')
elif args.type==3:
    print('Remapping/scaling rivers.')
    rvr_remap(args)
    print(f'Done river rivermapping.')
elif args.type==10:
    bdy_slc(args)
    print(f'Done extracting file as boundary input!')
else:
    sys.exit(f'ERROR: type must be -1, 0, 1, or 2, not {args.type}')
