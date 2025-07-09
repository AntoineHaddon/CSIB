"""
These functions allow remapping of initial conditions (either cold-start files, or NEMO restart files)

  ic_remap : remaps cold-start IC files
  rs_remap : remaps NEMO restart files

Written by: Jonathan Izett (2024/2025)
"""

from argparse import Namespace
import glob
import numpy as np
import os
import sys
import xarray as xr

# utilities for remapping
from remap_functions.remap_utils import addCoords, changeCoords, checkMeshFile, cleanTmp, fillMiss2, getZ, intZ, matchFileYear, remap, selBox, selYear, subproc

def ic_remap(args: Namespace):
    """
    Remap initial condition for cold start.

    Input:
        args - arguments passed to remap_canesm.py
    Output:
        Writes initial condition T&S files to args.outfile
    """

    # get month index as integer
    args.ic_ind=int(args.ic_ind)

    # loop through each variable
    print(f'Finding and processing data from {args.parent_name}_{args.parent_experiment}_{args.parent_ensemble}.')

    # get mesh file in correct format and find lat/lon bounds
    bN,bS=checkMeshFile(args.meshfile)

    # set output file name if not already
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
            # identify file based on desired year
            file,fileDates,_=matchFileYear(startYear,os.path.join(varPath,f'{vV}_Omon_{args.parent_name}_{args.parent_experiment}_{args.parent_ensemble}_gn_*.nc'))

        if file is not None:
            # extract desired year from file
            if len(file) > 1:
                # concatenate multiple files
                fstr=''
                for fle in file:
                    fstr+=f'{fle} '
                err=subproc(f"ncrcat -h {fstr} -O {outFile}_y{fdates}_{vV}.nc")
                selYear(startYear, startYear, f"{outFile}_y{fdates}_{vV}.nc", f"{outFile}_y{startYear}_{vV}.nc")
                cleanTmp(f"{outFile}_y{fdates}_{vV}.nc")
            else:
                selYear(startYear, startYear, file[0], f"{outFile}_y{startYear}_{vV}.nc")
            # get desired index/month
            err=subproc(f"ncks -h -d time,{int(args.ic_ind)},{int(args.ic_ind)} {outFile}_y{startYear}_{vV}.nc -O {outFile}_y{startYear}{args.ic_ind+1:02}_{vV}.nc")

            # subsample srcFile to be within +/- X degrees of southernmost point to speed things up
            selBox(-180, 180, np.nanmax([-90,bS-args.Xdeg]), np.nanmin([90,bN+args.Xdeg]), f"{outFile}_y{startYear}{args.ic_ind+1:02}_{vV}.nc", f"{outFile}_y{startYear}{args.ic_ind+1:02}_{vV}.sliced.tmp.nc")
            
            # fill any gaps in the sliced srcFile (two iterations)
            fillMiss2(f"{outFile}_y{startYear}{args.ic_ind+1:02}_{vV}.sliced.tmp.nc", f"{outFile}_y{startYear}{args.ic_ind+1:02}_{vV}.filled.tmp.nc")

            # remap to child grid
            remap("grd.tmp.nc", f"{outFile}_y{startYear}{args.ic_ind+1:02}_{vV}.filled.tmp.nc", f"{outFile}_y{startYear}{args.ic_ind+1:02}_{vV}.remapped.tmp.nc")

            # interpolate vertically (if not SSH or sea ice)
            if vV in ['zos','siconc','sithick']:
                # single level, so no interpolation required
                err=subproc(f"mv {outFile}_y{startYear}{args.ic_ind+1:02}_{vV}.remapped.tmp.nc {outFile.replace('VAR',vRep[vV]).replace('yYYYY',f'y{startYear}')}.nc")
            else:
                # interpolate vertically
                getZ(args.meshfile,outFile)
                intZ(f"{outFile}_y{startYear}{args.ic_ind+1:02}_{vV}.remapped.tmp.nc", f"{outFile.replace('VAR',vRep[vV]).replace('yYYYY',f'y{startYear}')}", f"{outFile}.onlyz.tmp.nc")

            # remove intermediate variable files
            cleanTmp(f"{outFile.replace('VAR','*')}*.tmp.nc*")
            cleanTmp(f"{outFile}*{vV}.nc")
    # # remove intermediate files
    cleanTmp(f"*.tmp.*nc*")

def rs_remap(args: Namespace):
    """
    Remap restart files for hot start.

    Input:
        args - arguments passed to remap_canesm.py
    Output:
        Writes restart[_ice,_trc] files to args.outfile
    """

    # loop through each variable
    print(f'Finding and processing data from {args.parent_name}_{args.parent_experiment}_{args.parent_ensemble}.')

    # get mesh file in correct format and find lat/lon bounds
    bN,bS=checkMeshFile(args.meshfile)

    # set output file name if not already
    if args.outfile is None:
        outFile='TYPE'
    else:
        outFile=args.outfile

    # find directory containing restart files
    filePath=np.array(sorted(glob.glob(os.path.join(args.parent_path,f'mc_{args.parent_name}_{min(args.years):04}_m*_nemors.*'))))
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
        rsFile=glob.glob(os.path.join(filePath[0],f"{args.parent_name}*_{tT}.nc"))[0]
        changeCoords(rsFile,f"{outFile.replace('TYPE',tT)}.atts.tmp.nc")

        # loop through variables in file and assign coordinates as appropriate
        geoVars=[]; scalars=[]
        with xr.open_dataset(f"{outFile.replace('TYPE',tT)}.atts.tmp.nc") as attFle:
            for vV in attFle.variables.keys():
                if vV not in ['nav_lon','nav_lat'] and len(np.shape(attFle.variables[vV]))>2:
                    # keep track of variables with geographic coordinates to modify metadata
                    geoVars.append(vV)
                elif len(np.shape(attFle.variables[vV]))==0:
                    # keep track of scalars to add back to file at end
                    scalars.append(vV)
        for gV in geoVars:
            # add coordinates to each variable
            addCoords(gV,f"{outFile.replace('TYPE',tT)}.atts.tmp.nc")

        # subsample srcFile to be within +/- X degrees of southernmost point to speed things up
        selBox(-180,180,np.nanmax([-90,bS-args.Xdeg]),np.nanmin([90,bN+args.Xdeg]),f"{outFile.replace('TYPE',tT)}.atts.tmp.nc",f"{outFile.replace('TYPE',tT)}.sliced.tmp.nc")
        
        # fill any gaps in the sliced srcFile (two iterations)
        fillMiss2(f"{outFile.replace('TYPE',tT)}.sliced.tmp.nc", f"{outFile.replace('TYPE',tT)}.filled.tmp.nc")

        # remap to child grid
        remap("grd.tmp.nc", f"{outFile.replace('TYPE',tT)}.filled.tmp.nc", f"{outFile.replace('TYPE',tT)}.remapped.tmp.nc")

        # interpolate vertically if appropriate
        if args.outfile is None and tT != 'restart':
            tSuff='_in'
        else:
            tSuff=''
        # TODO: Interpolate vertical levels...doesn't seem to work properly?! For now, leave out since grids all eORCA
        # with same vertical resolution.
        # if tT == 'restart_ice':
        #     # No vertical interpolation for ice files
        #     err=subproc(f"mv {outFile.replace('TYPE',tT)}.remapped.tmp.nc {outFile.replace('TYPE',tT)}.tmp.zmap.nc")
        # else:
        #     # interpolate vertical levels
        #     getZ(args.meshfile,outFile)
        #     err=subproc(f"cdo -w --no_history intlevelx$(cdo -s showlevel {outFile}.onlyz.tmp.nc | tr ' ' ',') {outFile.replace('TYPE',tT)}.remapped.tmp.nc {outFile.replace('TYPE',tT)}.tmp.zmap.nc")

        # add scalars back to file (they get dropped in cdo commands)
        err=subproc(f"cp {outFile.replace('TYPE',tT)}.remapped.tmp.nc {outFile.replace('TYPE',tT)}{tSuff}.nc")
        for sS in scalars:
            err=subproc(f"ncks -h -v {sS} {outFile.replace('TYPE',tT)}.atts.tmp.nc -A {outFile.replace('TYPE',tT)}{tSuff}.nc")

        if tT == 'restart' and int(args.flip_vel)==1:
            fixVelocities(f"{outFile.replace('TYPE',tT)}{tSuff}.nc",f"{outFile.replace('TYPE',tT)}{tSuff}.nc.velflip.nc")
            err=subproc(f"rm -f {outFile.replace('TYPE',tT)}{tSuff}.nc ; mv {outFile.replace('TYPE',tT)}{tSuff}.nc.velflip.nc {outFile.replace('TYPE',tT)}{tSuff}.nc")

        # remove intermediate variable files
        cleanTmp(f"{outFile.replace('TYPE',tT)}*.tmp.nc*")
    # remove remaining intermediate files
    cleanTmp("*.tmp.*nc*")

    return