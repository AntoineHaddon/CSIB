"""
These functions take surface forcing files (atmospheric and salinity restoring)
and slice/remap the files for use in NEMO.

    Functions:
        frc_slice - gets desired time range from atmospheric forcing
        sos_remap - takes surface salinity and remaps it to the child grid

Written by: Jonathan Izett (2024/2025)

"""

from argparse import Namespace
import glob
import numpy as np
import os
import sys
from typing import Tuple
import xarray as xr

# utilities for remapping
from remap_functions.remap_utils import checkMeshFile, cleanTmp, fillMiss2, matchFileYear, parseArgYears, remap, selBox, selYear, subproc

def frc_slice(args: Namespace) -> Tuple[int, int]:
    """
    Extracts the desired time range from atmospheric forcing file(s)
    to be used as input to NEMO with associated weighting files
    (generated separately).

    Inputs:
        args - arguments to remap_canesm.py
    Outputs:
        Writes forcing files on original grid with appropriate time slices.
        Returns number of files processed successfully and number of files expected.
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

    years=parseArgYears(args)

    if args.forcing=='OMIP':
        print(f'Finding OMIP forcing.')
    else:
        print(f'Finding forcing from {args.parent_name}_{args.parent_experiment}_{args.parent_ensemble}.')
    fexpect=0; fcount=0
    for iV,vV in enumerate(vars2process):
        if args.forcing=='OMIP':
            varPath='/space/hall6/sitestore/eccc/crd/ccrn/users/rdy001/forcing/corev2-ciaf'
            # find all files that match format
            fPattern=os.path.join(varPath,f'{vV}.1948-2009.23OCT2012.nc')
        else:
            varPath=os.path.join(args.parent_path,args.parent_experiment,args.parent_ensemble,f'{args.mor}/{vV}/gn/v20190429/')
            # find all files that match format
            fPattern=os.path.join(varPath,f'{vV}_{args.mor}_{args.parent_name}_{args.parent_experiment}_{args.parent_ensemble}_gn_*.nc')

        # identify file based on desired year
        for year in years:
            fexpect+=1
            if (args.iaf_year_offset is not None) and (args.iaf_loop_year is not None):
                fy=year + int(args.iaf_year_offset)
                yd=int(args.iaf_loop_year)-int(args.iaf_year_offset)
                yr=fy-int((year-1)/yd)*yd
                file=matchFileYear(yr,fPattern)[0]
            else:
                # find matching file
                file=matchFileYear(year,fPattern)[0]
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
                    err=subproc(f"ncrcat -h {fstr.strip(' ')} -O {outFile}.{vV}.concat.tmp.nc")
                else:
                    err=subproc(f"ln -sf {file[0]} {outFile}.{vV}.concat.tmp.nc")
                selYear(yr, yr, f"{outFile}.{vV}.concat.tmp.nc", f"{outFile.replace('VAR',vRep[vV]).replace('yYYYY',f'y{year:04}')}.nc")
                cleanTmp(f"{outFile}.{vV}.concat.tmp.nc")
                # fill missing points (e.g., in raw OMIP files)
                if args.forcing=='OMIP':
                    # need to rename lat/lon
                    err=subproc(f"ncrename -h -v LON,nav_lon -v LAT,nav_lat {outFile.replace('VAR',vRep[vV]).replace('yYYYY',f'y{year:04}')}.nc -O {outFile.replace('VAR',vRep[vV]).replace('yYYYY',f'y{year:04}')}.nc")
                # if OMIP precipitation, need to have snow and rain+snow variables
                    if vV == 'ncar_precip':
                        err=subproc(f"ncap2 -O -s 'PRECIP=RAIN+SNOW' {outFile.replace('VAR',vRep[vV]).replace('yYYYY',f'y{year:04}')}.nc {outFile.replace('VAR',vRep[vV]).replace('yYYYY',f'y{year:04}')}.nc")
                fcount+=1
    return fcount,fexpect

def frc_remap(args: Namespace) -> Tuple[int, int]:
    """
    Remap atmospheric forcing from CanESM grid to desired grid.

    Inputs:
        args - args passed to remap_canesm.py
    Output:
        Writes atmospheric forcing files for NEMO.
    """
    print(f'Finding and processing data from {args.parent_name}_{args.parent_experiment}_{args.parent_ensemble}.')

    # get mesh file in correct format and find lat/lon bounds
    bN,bS=checkMeshFile(args.meshfile)

    if args.outfile is None:
        outFile=f'VAR_atm_yYYYY'
    else:
        outFile=args.outfile

    years=parseArgYears(args)

    # variables to remap
    if args.parent_name=='omip':
      fVars=['u_10','v_10','ncar_rad','t_10','q_10','ncar_precip_mod','slp']
    else:
      fVars=['uas','vas','pr','prsn','rlds','rsds','tas','ps']

    # loop through desired years
    fcount=0
    for year in years:
      # loop through each variable 
      for vV in fVars:
          if args.parent_name=='omip':
            omip=True
            if '025' in args.parent_experiment: 
                file=f"{vV}_omip025_y{year:04}.nc"
                acFile=f"urdy_eorca025_corev2_{vV}_nomask_y{year:04}.nc"
            else:
                file=f"{vV}_omip1_y{year:04}.nc"
                acFile=f"urdy_eorca1_corev2_{vV}_nomask_js_y{year:04}.nc"
            subproc(f"access {file} {acFile}")
            file=[file]
            startYear=year
          else:
            omip=False
            varPath=os.path.join(args.parent_path,args.parent_experiment,args.parent_ensemble,f'{args.mor}/{vV}/gn/v20190429/')
            # find all files that match format
            fPattern=os.path.join(varPath,f'{vV}_{args.mor}_{args.parent_name}_{args.parent_experiment}_{args.parent_ensemble}_gn_*.nc')

            # identify file based on desired year
            if (args.iaf_year_offset is not None) and (args.iaf_loop_year is not None):
                fy=year + int(args.iaf_year_offset)
                yd=int(args.iaf_loop_year)-int(args.iaf_year_offset)
                yr=fy-int((year-1)/yd)*yd
                file=matchFileYear(yr,fPattern)[0]
            else:
                # find matching file
                file=matchFileYear(year,fPattern)[0]
                yr=year
            fdates=f'{yr}-{year}'
            startYear=yr

          if file is None:
              print(f'No {vV} files found for {year} ({yr}).')
          else:
              # get desired time from file
              # extract desired year from file
              if len(file) > 1:
                  # concatenate multiple files
                  fstr=''
                  for fle in file:
                      fstr+=f'{fle} '
                      print(fle)
                  print("")
                  err=subproc(f"ncrcat -h {fstr}-O {outFile}_y{fdates}_{vV}.nc")
                  selYear(startYear, startYear, f"{outFile}_y{fdates}_{vV}.nc", f"{outFile}_y{year:04}_{vV}.nc")
                  cleanTmp(f"{outFile}_y{fdates}_{vV}.nc")
              else:
                  selYear(startYear, startYear, file[0], f"{outFile}_y{year:04}_{vV}.nc")
          
              # subsample srcFile to be within +/- X degrees of southernmost point to speed things up
              selBox(-180, 180, bS-args.Xdeg, 90, f"{outFile}_y{year:04}_{vV}.nc", f"{outFile}_y{year:04}_{vV}.sliced.tmp.nc")
          
              # fill any gaps in the sliced srcFile (two iterations)
              fillMiss2(f"{outFile}_y{year:04}_{vV}.sliced.tmp.nc", f"{outFile}_y{year:04}_{vV}.filled.tmp.nc")
  
              # remap to child grid
              remap("grd.tmp.nc", f"{outFile}_y{year:04}_{vV}.filled.tmp.nc", f"{outFile.replace('VAR',vV).replace('yYYYY',f'y{year:04}')}.nc")
              # remove temporary files
              cleanTmp(f"{outFile}_y{year:04}_{vV}.*nc")
              fcount+=1
              
              if omip:
                subproc(f"rm -f {file[0]}")

    # remove intermediate files
    cleanTmp(f"*.tmp.nc")

    return fcount, len(years)*len(fVars)

def sos_remap(args: Namespace) -> None:
    """
    Remap surface salinity for sssr.

    Inputs:
        args - args passed to remap_canesm.py
    Output:
        Writes sssr file for NEMO.
    """
    print(f'Finding and processing data from {args.parent_name}_{args.parent_experiment}_{args.parent_ensemble}.')

    # get mesh file in correct format and find lat/lon bounds
    bN,bS=checkMeshFile(args.meshfile)

    if args.outfile is None:
        outFile=f'sss_data_yYYYY.nc'
    else:
        outFile=args.outfile

    years=parseArgYears(args)

    # hard code for now...in future make it flexible
    vV='so'
    varPath=os.path.join(args.parent_path,args.parent_experiment,args.parent_ensemble,f'Omon/{vV}/gn/v20190429/')
    # find all files that match format
    fPattern=os.path.join(varPath,f'{vV}_Omon_{args.parent_name}_{args.parent_experiment}_{args.parent_ensemble}_gn_*.nc')

    # identify file based on desired year
    for year in years:
        if (args.iaf_year_offset is not None) and (args.iaf_loop_year is not None):
            fy=year + int(args.iaf_year_offset)
            yd=int(args.iaf_loop_year)-int(args.iaf_year_offset)
            yr=fy-int((year-1)/yd)*yd
            file=matchFileYear(yr,fPattern)[0]
        else:
            # find matching file
            file=matchFileYear(year,fPattern)[0]
            yr=year
        fdates=f'{yr}-{year}'
        startYear=yr

        if file is None:
            print(f'No {vV} files found for {year} ({yr}).')
        else:
            # get desired time from file
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
        
            # subsample srcFile to be within +/- X degrees of southernmost point to speed things up
            selBox(-180, 180, bS-args.Xdeg, 90, f"{outFile}_y{startYear}_{vV}.nc", f"{outFile}_y{startYear}{args.ic_ind+1:02}_{vV}.sliced.tmp.nc")
        
            # only keep surface value
            err=subproc(f"ncks -h -d lev,0,0 {outFile}_y{startYear}{args.ic_ind+1:02}_{vV}.sliced.tmp.nc -O {outFile}_y{startYear}{args.ic_ind+1:02}_{vV}.sliced.tmp.nc")

            # fill any gaps in the sliced srcFile (two iterations)
            fillMiss2(f"{outFile}_y{startYear}{args.ic_ind+1:02}_{vV}.sliced.tmp.nc", f"{outFile}_y{startYear}{args.ic_ind+1:02}_{vV}.filled.tmp.nc")

            # remap to child grid
            remap("grd.tmp.nc", f"{outFile}_y{startYear}{args.ic_ind+1:02}_{vV}.filled.tmp.nc", f"{outFile.replace('VAR',vV).replace('yYYYY',f'y{year:04}')}")
            # remove temporary files
            cleanTmp(f"{outFile}*_y*_{vV}*nc")

    # remove intermediate files
    cleanTmp(f"*.tmp.nc")

