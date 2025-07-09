"""
This function remaps river output from a parent grid to
a child grid. It will need refining!!

For now, it takes the CMORized river output and either
1) slices it in time to match the desired run (grids are the same)
or 2) scales an existing river file to match the total freshwater
input of the parent run, keeping the ratio of river discharge in
the template file.

Written by: Jonathan Izett (2024/2025)
"""

from argparse import Namespace
import glob
import numpy as np
import os
import sys
import xarray as xr

from remap_functions.remap_utils import cellAreas, cleanTmp, matchFileYear, parseArgYears, selYear, subproc

def rvr_remap(args: Namespace):
    """
    (Very) simple river remapping (simply scales one file to match another).
    WARNING: THIS IS A VERY CRUDE APPROXIMATION!!!

    Inputs:
        args - arguments passed to remap_canesm.py
    Outputs:
        Writes a river forcing file to NEMO.

    """
    # calculate total disharge in kg/s
    # calculate ratio of total discharge
    # scale ratio of total discharge by ratio of water areas
    # scale discharge in second file by those values

    if args.outfile is None:
        outFile=f'rvr_cantods025'
    else:
        outFile=args.outfile
    
    years=parseArgYears(args)

    print(f'Finding river inputs forcing from {args.parent_name}_{args.parent_experiment}_{args.parent_ensemble}.')
    # find all files that match format
    rvrPath=os.path.join(args.parent_path,args.parent_experiment,args.parent_ensemble,f'Omon/friver/gn/v20190429/')
    # find all files that match format
    fPattern=os.path.join(rvrPath,f'friver_Omon_{args.parent_name}_{args.parent_experiment}_{args.parent_ensemble}_gn_*.nc')

    for year in years:
        if args.iaf_year_offset is not None and args.iaf_loop_year is not None:
            fy=year + int(args.iaf_year_offset)
            yd=int(args.iaf_loop_year)-int(args.iaf_year_offset)
            yr=fy-int((year-1)/yd)*yd
        else:
            yr=year
        file0=matchFileYear(yr,fPattern)[0]
        if file0 is None:
            print(f'No river file found for {year} ({yr}).')
        else:
            fr1=True    # flag if file contains friver (instead of runoff)
            # concatenate multiple files
            if len(file0) > 1:
                fstr=''
                for fle in file0:
                    fstr+=f'{fle} '
                err=subproc(f"ncrcat -h {fstr} -O {outFile}.concat.tmp.nc")
            else:
                err=subproc(f"ln -sf {file0[0]} {outFile}.concat.tmp.nc")
            # if only one file given, simply slice for correct dates
            if (args.rvr is None):
                print(f'Slicing river file - {year} ({yr})')
                selYear(yr, yr, f"{outFile}.concat.tmp.nc", f"{outFile.replace('yYYYY',f'y{year:04}')}.nc")
            # otherwise, remap
            else:
                print(f'Remapping {file0[0]} - {year} ({yr})') 
                # get river scaling and apply to new file
                # sum up all values in file 0 and file 1
                if 'flist1' not in locals():
                    flist1=np.array(sorted(glob.glob(args.rvr)))
                if len(flist1)==0:
                    print(args.rvr)
                    print('No files for remapping rivers!')
                else:
                    file1=matchFileYear(yr,args.rvr)[0]
                    exact=True
                    if file1 is None:
                        print('No exact file found for remapping. Instead, taking first file found.')
                        file1=[flist1[0]]
                        exact=False
                    if len(file1) > 1:
                        fstr=''
                        for fle in file1:
                            fstr+=f'{fle} '
                        err=subproc(f"ncrcat -h {fstr} -O {outFile}.concat1.tmp.nc")
                    else:
                        err=subproc(f"ln -sf {file1[0]} {outFile}.concat1.tmp.nc")
                    # get time-sliced files
                    selYear(yr, yr, f"{outFile}.concat.tmp.nc", f"{outFile.replace('yYYYY',f'y{year:04}')}.tmp0.nc")
                    # get single year from target file (may not match yr)
                    syr=os.path.basename(args.rvr).split('_')[-1].split('-')[0][0:4]
                    try:
                        selYear(syr, syr, f"{outFile}.concat1.tmp.nc", f"{outFile.replace('yYYYY',f'y{year:04}')}.tmp1.nc")
                    except:
                        err=subproc(f"mv {outFile}.concat1.tmp.nc {outFile.replace('yYYYY',f'y{year:04}')}.tmp1.nc")
                    if not os.path.isfile(f"{outFile.replace('yYYYY',f'y{year:04}')}.tmp1.nc"):
                        err=subproc(f"ln -sf {outFile}.concat1.tmp.nc {outFile.replace('yYYYY',f'y{year:04}')}.tmp1.nc")
                    area0=cellAreas(args.parent_grid)
                    area1=cellAreas(args.meshfile)
                    with xr.open_dataset(f"{outFile.replace('yYYYY',f'y{year:04}')}.tmp0.nc") as rvr0:
                        with xr.open_dataset(f"{outFile.replace('yYYYY',f'y{year:04}')}.tmp1.nc") as rvr1:
                            rr0=rvr0.friver
                            try:
                                rr1=rvr1.friver
                                fr1=True
                            except:
                                try:
                                    rr1=rvr1.runoff[:,0:-1,1:-1]
                                except:
                                    rr1=rvr1.runoff[0:-1,1:-1]
                                fr1=False
                            if len(rr0)==12 and len(rr1)==12:
                                # both files are monthly, can simply sum
                                sum0=np.nansum(rr0*area0)
                                sum1=np.nansum(rr1*area1)
                            else:
                                # files are not the same time span, need to average
                                # in time before taking the sum to get total volume
                                if len(rr0)==12:
                                    sum0=np.nansum(np.nanmean(rr0,axis=0)*area0)
                                else:
                                    sum0=np.nansum(rr0*area0)
                                if len(rr1)==12:
                                    sum1=np.nansum(np.nanmean(rr1,axis=0)*area1)
                                else:
                                    sum1=np.nansum(rr1*area1)
                            # get ratio
                            ratio=sum0/sum1
                    # apply ratio to files
                    if fr1:
                        # err=subroc(f"cdo expr,\"friver={ratio}\*friver\" {outFile.replace('yYYYY',f'y{year:04}')}.tmp1.nc {outFile.replace('yYYYY',f'y{year:04}')}.nc")
                        err=subproc(f"ncap2 -s friver={ratio}*friver {outFile.replace('yYYYY',f'y{year:04}')}.tmp1.nc -O {outFile.replace('yYYYY',f'y{year:04}')}.nc")
                    else:
                        # err=subproc(f"cdo expr,\"runoff={ratio}\*runoff\" {outFile.replace('yYYYY',f'y{year:04}')}.tmp1.nc {outFile.replace('yYYYY',f'y{year:04}')}.nc")
                        err=subproc(f"ncap2 -s runoff={ratio}*runoff {outFile.replace('yYYYY',f'y{year:04}')}.tmp1.nc -O {outFile.replace('yYYYY',f'y{year:04}')}.nc")
                    # if 'friver' in rvr0.keys():
                    #             rr0=rvr0.friver
                    #         else:
                    #             rr0=rvr0.runoff
                    #         if len(np.shape(rr0))==3:
                    #             sum0=np.nansum(np.nanmean(rr0,axis=0)*area0)
                    #         else:
                    #             sum0=np.nansum(rr0*area0)
                    #         if 'friver' in rvr1.keys():
                    #             rr1=rvr1.friver
                    #             fr1=True
                    #         else:
                    #             rr1=rvr1.runoff
                    #             fr1=False
                    #         if len(np.shape(rr1))==3:
                    #             sum1=np.nansum(np.nanmean(rr1,axis=0)*area1)
                    #         else:
                    #             sum0=np.nansum(rr1*area1)
                    #         # get ratio
                    #         ratio=sum0/sum1
                    # # apply ratio to files
                    # if fr1:
                    #     subprocess.run(f"cdo expr,'friver={ratio}*friver' {outFile.replace('yYYYY',f'y{year:04}')}.tmp1.nc {outFile.replace('yYYYY',f'y{year:04}')}.nc",shell=True)
                    # else:
                    #     subprocess.run(f"cdo expr,'runoff={ratio}*runoff' {outFile.replace('yYYYY',f'y{year:04}')}.tmp1.nc {outFile.replace('yYYYY',f'y{year:04}')}.nc",shell=True)
                cleanTmp(f"{outFile.replace('yYYYY',f'*')}.tmp*.nc")
            cleanTmp(f"rm -f {outFile}*concat*tmp.nc")
            
            # add river mask to file
            with xr.open_dataset(f"{outFile.replace('yYYYY',f'y{year:04}')}.nc") as rvr:
                if fr1:
                    if len(rvr.friver)==12:
                        rflag=np.squeeze(np.nansum(rvr.friver,axis=0))>0    
                    else:
                        rflag=np.squeeze(rvr.friver)>0    
                else:
                    if len(rvr.runoff)==12:
                        rflag=np.squeeze(np.nansum(rvr.runoff,axis=0))>0
                    else:
                        rflag=np.squeeze(rvr.runoff)>0
                rmask=np.full(np.shape(rflag),0.0)
                rmask[rflag]=0.5
            rset=xr.Dataset.from_dict({'riv_mask':{'dims':('y','x'),'data':rmask}})
            
            # write to temporary netCDF file
            rset.to_netcdf(f"{outFile.replace('yYYYY',f'y{year:04}')}.rmask.nc")
            # add to original file
            err=subproc(f"cdo --no_history merge {outFile.replace('yYYYY',f'y{year:04}')}.nc {outFile.replace('yYYYY',f'y{year:04}')}.rmask.nc {outFile.replace('yYYYY',f'y{year:04}')}.masked.nc")
            # delete temporary file
            cleanTmp(f"{outFile.replace('yYYYY',f'y{year:04}')}.rmask.nc")
    return