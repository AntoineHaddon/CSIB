"""
These utility functions are called by cantods_rtd (and submodules)
to facilitate calculation of diagnostics specific to CanTODS.
"""
#TODO: CLEAN UP
#TODO: write in format/link to RTD browser?? (would need to update RTD browser's variables...)
from argparse import Namespace
import datetime as dt
import gc
import glob
import matplotlib.pyplot as plt
from matplotlib import path
import numpy as np
import os
import pandas as pd
import subprocess
import time
import warnings
import xarray as xr

def splitCMD(cmd: str) -> List[str]:
    """
    Split a command, preserving values in quotations.

    Input:
        cmd - string of desired command
    Output:
        Returns string split into a list, preserving quotation marks
    """
    if "\"" in cmd or "(" in cmd:
        # remove leading/trailing spaces
        cmd=cmd.strip(" ")
        # split, preserving brackets and quotation marks
        nb_brackets=0; nb_quotes=0
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
            elif cC==" " and nb_brackets==0 and nb_quotes==0:
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

def subscript(cmd: str) -> None:
    """
    Create a temporary shell script and run it.

    Input:
        cmd - string containing command
    """

    with open("tmp.sh","w") as shFile:
        shFile.write(f"#!/bin/bash\n\nsource ~/.profile\n\n{cmd}")
    err=subproc("chmod +x tmp.sh")
    err=subproc("tmp.sh")          # run temporary shell script
    err=subproc("rm -f tmp.sh")    # remove temporary shell script

    return

def nemoFileType(var: str) -> str:
    """
    This function takes a nemo variable and
    determines what type of file it is in.
        diaptr
        icemod
        grid_T
        grid_U
        grid_V
        grid_W
    Note: many variables still need to be added!!
    """

    if var in ['thetao','so']:
        filetype='grid_T'
    elif var in ['siconc','sithick']:
        filetype='icemod'
    elif var in ['uo']:
        filetype='grid_U'
    elif var in ['vo']:
        filetype='grid_V'
    elif var in ['wo']:
        filetype='grid_W'
    else:
        sys.exit(f'Unsure where to find {var}!')

    return filetype

def cantodsFiles(runid: str, nemoVar: str, args: Namespace) -> Tuple[List[str], Namespace]:
    """
    Get list of files that fall within the desired time range.

    Inputs:
        runid - name of run
        nemoVar - variable of interest
        args - args passed to cantods_rtd.py
    Outputs:
        flist - list of files that match search criteria
        args - updated args
    """

    # NEMO file type
    ftype=nemoFileType(nemoVar)

    # search for all files that match
    if len(list(args.runpath)) == 1:
        globdir=os.path.join(args.runpath[0],runid,'data')
    else:
        globdir=os.path.join(np.array(args.runpath)[np.array(args.runid)==runid][0],runid,'data')
    # flist=glob.glob(f'{globdir}/mc_{runid}_*_m01_5d_{ftype.lower()}.nc.*')
    flist=glob.glob(f'{globdir}/mc_{runid}_*_m01_1m_{ftype.lower()}.nc.*')

    # identify the year in each file name
    fyear=[int(os.path.basename(fle).split(f'{runid}_')[1].split('_m')[0]) for fle in flist]

    # get start year if not passed
    if args.year0 is None:
        # assume it starts on 1 January
        args.year0=np.nanmin(fyear)
    
    # only keep files that are within the desired range
    fyear=np.array(fyear)-args.year0
    fGood=np.full(len(fyear),False,dtype='bool')
    fGood[np.logical_and(fyear>=args.xlim[0],fyear<=args.xlim[1])]=True

    # only keep files with desired file version
    if args.version is None:
        # loop through files to find latest file version
        ver=0
        for fF in np.array(flist)[fGood]:
            ver=np.nanmax([ver,int(fF.split('.nc.')[-1])])
        args.version=f'{ver:03}'
    for iF,fF in enumerate(flist):
        if fGood[iF]:
            if f'.nc.{args.version}' in os.path.basename(fF):
                fGood[iF]=True
            else:
                fGood[iF]=False
    
    # sort files by date
    flist=sorted(list(np.array(flist)[fGood]))

    return flist,args

def relativeYear(dset: xr.Dataset, args: Namespace) -> np.ndarray:
    """
    Convert from date in file to years since start date

    Inputs:
        dset - xarray Dataset containing model output
        args - arguments passed to cantods_rtd.py
    Output:
        fyears - years since start of model run
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore",category=RuntimeWarning)
        try:
            ftime=np.array(dset.indexes['time_counter'].values)
        except:
            ftime=np.array(dset.indexes['time_counter'].to_datetimeindex().values)
    
    d0=dt.datetime(args.year0,1,1)

    #TODO: handle leap years properly?
    fyears=np.full(len(ftime),np.nan)
    for iY,fT in enumerate(ftime):
        try:
          fyears[iY]=fT.year-args.year0 + (fT.month-1)/12 + (fT.day-1)/365. + (fT.hour-1)/(24*365.) + (fT.minute-1)/(60*24*365.) + (fT.second-1)/(60*60*24*365.)    # fraction by seconds in that minute
        except:
          # depending on the python version, numpy date object may not have access to year, etc. Instead, need to convert to string and parse that
          fS=str(fT)
          fY=int(fS.split('-')[0])
          fM=int(fS.split('-')[1])
          fD=int(fS.split('-')[2].split('T')[0])
          fH=int(fS.split('T')[1].split(':')[0])
          fm=int(fS.split(':')[1])
          fs=float(fS.split(':')[2])
          fyears[iY]=fY-args.year0 + (fM-1.0)/12.0 + (fD-1.0)/365. + (fH-1.0)/(24*365.) + (fm-1.0)/(60*24*365.) + (fs-1.0)/(60*60*24*365.)    # fraction by seconds in that minute

    return fyears

def meshGrid(meshfiles: Optional[List[str]], runid: str, olev: int, args: Namespace) -> Dict[str, np.ndarray]:
    """
    Get the grid cell areas from the meshfile.

    Inputs:
        meshfiles - list of files containing grid info.
        runid - run name
        olev - desired level
        args - arguments passed to cantods_rtd.py
    Outputs:
        Dictionary containing grid cell masks, area, volume, and depth
    """

    if meshfiles is None:
        # locate appropriate meshfile for current grid
        # search for all files that match
        if len(list(args.runpath)) == 1:
            globdir=os.path.join(args.runpath[0],runid,'data')
        else:
            globdir=os.path.join(np.array(args.runpath)[np.array(args.runid)==runid][0],runid,'data')
        meshfiles=glob.glob(f'{globdir}/mc_{runid}_*_m01_mesh_mask.nc.*')

    # get grid info from meshfile; loop until one works
    for meshfile in meshfiles:
        try:
            with xr.open_dataset(meshfile) as meshmask:
                if olev==-1:
                    lev=np.shape(meshmask['tmask'])[1]+1
                else:
                    lev=olev+1

                e1t = meshmask['e1t'].values.squeeze() # y,x
                e2t = meshmask['e2t'].values.squeeze() # y,x
                try:
                    e3t = meshmask['e3t_0'].values.squeeze()[0,:,:,:] # z,y,x
                except:
                    e3t = meshmask['e3t_0'].values.squeeze() # z,y,x

                # land mask
                tmask = np.squeeze(np.array(meshmask['tmask'].values[0,:,:,:].squeeze()[0:lev,:,:]))

                # calulate grid area
                if lev==1:
                    gridArea = e1t*e2t*tmask
                else:
                    gridArea = e1t*e2t*tmask[0,:,:]
                
                # grid cell depths (horizontally uniform)
                depths=np.cumsum(e3t,axis=0)[:,0,0]

                # get grid volume
                e3t = np.array(e3t[0:lev,:,:])*tmask
                gridVol=np.tile(gridArea,(len(e3t),1,1))*e3t
                if lev==1:
                    hdep = e3t.copy()
                else:
                    hdep = np.sum(e3t, axis = 0)
                hdep[hdep==0]=np.nan

            return {'tmask':tmask,'gridArea':gridArea,'gridVol':gridVol,'hdep':hdep,'e3t':e3t,'depths':depths}
        except:
            pass

def regionMask(region: str, lon: np.ndarray, lat: np.ndarray) -> Optional[np.ndarray]:
    """
    Get a 2D mask defining a given region.
    Regions include:
        labsea - Labrador Sea
        canpac - Canadian Pacific Ocean
        canatl - Canadian Atlantic Ocean
        canarc - Canadian Arctic Ocean
        hbay   - Hudson Bay
        arctic - Arctic circle (>~ 66.5 N)
    """
    # Bounds for some geometrically specified regions
    if region.lower()=='labsea' or ('lab' in region.lower()):
        # LABRADOR SEA
                # Kulliq    Greenland NFLD    Kippon     Belle Isl. NE LEdge  Cape St Charles
        bounds=[[-64.5786, -43.91753, -52.45, -55.41667, -55.4293, -55.25,   -55.6215, -64.5786], #longitude
                [ 60.,      60.,       47.75,  51.6667,   51.9,     52.0333,  52.21667, 60.]]     # latitude
    elif region.lower()=='canpac' or ('pac' in region.lower()):
        # NORTHEAST NORTH PACIFIC OCEAN
                 # W State               S BC                  Yukon            Alaska              Panhandle        W of Stn P    S of Stn P
        bounds=[[-121.58785773894593, -121.52563946173123, -136.44599841272867, -152.2812542378697, -155.4204018657861, -145.5, -137.28310446004681, -121.58785773894593], # longitude
                [47.0,                49.45270791401021,   62.09039797199102,   62.31810681168754,  58.353439820448244, 50.0,   47.731407329459,     47.0]]                # latitude
    elif region.lower()=='canatl' or ('atl' in region.lower()):
        # NORTHWEST NORTH ATLANTIC OCEAN (Excluding Labrador Sea)
                       # Maine   NA                                         Greenland NFLD    Kippon     Belle Isl. NE LEdge  Cape St Charles Quebec Maine
        bounds=[[-69.2682391165346, -67.23317959216988, -43.91753, -43.91753, -52.45, -55.41667, -55.4293, -55.25,   -55.6215, -73.03521627003431, -69.2682391165346],          #longitude
                [44.43925382969355, 40.50195865876387, 40.50195865876387, 60.,       47.75,  51.6667,   51.9,     52.0333,  52.21667, 48.38081145574862, 44.43925382969355]]    # latitude
    elif region.lower()=='canarc':
        # Canadian Arctic
        bounds=[[-39.086275076348855, -39.086275076348855, -43.91753, -64.5786, -68.28738702218202, -76.36557811037876, -87.8757860794056, -159.8726925856259, -175.3028255881627, -175.3028255881627, -39.086275076348855],
              [90, 83.51414598039267, 60, 60, 57.64614783645243, 62.009526045642836, 66.17634712197676, 65.56495014602127, 66.30904793978928, 90, 90]]
    elif region.lower()=='hbay':
        # Hudson Bay
        bounds=[[-76.36557811037876, -74.90064255068444, -79.71400224682294, -98.1303349972659, -95.89805223963646, -87.8757860794056, -76.36557811037876],
                [62.009526045642836, 56.84073792181265, 50.22396544413005, 59.025044220488915, 63.57204432318876, 66.17634712197676, 62.009526045642836]]
    else:
        bounds=None
    
    if region.lower()=='arctic':
        # ARCTIC CIRCLE
        rmask=(lat>(66+34/60)).astype(int)
    elif bounds is None:
        # if not a defined region, don't mask
        return None
    else:
        bounds=np.array(bounds)
        # convert bounds to only have positive longitudes if lon only positive
        if np.nanmax(lon)>180 and np.nanmin(lon)>=0:
            bounds[0,bounds[0,:]<0]=360.+bounds[0,bounds[0,:]<0]
        
        # check which points fall within the boundary
        pth=path.Path(list(zip(bounds[0,:],bounds[1,:])))
        cntns = pth.contains_points(np.hstack((lon.flatten()[:,np.newaxis],lat.flatten()[:,np.newaxis]))).reshape(np.shape(lon))
        rmask=np.full(np.shape(cntns),np.nan)
        rmask[cntns]=1.0
    
    return rmask

def volMean(arr: np.ndarray, lev: int, mesh: Dict[str, np.ndarray], mask: np.ndarray) -> float:
    """
    Calculate a volume- or area-weighted average of an array.
    
    Inputs:
        arr - array over which to take the mean
        lev - vertical level at which to take the mean
        mesh - dictionary with area and/or volume of grid cells
        mask - mask to remove land
    Output:
        v_zhw - weighted mean of arr
    """

    A = arr.squeeze()
    if lev == -1:
        A = A*mesh['tmask']*np.tile(mask,(len(A),1,1))
        # weight by volume
        V_zhw = np.nansum(A*mesh['gridVol'])/np.nansum(mesh['gridVol']*np.tile(mask,(len(A),1,1)))
    elif lev != 0:
        A = A[0:(lev+1),:,:]*mesh['tmask']*np.tile(mask,(lev+1,1,1))
        V_zhw = np.nansum(A[0:(lev+1),:,:]*mesh['gridVol'])/np.nansum(mesh['gridVol']*np.tile(mask,(lev+1,1,1)))
    else:
        # no vertical weighting if only the surface
        V_zhw=np.nansum(A*mesh['gridArea']*mask)/np.nansum(mesh['gridArea']*mask)

    return V_zhw

def cantodsIce(vars: Dict[str, Dict[str, np.ndarray]], siconc: np.ndarray, sithick: np.ndarray, mesh: Dict[str, np.ndarray]) -> Dict[str, Dict[str, np.ndarray]]:
    """
    Calculate sea-ice area and volume from concentration and thickness.
    
    Inputs:
        vars - dictionary of sea-ice variables to be calculated
        siconc - sea ice concentration
        sithick - sea ice thickness
        mesh - dictionary containing grid cell area
    Output:
        vars - returns vars with desired variables calculated

    """

    for vK in vars.keys():
        # calculate ice area
        if 'siarea' in vK:
            siconc=(siconc.squeeze()/100.)*mesh['tmask']
            vars[vK]['data']=np.nansum(siconc*mesh['gridArea'],axis=(1,2))/1e6 # returns area in km2

        # calculate ice volume
        if 'sivol' in vK:
            sivol = siconc*mesh['gridArea']*sithick.squeeze()*mesh['tmask']
            vars[vK]['data']=np.nansum(sivol,axis=(1,2))/1e12 # returns area in 10^3 km3

    return vars

def rtdNetCDF(runid: str, timeseries: Dict[str, np.ndarray], var: str, args: argparse.Namespace) -> str:
    """
    Write/append data to the "RTD" netCDF file

    Inputs:
        runid - run name
        timseries - RTD timeseries
        var - variable to append
        args - arguments passed to cantods_rtd.py
    Outputs:
        ncFile - name of the file to which data were appended

    """
    ncFile=os.path.join(args.outdir,f'{runid}_{var}_timeseries.nc')
    ifFile=os.path.isfile(ncFile)
    if len(timeseries['data'])>0:
        ds=xr.Dataset.from_dict({'time':{"dims":('time'),'data':365.*timeseries['years'],"attrs": {"units": f"days since {args.year0}-01-01 00:00:00"}},
                                'year':{"dims":('time'),'data':timeseries['years'],'attrs':{'units':'simluation duration in years'}},
                                var:{"dims":('time'),'data':timeseries['data']}})
    if ifFile and not args.redo:
        if len(timeseries['data'])>0:
            # simply append to existing file (don't check times)
            # write temporary file with time as record dimension
            ds.to_netcdf(f'{ncFile}.tmp.nc', mode='w',unlimited_dims=['time'])
            # append along time dimension using nco
            err=subproc(f"ncrcat -h -H {ncFile} {ncFile}.tmp.nc -O {ncFile}")
            # delete temporary file
            err=subproc(f"rm -f {ncFile}.tmp.nc")
            print(f'\rAppended data to {ncFile}')
    else:
        if ifFile and args.redo:
            # delete existing file
            err=subproc(f"rm -f {ncFile}")
        if len(timeseries['data']) > 0:
            # write to new file with time as record dimension
            ds.to_netcdf(ncFile, mode='w',unlimited_dims=['time'])
            print(f'\rWrote to {ncFile}')
    
    return ncFile

def plotTimeseries(ncFile: str, varName: str, runid: str, args: Namespace) -> None:
    """
    Plot complete time series and link to public_html for viewing

    Inputs:
        ncFile - netCDF file containing timeseries data
        varName - variable name
        runid - run name
        args - arguments passed to cantods_rtd.py

    Outputs:
        Creates a plot showing the timeseries.
    """
    plt.close('all')
    with xr.open_dataset(ncFile,decode_times=False) as ds:
        plt.plot(ds['year'],ds[varName],'k.')
    plt.xlabel('Years Since Run Start')
    plt.ylabel(varName)
    plt.grid(visible=True)
    plt.tight_layout()

    pltFile=os.path.join(args.outdir,f"timeseries_{runid}_{varName}.png")
    plt.savefig(pltFile)
    if args.link:
        subscript(f"ln -fs {pltFile} /home/$(whoami)/public_html/CanTODS_diagnostics")
    plt.close('all')

def linkNotebook(args: Namespace) -> None:
    """
    Link the Jupyter ntoebook for plotting and run comparison.

    Inputs:
        args - arguments passed to cantods_rtd.py
    """

    for iR,runid in enumerate(args.runid):
        jptFile=f'cantods_diagnostics_{runid}.ipynb'
        if not os.path.isfile(f'{os.path.join(args.outdir,jptFile)}') or args.redo:
            if args.redo:
                err=subproc(f"rm -f {os.path.join(args.outdir,jptFile)}")
                err=subproc(f"rm -f /home/$(whoami)/public_html/CanTODS_diagnostics/{jptFile}")
            # copy Jupyter notebook template
            err=subproc(f"cp {os.path.join(args.source,'cantods_diagnostics.ipynb')} {os.path.join(args.outdir,jptFile)}")
            # update run name
            err=subproc(f'sed -i "s/RUNID/{runid}/g" {os.path.join(args.outdir,jptFile)}')
            # update where to find files
            sedOutDir=args.outdir.replace('/','\\/')
            err=subproc(f'sed -i "s/RTDPATH/{sedOutDir}/g" {os.path.join(args.outdir,jptFile)}')
            # update initial year
            err=subproc(f'sed -i "s/YEAR0/{args.year0}/g" {os.path.join(args.outdir,jptFile)}')

        if not os.path.isfile(f'/home/$(whoami)/public_html/CanTODS_diagnostics/{jptFile}'):
            # link to file
            err=subproc(f"ln -sf {os.path.join(args.outdir,jptFile)} /home/$(whoami)/public_html/CanTODS_diagnostics/{jptFile}")

    return