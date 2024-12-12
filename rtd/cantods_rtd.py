"""
Computes diagnostics of the most recent CanTODS run.
Requires the py3_analysis_v2 environment.

Can be called stand-alone (post-run), or at runtime from
maestro (in CanNemo/rtd/canesm_nemo_runtime_diag.sh)

Simple plots are created by this script (showing entire time
series in netCDF file), with more detailed plots able to be
viewed using the Jupyter notebook which is linked to public_html.

Written by J. G. Izett (2024)
"""

import argparse
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

# argument parser
# allows to be calculated at runtime, or a posteriori
parser=argparse.ArgumentParser(description='Calculate CanTODS diagnostics.')

# add arguments
parser.add_argument('-r','--runid',help='Run id(s) to analyze/compare. Can invoke multiple times.',action='append',required=True)
parser.add_argument('-p','--runpath',help='Path to directory where runs are stored (excluding runid). Can have multiple paths to match runids, or simply single path for all.',action='append',required=True)
parser.add_argument('-m','--meshfile',help='Domain mesh file(s). Can have multiple paths to match runids, or simply single path for all. If not passed, finds file that matches',action='append',default=None)
parser.add_argument('-A','--all',help='Set True (default) to run everything.',default=True)
parser.add_argument('-o','--outdir',help='Output directory for diagnostic netCDF files. Default is current directory.',default='./')
parser.add_argument('-R','--redo',help='If True (default is False), deletes original netCDF output files (if present) and creates fresh files. Default is to append.',default=False)
parser.add_argument('-I','--ice',help='Calculate sea-ice parameters (volume and area). Not used if -A is True. 0: volume and area (default), 1: volume only, 2: area only, other: do not calculate',default=0)
parser.add_argument('-Y','--phys',help='Calculate mean temperature and salnity (3D and surface), max MLD, and mean SSH. Not used if -A is True. 0: calculate all variables (default), 1: T only, 2: S only, 3: MLD only, 4: SSH only. Other: do not calculate',default=0)
parser.add_argument('-z','--depths',help='Depths (in m) for calculating variables (as applicable). Pass as list, e.g., 0,250,1000,all',default='0,250,1000,all')
parser.add_argument('-x','--xlim',help='Set time limit for calculating diagnostics. All series are forced to 0-time.',default=[-np.inf,np.inf])
parser.add_argument('-y','--year0',help='Start year for model run(s). Can either pass single value or use multiple times to match number of runids. If not passed, finds earliest date in files and assumes that to be the start.',default=None)
parser.add_argument('-g','--regions',help='List of regions in which to calculate means.',default='labsea,arctic,pnw')
parser.add_argument('-P','--plot',help='Logical flag to plot entire time series from netCDF file. Links plot to public_html.',default=False)
parser.add_argument('-s','--source',help='Source directory for file if not the working directory.',default='.')
parser.add_argument('-v','--version',help='File version. If not passed, looks for highest version.',default=None)

#############
# FUNCTIONS #
#############

def nemoFileType(var):
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

def cantodsFiles(runid,nemoVar,args):
    """
    Get list of files that fall within desired range.
    """

    # NEMO file type
    ftype=nemoFileType(nemoVar)

    # search for all files that match
    if len(list(args.runpath)) == 1:
        globdir=os.path.join(args.runpath[0],runid,'data')
    else:
        globdir=os.path.join(np.array(args.runpath)[np.array(args.runid)==runid][0],runid,'data')
    flist=glob.glob(f'{globdir}/mc_{runid}_*_m01_5d_{ftype.lower()}.nc.*')

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

#TODO: delete if working properly
# def time2datetime(times):
#   """
#   Convert times from one format (e.g., np.datetime) to datetime.datetime
#   using pandas
#   """
#   datetime = pd.to_datetime(times,unit='s').to_pydatetime()

#   return datetime

def relativeYear(dset,args):
    # convert from date in file to years since start date
    with warnings.catch_warnings():
        warnings.simplefilter("ignore",category=RuntimeWarning)
        try:
            ftime=np.array(dset.indexes['time_counter'].values)
        except:
            ftime=np.array(dset.indexes['time_counter'].to_datetimeindex().values)
        #TODO: delete if working properly
        # try:
        #     ftime=np.array(time2datetime(dset.indexes['time_counter'].values))
        # except:
        #     ftime=np.array(time2datetime(dset.indexes['time_counter'].to_datetimeindex().values))
    
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

def meshGrid(meshfiles,runid,olev,args):
    """
    Get the grid cell areas from the meshfile.
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
                e3t = meshmask['e3t_0'].values.squeeze()[0,:,:,:] # z,y,x

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

def regionMask(region,lon,lat):
    """
    Get a 2D mask defining a given region.
    """
    if region=='labsea' or ('lab' in region.lower()):
                # Kulliq    Greenland NFLD    Kippon     Belle Isl. NE LEdge  Cape St Charles
        bounds=[[-64.5786, -43.91753, -52.45, -55.41667, -55.4293, -55.25,   -55.6215, -64.5786], #longitude
                [ 60.,      60.,       47.75,  51.6667,   51.9,     52.0333,  52.21667, 60.]]     # latitude
    else:
        bounds=None
    
    if bounds is None:
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

def volMean(arr,lev,mesh,mask):
    """
    Calculate a volume-weighted average
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

def cantodsIce(vars,siconc,sithick,mesh):
    """
    Calculate sea-ice area and volume from concentration and thickness.
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

def rtdNetCDF(runid,timeseries,var,args):
    """
    Write/append data to file
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
            subprocess.run(f'ncrcat -h -H {ncFile} {ncFile}.tmp.nc -O {ncFile}',shell=True)
            # delete temporary file
            subprocess.run(f'rm -f {ncFile}.tmp.nc',shell=True)
            print(f'\rAppended data to {ncFile}')
    else:
        if ifFile and args.redo:
            # delete existing file
            subprocess.run(f'rm -f {ncFile}',shell=True)
        if len(timeseries['data']) > 0:
            # write to new file with time as record dimension
            ds.to_netcdf(ncFile, mode='w',unlimited_dims=['time'])
            print(f'\rWrote to {ncFile}')
    
    return ncFile

def calcIce(args):
    """
    Calculate sea-ice volume and area
    """
    print('Calculating ice.')

    # initialize dictionaries
    series={}; siregs=[]
    for iV,sivar in enumerate(['siarea','sivol']):
        if args.ice in [0,iV]:
            series[sivar]={}
            for reg in args.regions:
                if 'lab' in reg.lower() or 'arc' in reg.lower() or 'gos' in reg.lower():
                    # only proceeed for (potentially) ice-covered regions
                    if reg not in siregs: siregs.append(reg)
                    series[f'{sivar}_{reg}']={}

    # regional masks
    rmasks={}   
    # loop through each run id
    for iR,runid in enumerate(args.runid):
        # get mesh info
        if args.meshfile is not None:
            if len(list(meshfile)) == len(args.runid):
                meshSurf=meshGrid(args.meshfile[iR],runid,0,args)
            else:
                meshSurf=meshGrid(None,runid,0,args)
        else:
            meshSurf=meshGrid(args.meshfile,runid,0,args)

        print(f"\r  {f'Loaded mesh.':<50}",end='',flush=True)

        # look for sea-ice concentration as it will be needed for both area and volume
        # thickness is in the same file
        iflist,args=cantodsFiles(runid,'siconc',args)
        
        # load sea-ice parameter(s)
        if len(iflist) > 0:
            # create empty arrays and assign values
            for iV,var in enumerate(['siarea','sivol']):
                if args.ice in [0,iV] and runid not in series[var].keys():
                    series[var][runid]={'years':np.array([]),'data':np.array([])}
                    for reg in siregs:
                        series[f'{var}_{reg}'][runid]={'years':np.array([]),'data':np.array([])}
            
            # load from files
            for ifl in iflist:
                print(f"\r  {f'Loading: {os.path.basename(ifl)}':<50}",end='',flush=True)
                with xr.open_dataset(ifl) as ds:
                    # get regional masks for domain
                    if len(list(rmasks.keys()))==0 and len(siregs) > 0:
                        rlon=ds['nav_lon'].values
                        rlat=ds['nav_lat'].values
                    for reg in siregs:
                        if reg not in rmasks.keys():
                            rmasks[reg]=regionMask(reg,rlon,rlat)
                    if 'rlon' in locals():
                        del(rlon); del(rlat); gc.collect()

                    # convert time to years since run start date
                    years=relativeYear(ds,args)
                    tID=range(len(years))

                    # restrict time based on x-limits (if set)
                    xlim=np.full(len(years),True,dtype='bool')
                    xlim[years<args.xlim[0]]=False
                    xlim[np.floor(years)>args.xlim[1]]=False
                    for yr in np.unique(years):
                        if yr in series[var][runid]['years']:
                            xlim[years==yr]=False
                    tID=np.array(tID)[xlim]
                    years=years[xlim]

                    # create empty arrays assign values
                    lvars={}; rlvars={}
                    for iV,var in enumerate(['siarea','sivol']):
                        if args.ice in [0,iV]:
                            lvars[var]={'years':years}
                            for reg in siregs:
                                if rmasks[reg] is not None:
                                    rlvars[f'{var}_{reg}']={'years':years}

                    # calculate area and/or volume
                    lvars=cantodsIce(lvars,ds['siconc'],ds['sithick'],meshSurf)
                    for reg in siregs:
                        if rmasks[reg] is not None:
                            rltmp={}
                            for iV,var in enumerate(['siarea','sivol']):
                                if args.ice in [0,iV]:
                                    rltmp[var]=rlvars[f'{var}_{reg}'].copy()
                            rltmp=cantodsIce(rltmp,rmasks[reg]*ds['siconc'],rmasks[reg]*ds['sithick'],meshSurf)
                            for iV,var in enumerate(['siarea','sivol']):
                                rlvars[f'{var}_{reg}']['data']=rltmp[var]['data'].copy()
                            del(rltmp)
                    for iV,var in enumerate(['siarea','sivol']):
                        if args.ice in [0,iV]:
                            lvars[var]['data']=lvars[var]['data'][tID]
                            for reg in siregs:
                                if rmasks[reg] is not None:
                                    rlvars[f'{var}_{reg}']['data']=rlvars[f'{var}_{reg}']['data'][tID]

                # append to time series
                for iV,var in enumerate(['siarea','sivol']):
                    if args.ice in [0,iV]:
                        series[var][runid]['years']=np.append(series[var][runid]['years'],lvars[var]['years'])
                        series[var][runid]['data']=np.append(series[var][runid]['data'],lvars[var]['data'])
                        for reg in siregs:
                            if rmasks[reg] is not None:
                                series[f'{var}_{reg}'][runid]['years']=np.append(series[f'{var}_{reg}'][runid]['years'],rlvars[f'{var}_{reg}']['years'])
                                series[f'{var}_{reg}'][runid]['data']=np.append(series[f'{var}_{reg}'][runid]['data'],rlvars[f'{var}_{reg}']['data'])

            # sort in time
            for iV,var in enumerate(['siarea','sivol']):
                if args.ice in [0,iV]:
                    iY=np.argsort(np.array(series[var][runid]['years']))
                    series[var][runid]['years']=np.array(series[var][runid]['years'])[iY]
                    series[var][runid]['data']=np.array(series[var][runid]['data'])[iY]
                    for reg in siregs:
                        if rmasks[reg] is not None:
                            series[f'{var}_{reg}'][runid]['years']=np.array(series[f'{var}_{reg}'][runid]['years'])[iY]
                            series[f'{var}_{reg}'][runid]['data']=np.array(series[f'{var}_{reg}'][runid]['data'])[iY]

            # save timeseries to file
            print(f"\r  {' ':<50}",end='',flush=True) # clear line
            for iV,var in enumerate(['siarea','sivol']):
                if args.ice in [0,iV]:
                    # only keep unique values
                    if len(series[var][runid]['data'])>0:
                        series[var][runid]['years'],iU=np.unique(series[var][runid]['years'],return_index=True)
                        series[var][runid]['data']=series[var][runid]['data'][iU]
                    ncF=rtdNetCDF(runid,series[var][runid],var,args)

                    # if redo, delete existing plot, even if not re-plotting
                    if args.redo:
                        pltFile=os.path.join(args.outdir,f"timeseries_{runid}_{var}.png")
                        if os.path.isfile(pltFile):
                            subprocess.run(f'rm -f {pltFile}',shell=True)
                            subprocess.run(f'rm -f /home/$(whoami)/public_html/CanTODS_diagnostics/{os.path.basename(pltFile)}',shell=True)
                    # plot time series if desired
                    if args.plot and os.path.isfile(ncF):
                        plotTimeseries(ncF,var,runid)
                
                    # plot individual regions
                    for reg in siregs:
                        if rmasks[reg] is not None:
                            # only keep unique values
                            if len(series[f'{var}_{reg}'][runid]['data']):
                                series[f'{var}_{reg}'][runid]['years'],iU=np.unique(series[f'{var}_{reg}'][runid]['years'],return_index=True)
                                series[f'{var}_{reg}'][runid]['data']=series[f'{var}_{reg}'][runid]['data'][iU]
                            ncF=rtdNetCDF(runid,series[f'{var}_{reg}'][runid],f'{var}_{reg}',args)

                            # if redo, delete existing plot, even if not re-plotting
                            if args.redo:
                                pltFile=os.path.join(args.outdir,f"timeseries_{runid}_{var}_{reg}.png")
                                if os.path.isfile(pltFile):
                                    subprocess.run(f'rm -f {pltFile}',shell=True)
                                    subprocess.run(f'rm -f /home/$(whoami)/public_html/CanTODS_diagnostics/{os.path.basename(pltFile)}',shell=True)
                            # plot time series if desired
                            if args.plot and os.path.isfile(ncF):
                                plotTimeseries(ncF,f'{var}_{reg}',runid)                   

def calcPhys(args):
    """
    Calculate temperature, salinity, MLD, and/or SSH.
    """
    print('Calculating physical diagnostics.')

    # initialize dictionaries
    aVars=[]; series={}
    for vV in ['T','S','MLD','SSH']:
        if vV in ['MLD','SSH']:
            aVars.append(vV)
            series[vV]={}
            for reg in args.regions:
                aVars.append(f'{vV}_{reg}')
                series[f'{vV}_{reg}']={}
        else:
            for iZ,zZ in enumerate(args.depths):
                if np.isinf(zZ):
                    dpth=f'3D'
                else:
                    dpth=f'{int(zZ):04}m'
                aVars.append(f'{vV}_{dpth}')
                series[f'{vV}_{dpth}']={}
                for reg in args.regions:
                    aVars.append(f'{vV}_{dpth}_{reg}')
                    series[f'{vV}_{dpth}_{reg}']={}

    # loop through each run id
    for iR,runid in enumerate(args.runid):
        # get mesh info
        if args.meshfile is not None:
            if len(list(meshfile)) == len(args.runid):
                meshAll=meshGrid(args.meshfile[iR],runid,-1,args)
                meshSurf=meshGrid(args.meshfile[iR],runid,0,args)
            else:
                meshAll=meshGrid(None,runid,-1,args)
                meshSurf=meshGrid(None,runid,0,args)
        else:
            meshAll=meshGrid(args.meshfile,runid,-1,args)
            meshSurf=meshGrid(args.meshfile,runid,0,args)
        levZ=[]
        for iZ,var in enumerate(aVars):
            if 'MLD' in var or 'SSH' in var:
                levZ.append(np.nan)
            else:
                if '3D' in var:
                    levZ.append(-1)
                else:
                    levZ.append(np.argmin(np.abs(meshAll['depths']-int(var.split('_')[1].split('m')[0]))))
        print(f"\r  {f'Loaded mesh.':<50}",end='',flush=True)

        # get list of files to load using glob
        tflist,args=cantodsFiles(runid,'thetao',args)
                
        # load temperature/salinity
        if len(tflist) > 0:
            # create empty arrays to assign values
            for iV,var in enumerate(aVars):
                if (args.phys==0 or ('T' in var and args.phys==1) or ('S_' in var and args.phys==2) or ('MLD' in var and args.phys==3) or ('SSH' in var and args.phys==4)) and (runid not in series[var].keys()):
                    series[var][runid]={'years':np.array([]),'data':np.array([])}
            
            # regional masks
            rmasks={}

            # load from files
            for tfl in tflist:
                print(f"\r  {f'Loading: {os.path.basename(tfl)}.':<50}",end='',flush=True)
                with xr.open_dataset(tfl) as ds:
                    # get regional masks for domain
                    if len(list(rmasks.keys()))==0 and len(args.regions) > 0:
                        rlon=ds['nav_lon'].values
                        rlat=ds['nav_lat'].values
                        rmasks['domain']=np.ones_like(rlon)
                    for reg in args.regions:
                        if reg not in rmasks.keys():
                            rmasks[reg]=regionMask(reg,rlon,rlat)
                    if 'rlon' in locals():
                        del(rlon); del(rlat); gc.collect()

                    # convert time to years since run start date
                    years=relativeYear(ds,args)
                    tID=range(len(years))

                    # restrict time based on x-limits (if set)
                    xlim=np.full(len(years),True,dtype='bool')
                    xlim[years<args.xlim[0]]=False
                    xlim[np.floor(years)>args.xlim[1]]=False
                    tID=np.array(tID)[xlim]
                    years=years[xlim]

                    # create empty arrays assign values
                    lvars={}
                    for iV,var in enumerate(aVars):
                        try:
                            reg=var.split('_')[2]
                        except:
                            reg=var.split('_')[-1]
                            if ('3D' in reg) or (reg==var) or ('m' in reg):
                                reg='domain'
                        if rmasks[reg] is not None:
                            if args.phys==0 or ('T' in var and args.phys==1) or ('S_' in var and args.phys==2) or ('MLD' in var and args.phys==3) or ('SSH' in var and args.phys==4):
                                lvars[var]={'years':years,'data':np.full(len(years),np.nan)}

                    # calculate means at each time step
                    # Loop through time because not doing so leads to crash given large annual files...at least I assume that is why it crashed...
                    for tcount,tid in enumerate(tID):
                        print(f"\r  {f'Loading: {os.path.basename(tfl)}. {tcount+1:02}/{len(years)}':<50}",end='',flush=True)
                        for iV,var in enumerate(aVars):
                            try:
                                reg=var.split('_')[2]
                            except:
                                reg=var.split('_')[-1]
                                if ('3D' in reg) or (reg==var) or ('m' in reg):
                                    reg='domain'
                            if rmasks[reg] is not None:
                                if 'T' in var and (args.phys in [0,1]):
                                    if '3D' in var:
                                        lvars[var]['data'][tcount]=volMean(ds['thetao'].isel(time_counter=tid).values,-1,meshAll,rmasks[reg])
                                    else:
                                        lvars[var]['data'][tcount]=volMean(ds['thetao'].isel(time_counter=tid,deptht=levZ[iV]).values,0,meshSurf,rmasks[reg])
                                elif 'S_' in var and (args.phys in [0,2]):
                                    if '3D' in var:
                                        lvars[var]['data'][tcount]=volMean(ds['so'].isel(time_counter=tid).values,-1,meshAll,rmasks[reg])
                                    else:
                                        lvars[var]['data'][tcount]=volMean(ds['so'].isel(time_counter=tid,deptht=levZ[iV]).values,0,meshSurf,rmasks[reg])
                                elif 'MLD' in var and (args.phys in [0,3]):
                                    lvars[var]['data'][tcount]=np.nanmax((ds['mlotst'].isel(time_counter=tid).values*rmasks[reg]).ravel())
                                elif 'SSH' in var and (args.phys in [0,4]):
                                    lvars[var]['data'][tcount]=volMean(ds['zos'].isel(time_counter=tid).values,0,meshSurf,rmasks[reg])
                        
                # append to time series
                for iV,var in enumerate(aVars):
                    if (args.phys==0 or ('T' in var and args.phys==1) or ('S_' in var and args.phys==2) or ('MLD' in var and args.phys==3) or ('SSH' in var and args.phys==4)) and var in lvars.keys():
                        series[var][runid]['years']=np.append(series[var][runid]['years'],lvars[var]['years'])
                        series[var][runid]['data']=np.append(series[var][runid]['data'],lvars[var]['data'])

            # sort in time
            for iV,var in enumerate(aVars):
                if args.phys==0 or ('T' in var and args.phys==1) or ('S_' in var and args.phys==2) or ('MLD' in var and args.phys==3) or ('SSH' in var and args.phys==4):
                    if len(series[var][runid]) > 0:
                        iY=np.argsort(np.array(series[var][runid]['years']))
                        series[var][runid]['years']=np.array(series[var][runid]['years'])[iY]
                        series[var][runid]['data']=np.array(series[var][runid]['data'])[iY]

            # save timeseries to file
            print(f"\r  {' ':<50}",end='',flush=True) # clear line
            for iV,var in enumerate(aVars):
                if (args.phys==0 or ('T' in var and args.phys==1) or ('S_' in var and args.phys==2) or ('MLD' in var and args.phys==3) or ('SSH' in var and args.phys==4)) and (var in series):
                    if len(series[var][runid]['data']) > 0:
                        # only keep unique values
                        series[var][runid]['years'],iU=np.unique(series[var][runid]['years'],return_index=True)
                        series[var][runid]['data']=series[var][runid]['data'][iU]
                    ncF=rtdNetCDF(runid,series[var][runid],var,args)

                    # delete existing plots if redo (even if not replotting)
                    if args.redo:
                        pltFile=os.path.join(args.outdir,f"timeseries_{runid}_{var}.png")
                        if os.path.isfile(pltFile):
                            subprocess.run(f'rm -f {pltFile}',shell=True)
                            subprocess.run(f'rm -f /home/$(whoami)/public_html/CanTODS_diagnostics/{os.path.basename(pltFile)}',shell=True)
                    # plot time series if desired
                    if args.plot and os.path.isfile(ncF):
                        plotTimeseries(ncF,var,runid)

def plotTimeseries(ncFile,varName,runid):
    """
    Plot complete time series and link to public_html for viewing
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
    subprocess.run(f'ln -fs {pltFile} /home/$(whoami)/public_html/CanTODS_diagnostics',shell=True)
    plt.close('all')

#################
# END FUNCTIONS #
#################

start=time.time()
print('\nStarting CanTODS diagnostics')

# parse arguments
args=parser.parse_args()

# override ice and temp values if calculating/plotting all
if args.all!=True:
    args.all=np.logical_or(args.all=='True',args.all=='1')
if args.all:
    args.ice=0
    args.phys=0
else:
    args.ice=int(args.ice)
    args.phys=int(args.phys)

# ensure run name not in runpath
for iN,rN in enumerate(args.runid):
    if rN in args.runpath[iN]:
        args.runpath[iN]=args.runpath[iN].replace(f'/{rN}','').replace(f'/data','')

# set xlim to a list
args.xlim=[int(args.xlim.split(',')[0].replace('[','')),int(args.xlim.split(',')[1].replace(']',''))]

# set depths to a list
args.depths=[float(depth) for depth in args.depths.lower().replace('[','').replace(']','').replace('all','inf').replace('3d','inf').split(',')]

# get regions as a list
args.regions=[reg for reg in args.regions.lower().replace('[','').replace(']','').split(',')]

# ensure redo is logical
args.redo=np.logical_or(args.redo=='True',args.redo=='1')

# ensure plot is logical
args.plot=np.logical_or(args.plot=='True',args.plot=='1')

# start year as integer
if args.year0 is not None:
    args.year0=int(args.year0)

# create output directory and link in public_html
args.outdir=os.path.abspath(args.outdir)    # ensure absolute path
os.makedirs(args.outdir,exist_ok=True)
subprocess.run(f'mkdir -p /home/$(whoami)/public_html/CanTODS_diagnostics/',shell=True)

# Sea Ice Calculations
if args.ice in range(3):
    calcIce(args)

# Temperature and Salinity Calculations
if args.phys in range(5):
    calcPhys(args)

# link jupyter notebook to public_html if not already
for iR,runid in enumerate(args.runid):
    jptFile=f'cantods_diagnostics_{runid}.ipynb'
    if not os.path.isfile(f'{os.path.join(args.outdir,jptFile)}') or args.redo:
        if args.redo:
            subprocess.run(f'rm -f {os.path.join(args.outdir,jptFile)}',shell=True)
            subprocess.run(f'rm -f /home/$(whoami)/public_html/CanTODS_diagnostics/{jptFile}',shell=True)
        # copy Jupyter notebook template
        subprocess.run(f"cp {os.path.join(args.source,'cantods_diagnostics.ipynb')} {os.path.join(args.outdir,jptFile)}",shell=True)
        # update run name
        subprocess.run(f'sed -i "s/RUNID/{runid}/g" {os.path.join(args.outdir,jptFile)}',shell=True)
        # update where to find files
        sedOutDir=args.outdir.replace('/','\\/')
        subprocess.run(f'sed -i "s/RTDPATH/{sedOutDir}/g" {os.path.join(args.outdir,jptFile)}',shell=True)
        # update initial year
        subprocess.run(f'sed -i "s/YEAR0/{args.year0}/g" {os.path.join(args.outdir,jptFile)}',shell=True)

    if not os.path.isfile(f'/home/$(whoami)/public_html/CanTODS_diagnostics/{jptFile}'):
        # link to file
        subprocess.run(f'ln -sf {os.path.join(args.outdir,jptFile)} /home/$(whoami)/public_html/CanTODS_diagnostics/{jptFile}',shell=True)

end = time.time()
print(f'\n\nFinished CanTODS diagnostics.\nElapsed time:\n{(end-start)/60} minutes')
