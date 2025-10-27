import gc
import numpy as np
import os
import xarray as xr

# utilities for remapping
from cantods_diag.ctds_rtd_utils import cantodsFiles, meshGrid, plotTimeseries, regionMask, relativeYear, rtdNetCDF, volMean

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
            if len(list(args.meshfile)) == len(args.runid):
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
                    try:
                        levZ.append(np.argmin(np.abs(meshAll['depths']-int(var.split('_')[1].split('m')[0]))))
                    except:
                        pass
        print(f"\r  {f'Loaded mesh.':<75}",end='',flush=True)

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
                print(f"\r  {f'Loading: {os.path.basename(tfl)}.':<75}",end='',flush=True)
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
                        print(f"\r  {f'Loading: {os.path.basename(tfl)}. {tcount+1:02}/{len(years)}':<75}",end='',flush=True)
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
                            err=subproc(f"rm -f {pltFile}")
                            err=subproc(f"rm -f /home/$(whoami)/public_html/CanTODS_diagnostics/{os.path.basename(pltFile)}")
                    # plot time series if desired
                    if args.plot and os.path.isfile(ncF):
                        plotTimeseries(ncF,var,runid,args)