
import gc
import numpy as np
import os
import xarray as xr

# utilities for remapping
from cantods_diag.ctds_rtd_utils import cantodsFiles, cantodsIce, meshGrid, plotTimeseries, regionMask, relativeYear, rtdNetCDF, subproc

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
                if ('lab' in reg.lower()) or ('arc' in reg.lower()) or (reg.lower()=='arctic'):
                    # only proceeed for (potentially) ice-covered regions
                    if reg not in siregs: siregs.append(reg)
                    series[f'{sivar}_{reg}']={}

    # regional masks
    rmasks={}   
    # loop through each run id
    for iR,runid in enumerate(args.runid):
        # get mesh info
        if args.meshfile is not None:
            if len(list(args.meshfile)) == len(args.runid):
                meshSurf=meshGrid(args.meshfile[iR],runid,0,args)
            else:
                meshSurf=meshGrid(None,runid,0,args)
        else:
            meshSurf=meshGrid(args.meshfile,runid,0,args)

        print(f"\r  {f'Loaded mesh.':<75}",end='',flush=True)

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
                print(f"\r  {f'Loading: {os.path.basename(ifl)}':<75}",end='',flush=True)
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
            print(f"\r  {' ':<75}",end='',flush=True) # clear line
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
                            err=subproc(f"rm -f {pltFile}")
                            err=subproc(f"rm -f /home/$(whoami)/public_html/CanTODS_diagnostics/{os.path.basename(pltFile)}")
                    # plot time series if desired
                    if args.plot and os.path.isfile(ncF):
                        plotTimeseries(ncF,var,runid,args)
                
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
                                    err=subproc(f"rm -f {pltFile}")
                                    err=subproc(f"rm -f /home/$(whoami)/public_html/CanTODS_diagnostics/{os.path.basename(pltFile)}")
                            # plot time series if desired
                            if args.plot and os.path.isfile(ncF):
                                plotTimeseries(ncF,f'{var}_{reg}',runid,args) 