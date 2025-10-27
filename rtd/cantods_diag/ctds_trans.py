import numpy as np
import os
import xarray as xr

# utilities for remapping
from cantods_diag.ctds_rtd_utils import cantodsFiles, meshGrid, plotTimeseries, relativeYear, rtdNetCDF

def calcTransports(args):
    """
    Calculate volume transports through Bering and Davis Straits
    """
    if args.trans==0:
        tlist=[1,2,3,4]
    else:
        try:
            tlist=list(args.trans)
        except:
            tlist=[args.trans]
    tpvars=[None,'vol-Bering','vol-Davis','vol-CC','vol-CUC']

    series={}
    for iTrans in tlist:
        series[tpvars[iTrans]]={}
        
    # get list of files to load using glob
    for iR,runid in enumerate(args.runid):
        tpflist,args=cantodsFiles(runid,'vo',args)
        if len(tpflist)>0:
            for iTrans in tlist:
                series[tpvars[iTrans]][runid]={'years':np.array([]),'data':np.array([])}

            # load from files
            for tpfl in tpflist:
                print(f"\r  {f'Loading: {os.path.basename(tpfl)}.':<75}",end='',flush=True)
                mmFile=tpfl.replace('1m_grid_v','mesh_mask')
                # Get mask file. TODO: make use of access
                while not os.path.isfile(mmFile):
                    mmInt=int(mmFile[-3::])-1
                    if mmInt == 0:
                      break
                    else:
                      mmFile=mmFile.replace(mmFile[-3::],f'{mmInt:03}')

                with xr.open_dataset(mmFile) as mm:
                #with xr.open_dataset(tpfl.replace('1m_grid_v','mesh_mask')) as mm:
                    with xr.open_dataset(tpfl) as ds:
                        # convert time to years since run start date
                        years=relativeYear(ds,args)
                        tID=range(len(years))

                        # restrict time based on x-limits (if set)
                        xlim=np.full(len(years),True,dtype='bool')
                        xlim[years<args.xlim[0]]=False
                        xlim[np.floor(years)>args.xlim[1]]=False
                        tID=np.array(tID)[xlim]
                        years=years[xlim]

                        # file dimensions
                        tdim=ds.sizes['time_counter']
                        zdim=ds.sizes['depthv']

                        for iTrans in tlist:                            
                            tpvar=tpvars[iTrans]
                            print(f"\r  {f'Calculating: {tpvar}.':<75}",end='',flush=True)
                            if iTrans==1:
                                # Bering St
                                lon0=-170.6         # westernmost point of strait
                                lon1=-167.8         # easternmost point of strait
                                lat=65.+45./60.     # latitude of strait
                                dmin=0
                                dmax=500
                            elif iTrans==2:
                                # Davis St
                                lon0=-61.52   # westernmost point of strait
                                lon1=-53.66   # easternmost point of strait
                                lat=66.65     # latitude of strait
                                dmin=0
                                dmax=2700
                            elif iTrans==3:
                                # Cal. Current   # Newport Line
                                lon0=-126
                                lon1=-123
                                lat=44.65
                                dmin=0
                                dmax=150
                            elif iTrans==4:
                                # Cal. Undercurrent   # Newport Line
                                lon0=-126
                                lon1=-123
                                lat=44.65
                                dmin=150
                                dmax=2000

                            # find longitudes that falls within the desired range (get a mask)
                            dlon=ds.nav_lon.values
                            if lon0 > lon1:
                                iLon=np.logical_or(dlon>=lon0,dlon<=lon1)
                            else:
                                iLon=np.logical_and(dlon>=lon0,dlon<=lon1)
                            # determine which horizontal cells fall within the desired range
                            iX=np.where(np.sum(iLon,axis=0)>0)[0]
                            # get desired depth indices
                            kk=np.where(np.logical_and(dmin<=ds['depthv'],ds['depthv']<=dmax).values)[0]

                            # create empty array for assigning transport
                            ldata=np.zeros(len(years))
                            for ii,ix in enumerate(iX):
                                # find the nearest latitude in that column
                                iy=np.argmin(np.abs(ds['nav_lat'].isel(x=ix).values-lat))
                                proceed=False
                                if lon0>lon1:
                                    if (dlon[iy,ix]>=lon0) or (dlon[iy,ix]<=lon1):
                                        proceed=True
                                else:
                                    if (dlon[iy,ix]>=lon0) and (dlon[iy,ix]<=lon1):
                                        proceed=True
                                if proceed and np.abs(ds['nav_lat'].isel(x=ix,y=iy).values-lat) < 0.2:
                                    if mm['vmask'].isel(nav_lev=kk[0],time_counter=0,x=ix,y=iy).values==1:
                                        print(f"\r  {f'Calculating: {tpvar} {ii+1}/{len(iX)} (y={int(np.nanmean(years))})':<75}",end='',flush=True)
                                        dprod=ds['vo'].isel(x=ix,y=iy,depthv=kk)*ds['e3v'].isel(x=ix,y=iy,depthv=kk)*mm['e1v'].isel(x=ix,y=iy)*(mm['vmask'].isel(time_counter=0,x=ix,y=iy,nav_lev=kk).values)
                                        ldata[:]+=dprod.sum(dim='depthv').values # volumetric flow in m3/s
                            # convert to Sv and restrict in time
                            ldata=ldata[xlim]*1e-6

                            # append to time series
                            series[tpvar][runid]['years']=np.append(series[tpvar][runid]['years'],years)
                            series[tpvar][runid]['data']=np.append(series[tpvar][runid]['data'],ldata)

            for iTrans in tlist:
                tpvar=tpvars[iTrans]
                print(f"\r  {f'Finalizing {tpvar}':<75}",end='',flush=True)
                # sort in time
                print(f"\r  {f'Sorting {tpvar} in time':<75}",end='',flush=True)
                if len(series[tpvar][runid]) > 0:
                    iY=np.argsort(np.array(series[tpvar][runid]['years']))
                    series[tpvar][runid]['years']=np.array(series[tpvar][runid]['years'])[iY]
                    series[tpvar][runid]['data']=np.array(series[tpvar][runid]['data'])[iY]

                # save timeseries to file
                print(f"\r  {f'Saving {tpvar}':<75}",end='',flush=True)
                if len(series[tpvar][runid]['data']) > 0:
                    # only keep unique values
                    series[tpvar][runid]['years'],iU=np.unique(series[tpvar][runid]['years'],return_index=True)
                    series[tpvar][runid]['data']=series[tpvar][runid]['data'][iU]
                ncF=rtdNetCDF(runid,series[tpvar][runid],tpvar,args)

                # delete existing plots if redo (even if not replotting)
                if args.redo:
                    pltFile=os.path.join(args.outdir,f"timeseries_{runid}_{tpvar}.png")
                    if os.path.isfile(pltFile):
                        err=subproc(f"rm -f {pltFile}")
                        err=subproc(f"rm -f /home/$(whoami)/public_html/CanTODS_diagnostics/{os.path.basename(pltFile)}")
                # plot time series if desired
                if args.plot and os.path.isfile(ncF):
                    print(f"\r  {f'Plotting {tpvar}':<75}",end='',flush=True)
                    plotTimeseries(ncF,tpvar,runid,args)
    
    print(f"\r  {f'Done calculating transports.':<75}",end='',flush=True)