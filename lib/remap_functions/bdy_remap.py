"""
These functions allow remapping of CMORized files to boundaries of a regional domain.

  bdy_remap : remaps/interpolate to boundaries
  bdy_slice : simply slices fields from original domain
              (e.g., in case where regional domain is just a subset)

Written by: Jonathan Izett (2024/2025)

"""

from argparse import Namespace
import glob
import numpy as np
import os
from scipy.interpolate import LinearNDInterpolator
import sys
import time
import xarray as xr

# utilities for remapping
from remap_functions.remap_utils import checkMeshFile, cleanTmp, fillMiss2, getZ, his2cmor, intZ, matchFileYear, parseArgYears, remap, selYear, subproc

def bdy_remap(args: Namespace) -> int:
    """
    Remap CMORized file onto boundary.

    Inputs:
        args - args passed to remap_canesm.py
    Outputs:
        Writes boundary condition files and 
        returns number of files processed.
    """

    print(f'Finding and processing boundaries from {args.parent_name}_{args.parent_experiment}_{args.parent_ensemble}.')

    years=parseArgYears(args)

    # get output file name
    if args.outfile is None:
        # BDY gets replaced with north/south and YYYY gets replaced with the year
        outFile='obc_BDY_cantods025_yYYYY.nc'
    else:
        outFile=args.outfile

    # get boundaries from meshfile
    bN,bS=checkMeshFile(args.meshfile)
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
            if 'P' in bbase:
                # boundary max/min longitude (Note: may wrap around dateline or meridian)
                blon=np.mod(bdyFile.nav_lon.values,360.) ; blon2=blon.copy(); blon2[blon>180]=blon[blon>180]-360.
                bW=np.nanmin(blon[:,0])   # Western-most point of boundary   (   0 - 360)
                bE=np.nanmax(blon[:,-1])  # Eastern-most point of boundary   (   0 - 360)
                bW2=np.nanmin(blon2[:,0])   # Western-most point of boundary (-180 - 180)
                bE2=np.nanmax(blon2[:,-1])  # Eastern-most point of boundary (-180 - 180)
            else:
                bW=None; bE=None
    else:
        # assumed boundary from grid
        blist=['north','south']
        for bdy in blist:
            # sub-sample meshfile to only include the boundary and 10 rows around it
            if bdy=='south':
                err=subproc(f"ncks -h -d y,1,10 grd.tmp.nc -O bdy.south.nc")
            elif bdy=='north':
                err=subproc(f"ncks -h -d y,-11,-2 grd.tmp.nc -O bdy.north.nc")
            # get the maximum latitude extent of the boundaries to slice file (quicker processing)
            with xr.open_dataset(f'bdy.{bdy}.nc') as ds0:
                # boundary max/min latitude and longitude
                bN=np.nanmax([bN,np.nanmax(ds0['nav_lat'].values)])
                bS=np.nanmin([bS,np.nanmin(ds0['nav_lat'].values)])
                bE=None ; bW=None

    # loop through each variable and interpolate to the regional boundaries
    vars2interp=['thetao','so','uo','vo','zos']
    for iV,vV in enumerate(vars2interp):
        # get history file in format that matches CMORized files
        if int(args.his2cmor) != 1:
            varPath=os.path.join(args.parent_path,args.parent_experiment,args.parent_ensemble,f'Omon/{vV}/gn/v20190429/')
            # # find all files that match format
            # flist=np.array(sorted(glob.glob(os.path.join(varPath,f'{vV}_Omon_{args.parent_name}_{args.parent_experiment}_{args.parent_ensemble}_gn_*.nc'))))
        
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
            print(f'{vV} - {year:04}')
            if (args.iaf_year_offset is not None) and (args.iaf_loop_year is not None):
                fy=year + int(args.iaf_year_offset)
                yd=int(args.iaf_loop_year)-int(args.iaf_year_offset)
                yr=fy-int((year-1)/yd)*yd
                if int(args.his2cmor)==1:
                    file=his2cmor(vV,yr,outFile.replace('YYYY',f'{yr:04}'),args)
                    fdates=f'{yr}'
                    sameFile=False
                else:
                    file,fdates,sameFile=matchFileYear(yr,os.path.join(varPath,f'{vV}_Omon_{args.parent_name}_{args.parent_experiment}_{args.parent_ensemble}_gn_*.nc'),file0=file0)
            else:
                # find matching file
                if int(args.his2cmor)==1:
                    file=his2cmor(vV,year,outFile.replace('YYYY',f'{year:04}'),args)
                    fdates=f'{year}'
                    sameFile=False
                else:
                    file,fdates,sameFile=matchFileYear(year,os.path.join(varPath,f'{vV}_Omon_{args.parent_name}_{args.parent_experiment}_{args.parent_ensemble}_gn_*.nc'),file0=file0)
                yr=int(year)

            if (file is not None) and (not sameFile):
                file0=file
                            
                # subsample srcFile to be within +/- X degrees of the boundary
                if len(file) > 1:
                    fstr=''
                    for fle in file:
                        fstr+=f'{fle} '
                    err=subproc(f"ncrcat -h {fstr} -O {outFile}.{vV}.concat.tmp.nc")
                else:
                    err=subproc(f"ln -sf {file[0]} {outFile}.{vV}.concat.tmp.nc")
                with xr.open_dataset(file[0]) as src:
                    # find indices of region north and south of the boundary
                    try:
                        ilat=np.where(np.sum(np.logical_and(src['latitude']>=bS-args.Xdeg,src['latitude']<=bN+args.Xdeg),axis=1))[0]
                        latLen=len(src['latitude'].isel(i=0))
                    except:
                        ilat=np.where(np.logical_and(src['latitude']>=bS-args.Xdeg,src['latitude']<=bN+args.Xdeg))[0]
                        latLen=len(src['latitude'])
                    if bW is not None and bE is not None:
                        # get longitude limits as well
                        sln=np.mod(src['longitude'],360.)
                        if bE > bW:
                            # use positive longitudes
                            ilon=np.where(np.logical_and(sln>=bW-args.Xdeg,sln<=bE+args.Xdeg))[0]
                        else:
                            # use negative longitudes
                            sln[sln>180]=sln[sln>180]-360.
                            ilon=np.where(np.logical_and(sln>=bW-args.Xdeg,sln<=bE+args.Xdeg))[0]
                    else:
                        ilon=range(np.shape(src['longitude'])[1])
                    
                    if bW is not None and bE is not None:
                        # create index arrays in x and y
                        glt=src['latitude'].values; gln=np.mod(src['longitude'].values,360.)
                        if bW > bE:
                            gln[gln>180]=gln-360.
                        y,x=np.shape(glt)
                        Y=np.transpose(np.tile(range(y),(x,1)))
                        X=np.tile(range(x),(y,1))
                        # select only those values within the identified bounds
                        Y=Y[np.nanmin(ilat):np.nanmax(ilat)+1,:] ; X=X[np.nanmin(ilat):np.nanmax(ilat)+1,:]
                        glt=glt[np.nanmin(ilat):np.nanmax(ilat)+1,:] ; gln=gln[np.nanmin(ilat):np.nanmax(ilat)+1,:]
                        Y=Y[:,np.nanmin(ilon):np.nanmax(ilon)+1]; X=X[:,np.nanmin(ilon):np.nanmax(ilon)+1]
                        glt=glt[:,np.nanmin(ilon):np.nanmax(ilon)+1]; gln=gln[:,np.nanmin(ilon):np.nanmax(ilon)+1]
                        # interpolate the indices to the boundary coordinates
                        with xr.open_dataset(args.bdyFile) as grd:
                            # get coordinates
                            blt=grd['nav_lat'].values; bln=np.mod(grd['nav_lon'].values,360.)
                            if bW > bE:
                                bln[bln>180]=bln-360.
                            # create interpolants
                            yLin=LinearNDInterpolator(list(zip(gln.ravel(),glt.ravel())),Y.ravel())
                            xLin=LinearNDInterpolator(list(zip(gln.ravel(),glt.ravel())),X.ravel())
                            # interpolate to boundary coordinates
                            Yint=yLin(bln,blt)
                            Xint=xLin(bln,blt)
                            # round down from smallest and up from largest
                            ilat=[int(np.floor(np.nanmin(Yint))),int(np.ceil(np.nanmax(Yint)))]
                            ilon=[int(np.floor(np.nanmin(Xint))),int(np.ceil(np.nanmax(Xint)))]
                # add some extra latitude points if min and max the same, or only one point
                if (np.nanmax(ilat)-np.nanmin(ilat)) <= 1:
                    ilat=[np.nanmax([0,np.nanmin(ilat)-1]),np.nanmin([np.nanmax(ilat)+1,latLen])]
                if bW is not None and bE is not None:
                    err=subproc(f"ncks -h -d i,{np.nanmin(ilon)},{np.nanmax(ilon)} -d j,{np.nanmin(ilat)}.,{np.nanmax(ilat)}. {outFile}.{vV}.concat.tmp.nc -O {outFile}.{vV}.sliced.tmp.nc")
                    print(f"ncks -h -d i,{np.nanmin(ilon)},{np.nanmax(ilon)} -d j,{np.nanmin(ilat)}.,{np.nanmax(ilat)}. {outFile}.{vV}.concat.tmp.nc -O {outFile}.{vV}.sliced.tmp.nc")
                else:
                    err=subproc(f"ncks -h -d j,{np.nanmin(ilat)}.,{np.nanmax(ilat)}. {outFile}.{vV}.concat.tmp.nc -O {outFile}.{vV}.sliced.tmp.nc")
                    print(f"ncks -h -d j,{np.nanmin(ilat)}.,{np.nanmax(ilat)}. {outFile}.{vV}.concat.tmp.nc -O {outFile}.{vV}.sliced.tmp.nc")

                # subsample file to only include desired years
                selYear(mnY,mxY,f"{outFile}.{vV}.sliced.tmp.nc", f"{outFile}.{vV}.sliced2.tmp.nc")
                
                # fill any gaps in the sliced srcFile (two iterations)
                fillMiss2(f"{outFile}.{vV}.sliced2.tmp.nc",f"{outFile}.{vV}.filled.tmp.nc")

                # interpolate vertically (if not SSH)
                if vV == 'zos':
                    # interpolate filled src file to boundary points
                    err=subproc(f"mv {outFile}.{vV}.filled.tmp.nc {outFile}.{vV}.z.tmp.nc")
                else:
                    getZ(args.meshfile,outFile)
                    intZ(f"{outFile}.{vV}.filled.tmp.nc", f"{outFile}.{vV}.z.tmp.nc", f"{outFile}.onlyz.tmp.nc")

                # remap to meshfiles to boundaries
                for bdy in blist:
                    remap(f"bdy.{bdy}.nc", f"{outFile}.{vV}.z.tmp.nc", f"{outFile.replace('BDY',bdy)}_y{fdates}.{vV}.tmp.nc")

                # remove intermediate variable files
                cleanTmp(f"{outFile}.{vV}*.tmp.nc*")
            
            if file is not None:
                # slice files to match the desired date range
                for bdy in blist:
                    selYear(yr,yr,f"{outFile.replace('BDY',bdy)}_y{fdates}.{vV}.tmp.nc", f"{outFile.replace('BDY',bdy).replace('YYYY',f'{year}')}.{vV}.nc")
            else:
                print('No files found.')
        
        # remove intermediate files
        cleanTmp(f"{outFile.replace('BDY','*')}_y*-*.{vV}*")
    
    # remove all remaining intermediate files generated above
    cleanTmp(f"*.tmp.*nc")

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
                err=subproc(f"cdo merge {flist.strip()} {outFileBdy.replace('YYYY',f'{year:04}')}")
                # remove individual variable files
                err=subproc(f"rm -f {flist}")
            
        # rename variables, dimensions, etc. for NEMO
        dimNames={'time':'t','lev':'z','i':'x','j':'y'}
        varNames={'longitude':'nav_lon','latitude':'nav_lat','i':'x','j':'y','thetao':'votemper',
                'so':'vosaline','uo':'vozocrtx','vo':'vomecrty','zos':'sossheig'}
        fcount=0
        for fF in sorted(glob.glob(f"{outFileBdy.replace('YYYY','*')}")):
            # rename depth variable (shouldn't matter since all the same depth...)
            # and add attributes
            err=subproc(f"cdo --no_history setattribute,thetao@grid=T,so@grid=T,zos@grid=T,uo@grid=U,vo@grid=V {fF} {fF}2")
            err=subproc(f"ncrename -h -v .lev,deptht {fF}2 -O {fF}")
            err=subproc(f"rm -f {fF}2")

            # rename other variables and dimensions
            for iD,dD in enumerate(dimNames.keys()):
                err=subproc(f"ncrename -h -d .{dD},{dimNames[dD]} {fF} -O {fF}")
            for iV,vV in enumerate(varNames.keys()):
                err=subproc(f"ncrename -h -v .{vV},{varNames[vV]} {fF} -O {fF}")
            
            # count files
            fcount+=1

    # remove sliced meshfiles and intermediate files
    err=subproc(f"rm -f {outFile.replace('BDY',bdy).replace('YYYY',f'{year}')}.{vV}.nc")
    cleanTmp(f"{outFile.replace('BDY',bdy)}*.tmp.nc")
    for bdy in blist:
        err=subproc(f"rm -f bdy.{bdy}.nc")

    return fcount

def bdy_slc(args: Namespace) -> int:
    """
    Take global CMORized file and simply slice it in time and at a desired index
    to match the chosen boundary.
    
    Useful, e.g., when the child domain is simply a subset of the parent and no
    interpolation is required.
    
    Inputs:
        args - args passed to remap_canesm.py
    Outputs:
        Writes sliced boundary files and
        returns number of files processed.

    """
    
    print(f'Finding and processing boundaries from {args.parent_name}_{args.parent_experiment}_{args.parent_ensemble}.')

    years=parseArgYears(args)

    # get output file name
    if args.outfile is None:
        # BDY gets replaced with north/south and YYYY gets replaced with the year
        outFile='obc_BDY_cantods025_yYYYY.nc'
    else:
        outFile=args.outfile

    # extract boundary indices from bdyFile
        # identify boundary with file name
    bbase=os.path.basename(args.bdyFile).replace('.nc','')
    blist=[bbase]
    with xr.open_dataset(args.bdyFile) as bdyFile:
        # write to temporary file
        msh=xr.Dataset.from_dict(
            {'nav_lon':{'dims':('y','x'),'data':bdyFile.nav_lon.values,'attrs':{'_CoordinateAxisType':'Lon','units':'degrees_east'}},
            'nav_lat':{'dims':('y','x'),'data':bdyFile.nav_lat.values,'attrs':{'_CoordinateAxisType':'Lat','units':'degrees_north'}},
            'nav_ones':{'dims':('time_counter','y','x'),'data':np.ones((1,np.shape(bdyFile.nav_lat.values)[0],np.shape(bdyFile.nav_lat.values)[1])),'attrs':{'coordinates':'nav_lat nav_lon'}},
            'time_counter':{'dims':('time_counter'),'data':[0.]}})
        # write to temporary netCDF file
        msh.to_netcdf(f'bdy.{bbase}.nc')
    # get coordinates, include buffer zone, write to file
    with xr.open_dataset(args.bdyFile) as bdyFile:
        # get range of i and j coordinates, padding by 2
        # max. is restricted later to be no more than the size of the variable
        iMin=np.nanmax([np.nanmin(bdyFile.nbit.values-1),0])
        iMax=np.nanmax(bdyFile.nbit.values)+1
        jMin=np.nanmax([np.nanmin(bdyFile.nbjt.values-1),0])
        jMax=np.nanmax(bdyFile.nbjt.values)+1
    # loop through each variable and extract at the regional boundaries
    vars2interp=['thetao','so','uo','vo','zos']
    for iV,vV in enumerate(vars2interp):
        # get history file in format that matches CMORized files
        if int(args.his2cmor) != 1:
            varPath=os.path.join(args.parent_path,args.parent_experiment,args.parent_ensemble,f'Omon/{vV}/gn/v20190429/')
        
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
            print(f'{vV} - {year:04}')
            if (args.iaf_year_offset is not None) and (args.iaf_loop_year is not None):
                fy=year + int(args.iaf_year_offset)
                yd=int(args.iaf_loop_year)-int(args.iaf_year_offset)
                yr=fy-int((year-1)/yd)*yd
                if int(args.his2cmor)==1:
                    file=his2cmor(vV,yr,outFile.replace('YYYY',f'{yr:04}'),args)
                    fdates=f'{yr}'
                    sameFile=False
                else:
                    file,fdates,sameFile=matchFileYear(yr,os.path.join(varPath,f'{vV}_Omon_{args.parent_name}_{args.parent_experiment}_{args.parent_ensemble}_gn_*.nc'),file0=file0)
            else:
                # find matching file
                if int(args.his2cmor)==1:
                    file=his2cmor(vV,year,outFile.replace('YYYY',f'{year:04}'),args)
                    fdates=f'{year}'
                    sameFile=False
                else:
                    file,fdates,sameFile=matchFileYear(year,os.path.join(varPath,f'{vV}_Omon_{args.parent_name}_{args.parent_experiment}_{args.parent_ensemble}_gn_*.nc'),file0=file0)
                yr=int(year)

            if (file is not None) and (not sameFile):
                file0=file
                # slice file to match boundary file indices +/- 2
                # subsample srcFile to be within +/- X degrees of the boundary
                if len(file) > 1:
                    fstr=''
                    for fle in file:
                        fstr+=f'{fle} '
                    err=subproc(f"ncrcat -h {fstr} -O {outFile}.{vV}.concat.tmp.nc")
                else:
                    err=subproc(f"ln -sf {file[0]} {outFile}.{vV}.concat.tmp.nc")
                with xr.open_dataset(f"{outFile}.{vV}.concat.tmp.nc") as ccat:
                    # restrict length of variables
                    iMax2=np.nanmin([iMax,np.shape(ccat[vV].values)[3]])
                    jMax2=np.nanmin([jMax,np.shape(ccat[vV].values)[2]])
                err=subproc(f"ncks -h -d i,{iMin},{iMax2} -d j,{jMin},{jMax2} {outFile}.{vV}.concat.tmp.nc -O {outFile}.{vV}.sliced.tmp.nc")
                print(f"ncks -h -d i,{iMin},{iMax2} -d j,{jMin},{jMax2} {outFile}.{vV}.concat.tmp.nc -O {outFile}.{vV}.sliced.tmp.nc")

                # subsample file to only include desired years
                selYear(mnY,mxY,f"{outFile}.{vV}.sliced.tmp.nc", f"{outFile}.{vV}.sliced2.tmp.nc")
                
                # fill any gaps in the sliced srcFile (two iterations)
                fillMiss2(f"{outFile}.{vV}.sliced2.tmp.nc",f"{outFile}.{vV}.filled.tmp.nc")

                # interpolate vertically (if not SSH)
                if vV == 'zos':
                    # interpolate filled src file to boundary points
                    err=subproc(f"mv {outFile}.{vV}.filled.tmp.nc {outFile}.{vV}.z.tmp.nc")
                else:
                    getZ(args.meshfile,outFile)
                    intZ(f"{outFile}.{vV}.filled.tmp.nc", f"{outFile}.{vV}.z.tmp.nc", f"{outFile}.onlyz.tmp.nc")

                # remap to boundary file
                for bdy in blist:
                    remap(f"bdy.{bdy}.nc", f"{outFile}.{vV}.z.tmp.nc", f"{outFile.replace('BDY',bdy)}_y{fdates}.{vV}.tmp.nc")

                # remove intermediate variable files
                cleanTmp(f"{outFile}.{vV}*.tmp.nc*")
            
            if file is not None:
                # slice files to match the desired date range
                for bdy in blist:
                    selYear(yr,yr,f"{outFile.replace('BDY',bdy)}_y{fdates}.{vV}.tmp.nc", f"{outFile.replace('BDY',bdy).replace('YYYY',f'{year}')}.{vV}.nc")
            else:
                print('No files found.')
       
        # remove intermediate files
        cleanTmp(f"{outFile.replace('BDY','*')}_y*-*.{vV}*")
    
    # remove all remaining intermediate files generated above
    cleanTmp(f"*.tmp.*nc")

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
                err=subproc(f"cdo merge {flist.strip()} {outFileBdy.replace('YYYY',f'{year:04}')}")
                # remove individual variable files
                err=subproc(f"rm -f {flist}")
            
        # rename variables, dimensions, etc. for NEMO
        dimNames={'time':'t','lev':'z','i':'x','j':'y'}
        varNames={'longitude':'nav_lon','latitude':'nav_lat','i':'x','j':'y','thetao':'votemper',
                'so':'vosaline','uo':'vozocrtx','vo':'vomecrty','zos':'sossheig'}
        fcount=0
        for fF in sorted(glob.glob(f"{outFileBdy.replace('YYYY','*')}")):
            # rename depth variable (shouldn't matter since all the same depth...)
            # and add attributes
            err=subproc(f"cdo --no_history setattribute,thetao@grid=T,so@grid=T,zos@grid=T,uo@grid=U,vo@grid=V {fF} {fF}2")
            err=subproc(f"ncrename -h -v .lev,deptht {fF}2 -O {fF}")
            err=subproc(f"rm -f {fF}2")

            # rename other variables and dimensions
            for iD,dD in enumerate(dimNames.keys()):
                err=subproc(f"ncrename -h -d .{dD},{dimNames[dD]} {fF} -O {fF}")
            for iV,vV in enumerate(varNames.keys()):
                err=subproc(f"ncrename -h -v .{vV},{varNames[vV]} {fF} -O {fF}")
            
            # count files
            fcount+=1

    # remove sliced meshfiles and intermediate files
    err=subproc(f"rm -f {outFile.replace('BDY',bdy).replace('YYYY',f'{year}')}.{vV}.nc")
    cleanTmp(f"{outFile.replace('BDY',bdy)}*.tmp.nc")
    for bdy in blist:
        err=subproc(f"rm -f bdy.{bdy}.nc")

    return fcount


# NOTE: this function is old and may not be useful; at a minimum it needs refactoring
#       ...preserving for now in case it does come in use later
# def bdy_time_slc(args):
#     """
#     Take global field CMORized file and simply slice it in time to match the model run.
    
#     Useful, e.g., when wanting to run with a (pre-existing) weight file, rather than
#     pre-processing
    
#     Inputs:
#         args - args passed to remap_canesm.py
#     Outputs:
#         Writes sliced files to use with a boundary weight file.

#     """
#     print(f'Finding and processing boundaries from {args.parent_name}_{args.parent_experiment}_{args.parent_ensemble}.')

#     years=parseArgYears(args)

#     # get output file name
#     if args.outfile is None:
#         # BDY gets replaced with north/south and YYYY gets replaced with the year
#         outFile='obc_cantods025_yYYYY.nc'
#     else:
#         outFile=args.outfile

#     # loop through each variable and interpolate to the regional boundaries
#     vars2interp=['thetao','so','uo','vo','zos']
#     for iV,vV in enumerate(vars2interp):
#         varPath=os.path.join(args.parent_path,args.parent_experiment,args.parent_ensemble,f'Omon/{vV}/gn/v20190429/')

#         # find all files that match format
#         flist=np.array(sorted(glob.glob(os.path.join(varPath,f'{vV}_Omon_{args.parent_name}_{args.parent_experiment}_{args.parent_ensemble}_gn_*.nc'))))
        
#         # identify and process file based on desired year
#         # slow on first pass, but then much quicker on subsequent years if reading from same file
#         file0=''
#         for year in years:
#             print(f'{vV} - {year}')
#             file,fdates,sameFile=matchFileYear(year,flist,file0=file0)
#             if (file is not None):
#                 # concatenate files if necessary
#                 if len(file) > 1:
#                     fstr=''
#                     for fle in file:
#                         fstr+=f'{fle} '
#                     subprocess.run(f'ncrcat -h {fstr} -O {outFile}.{vV}.concat.tmp.nc',shell=True)
#                 else:
#                     subprocess.run(f'ln -s {file[0]} {outFile}.{vV}.concat.tmp.nc',shell=True)
                
#                 # subsample file to only include desired years
#                 subprocess.run(f"cdo --no_history selyear,{min(years)}/{max(years)} {outFile}.{vV}.concat.tmp.nc {outFile}.{vV}.sliced2.tmp.nc",shell=True)
                
#                 # fill any gaps in the sliced srcFile (two iterations)
#                 subprocess.run(f'cdo --no_history fillmiss2,2 {outFile}.{vV}.sliced2.tmp.nc {outFile}.{vV}.filled.tmp.nc',shell=True)

#                 # interpolate vertically (if not SSH)
#                 if vV == 'zos':
#                     # interpolate filled src file to boundary points
#                     subprocess.run(f"mv {outFile}.{vV}.filled.tmp.nc {outFile.replace('YYYY',f'{year}')}.{vV}.nc",shell=True)
#                 else:
#                     getZ(args.meshfile,outFile)
#                     subprocess.run(f"cdo --no_history intlevelx$(cdo -s showlevel {outFile}.onlyz.tmp.nc | tr ' ' ',') {outFile}.{vV}.filled.tmp.nc {outFile.replace('YYYY',f'{year}')}.{vV}.nc",shell=True)

#                 # remove intermediate files
#                 subprocess.run(f"rm -f {outFile}.*.tmp.nc*",shell=True)    

#     # remove all remaining intermediate files generated above
#     subprocess.run(f"rm -f {outFile}*.tmp.nc*",shell=True)

#     # now concatenate physical variables and rename currents
#     print(f'\rConcatenating files')
#     for year in years:
#         flist=''; ccount=0
#         for vV in vars2interp:
#             if os.path.isfile(f"{outFile.replace('YYYY',f'{year}')}.{vV}.nc"):
#                 flist+=f"{outFile.replace('YYYY',f'{year}')}.{vV}.nc "
#                 ccount+=1
#         if ccount > 0:
#             subprocess.run(f"cdo merge {flist} {outFile.replace('YYYY',f'{year}')}",shell=True)
#             # remove individual variable files
#             subprocess.run(f'rm -f {flist}',shell=True)
        
#     # rename variables, dimensions, etc. for NEMO
#     dimNames={'time':'t','lev':'z','i':'x','j':'y'}
#     varNames={'longitude':'nav_lon','latitude':'nav_lat','i':'x','j':'y','thetao':'votemper',
#             'so':'vosaline','uo':'vozocrtx','vo':'vomecrty','zos':'sossheig'}
#     fcount=0
#     for fF in sorted(glob.glob(f"{outFile.replace('YYYY','*')}")):
#         # rename depth variable (shouldn't matter since all the same depth...)
#         # and add attributes
#         subprocess.run(f'cdo --no_history setattribute,thetao@grid=T,so@grid=T,zos@grid=T,uo@grid=U,vo@grid=V {fF} {fF}2',shell=True)
#         subprocess.run(f'ncrename -h -v .lev,deptht {fF}2 -O {fF}',shell=True)
#         subprocess.run(f'rm -f {fF}2',shell=True)

#         # rename other variables and dimensions
#         for iD,dD in enumerate(dimNames.keys()):
#             subprocess.run(f'ncrename -h -d .{dD},{dimNames[dD]} {fF} -O {fF}',shell=True)
#         for iV,vV in enumerate(varNames.keys()):
#             subprocess.run(f'ncrename -h -v .{vV},{varNames[vV]} {fF} -O {fF}',shell=True)
        
#         # count files
#         fcount+=1

#     # remove sliced meshfiles and intermediate files
#     subprocess.run(f"rm -f {outFile.replace('YYYY','*')}.*.nc",shell=True)

#     return fcount

