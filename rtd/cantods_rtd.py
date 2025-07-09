"""
Computes diagnostics of the most recent CanTODS run.
Requires the py3_analysis_v3 environment.

Can be called stand-alone (post-run), or at runtime from
maestro (in CanNemo/rtd/canesm_nemo_runtime_diag.sh)

Simple plots are created by this script (showing entire time
series in netCDF file), with more detailed plots able to be
viewed using the Jupyter notebook which is linked to public_html.

Written by J. G. Izett (2024)
"""

import argparse
import json
import numpy as np
import os
import subprocess
import time

# import specific remapping functions
from cantods_diag.ctds_rtd_utils import linkNotebook, subproc
from cantods_diag.ctds_ice import calcIce
from cantods_diag.ctds_phys import calcPhys
from cantods_diag.ctds_trans import calcTransports

# argument parser; allows diagnostics to be calculated at runtime, or a posteriori
parser=argparse.ArgumentParser(description='Calculate CanTODS diagnostics.')

# main arguments for parser
parser.add_argument('-r','--runid',help='Run id(s) to analyze/compare. Can invoke multiple times.',action='append')
parser.add_argument('-p','--runpath',help='Path to directory where runs are stored (excluding runid). Can have multiple paths to match runids, or simply single path for all.',action='append')
parser.add_argument("-j","--json",help="Path to JSON file for run-specific options. Default: remap_canesm.json",default="cantods_rtd.json")
parser.add_argument("-J","--jout",help="Specify path to write a new JSON file with options used to run here. Default: do not write.",default="None")

# additional arguments to overwrite JSON
parser.add_argument('-m','--meshfile',help='Domain mesh file(s). Can have multiple paths to match runids, or simply single path for all. If not passed, finds file that matches',action='append',default=None)
parser.add_argument('-A','--all',help='Set True to run everything.',default=True)
parser.add_argument('-o','--outdir',help='Output directory for diagnostic netCDF files. Default is current directory.',default='./')
parser.add_argument('-R','--redo',help='If True, deletes original netCDF output files (if present) and creates fresh files, otherwise appends to existing files.',default=False)
parser.add_argument('-I','--ice',help='Calculate sea-ice parameters (volume and area). Not used if -A is True. 0: volume and area, 1: volume only, 2: area only, other: do not calculate',default=0)
parser.add_argument('-Y','--phys',help='Calculate mean temperature and salnity (3D and surface), max MLD, and mean SSH. Not used if -A is True. 0: calculate all variables, 1: T only, 2: S only, 3: MLD only, 4: SSH only. Other: do not calculate',default=0)
parser.add_argument('-T','--trans',help='Calculate volume transports through straits. Not used if -A is True. 0: calculate all transports, 1: Bering only, 2: Davis N only, 3: Davis S only, Other: do not calculate',default=0)
parser.add_argument('-z','--depths',help='Depths (in m) for calculating variables (as applicable). Pass as list, e.g., 0,250,1000,all',default='0,250,1000,all')
parser.add_argument('-x','--xlim',help='Set time limit for calculating diagnostics. All series are referenced to model run time (starting from year 0).',default=[-np.inf,np.inf])
parser.add_argument('-y','--year0',help='Start year for model run(s). Can either pass single value or use multiple times to match number of runids. If not passed, finds earliest date in files and assumes that to be the start.',default=None)
parser.add_argument('-g','--regions',help='List of regions in which to calculate means.',default='labsea,arctic,canpac,canatl')
parser.add_argument('-P','--plot',help='Logical flag to plot entire time series from netCDF file. Links plot to public_html.',default=False)
parser.add_argument('-L','--link',help='Logical flag to link to public_html.',default=True)
parser.add_argument('-s','--source',help='Source directory for file if not the working directory.',default='.')
parser.add_argument('-v','--version',help='File version. If not passed, looks for highest version.',default=None)

# parse arguments; convert to dictionary
prs=vars(parser.parse_args())

# read JSON file
if (prs['json'] is not None) and (prs['json'].lower()!="none") and os.path.isfile(prs['json']):
    with open(prs['json'],'r') as jFile:
        args=json.load(jFile)
    # replace JSON file values with value from command line argument if present
    for jJ in args.keys():
        if (jJ in prs.keys()) and (prs[jJ] is not None):
            try:
                # make sure string is not 'none'
                if prs[jJ].lower()!='none':
                    args[jJ]=prs[jJ]
            except:
                args[jJ]=prs[jJ]
    for jJ in ['json','jout']:
        args[jJ]=prs[jJ]    
else:
    # no JSON file, so just assign values from argument parser
    args=prs

# write JSON file with all options if desired
if (args['jout'] is not None) and (args['jout'].lower() != "none"):
    with open(args['jout'],'w') as jOut:
        args.pop('jout')
        args.pop('json')
        json.dump(args,jOut)

# convert args back to namespace object
args=argparse.Namespace(**args)

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
    args.trans=0
else:
    args.ice=int(args.ice)
    args.phys=int(args.phys)
    args.trans=int(args.trans)

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
args.link=np.logical_or(args.link=='True',args.link=='1')

# start year as integer
if args.year0 is not None:
    args.year0=int(args.year0)

# create output directory and link in public_html
args.outdir=os.path.abspath(args.outdir)    # ensure absolute path
os.makedirs(args.outdir,exist_ok=True)
if args.link:
    err=subproc(f"mkdir -p /home/$(whoami)/public_html/CanTODS_diagnostics/")

# Sea Ice Calculations
if args.ice in range(3):
    calcIce(args)

# Temperature and Salinity Calculations
if args.phys in range(5):
    calcPhys(args)

# Calculate volume transports
if args.trans in range(5):
    calcTransports(args)

# link jupyter notebook to public_html if not already
if args.link:
    linkNotebook(args)

end = time.time()
print(f'\n\nFinished CanTODS diagnostics.\nElapsed time:\n{(end-start)/60} minutes')
