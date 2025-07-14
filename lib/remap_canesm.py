"""
This script provides the wrapper to various remapping functions written to take
CMORIZED outputs and map them to a different domain. While the script was
originally written to take CMORized output from CanESM5 (hence the name), it should
work for any CMORized data, and has been tested on different re-analysis datasets.

Remapping (where appropriate) is done using CDO. Quantities may not be conserved.

This script can be run from command line using command line arguments, a JSON file,
or both. (Command line arguments overwrite JSON).

Specific functions include:
    rs_remap - remap NEMO restart files for use on a different domain
    ic_remap - remap NEMO initial condition files.
    bdy_slice - simply slice a file to extract boumdary forcing
                (i.e., where a child domain is simply a subset of the parent domain)
    bdy_remap - remap parent output onto a child boundary
    frc_slice - slice atmospheric forcing in time (to be used with weight files)
    sos_remap - remap sea-surface salinity for use with ln_ssr != 0
    rvr_remap - at the moment, a very rough river remapping that simply scales
                total volume to match the parent run

Written by Jonathan Izett (2024/2025)
"""

#TODO: Proper handling of rivers
#TODO: CDO commands can likely be made more efficient, rather than creating multiple intermediate files (but it is easier for debugging!)

import argparse
import json
import os
import sys
import time

# import specific remapping functions
from remap_functions.bdy_remap import bdy_remap, bdy_slc
from remap_functions.frc_remap import frc_slice, sos_remap
from remap_functions.ic_remap import ic_remap,rs_remap
from remap_functions.remap_utils import calc_nemo_chunk_dates
from remap_functions.rvr_remap import rvr_remap

# argument parser to get type of remapping to perform and JSON file
# Parser allows to be calculated at runtime, or offline
parser=argparse.ArgumentParser(description='Remap canesm outputs. Reguires the type of remapping be specified. Other options can either be read from a JSON file or command line. Command line overwrites options read from JSON file.')

# Main arguments for parser
parser.add_argument("-t","--type",help="Type of files to produce. Either:\n  -1 - nemo configuration info\n  0 - Initial conditions (default)\n  1 - Boundary conditions\n  2 - Atmospheric forcing (does not remap; only gets correct file time)\n  3 - Rivers (currently simple scaling)\n  4 - Surface salinity",type=int)
parser.add_argument("-j","--json",help="Path to JSON file for run-specific options. Default: remap_canesm.json",default="remap_canesm.json")
parser.add_argument("-J","--jout",help="Specify path to write a new JSON file with options used to run here. Default: do not write.",default="None")

# additional arguments to overwrite JSON
parser.add_argument('-y','--years',help='Years to process. It two years passed, processes years in range(year1,year2+1). If not two years, processes as a list.',action='append',type=int)
parser.add_argument('-a','--iaf_loop_year',help='Reference year for loop if offset for cyclical forcing.')
parser.add_argument('-A','--iaf_year_offset',help='Year offset if wanting to use cyclical forcing.')
parser.add_argument('-P','--parent_path',help='Full path to files containing model output from parent run, ending before the ensemble identifier.\ne.g., /fs/site5/eccc/crd/ccrn/model_output/CMIP6/final/CMIP6/CMIP/CCCma/CanESM5')
parser.add_argument('-p','--parent_name',help='Name of parent run to help identify output files. Required for finding files.\ne.g., CanESM5')
parser.add_argument('-e','--parent_ensemble',help='Ensemble identifier for parent run.')
parser.add_argument('-x','--parent_experiment',help='Experiment (e.g., piControl) for parent run.')
parser.add_argument('-g','--parent_grid',help='Parent meshfile. Only needed if remapping rivers.')
parser.add_argument('-o','--outfile',help='Prefix for output file. If None, default name is created based on type of file being produced.')
parser.add_argument('-m','--meshfile',help='Meshfile for remapping.')
parser.add_argument('-B','--bdyFile',help='Map to a specific boundary coordinate file (if type=1), e.g., to the Med. Otherwise maps to presumed north/south boudary.')
parser.add_argument('-F','--forcing',help='Type of forcing CanESM (default) or OMIP')
parser.add_argument('-i','--ic_ind',help='Index in file of desired time for initial condition.')
parser.add_argument('-H','--hot_start',help='If remapping IC, generate hotstart file from NEMO restart.')
parser.add_argument('-X','--Xdeg',help='Slice to contain X degrees either side of the boundary. Should be larger than the parent grid resolution to guarantee border.',type=float)
parser.add_argument('-R','--run_start_year',help='Run start year if type==-1')
parser.add_argument('-r','--run_start_month',help='Run start month if type==-1')
parser.add_argument('-l','--loop',help='Sequencer loop if type==-1')
parser.add_argument('-f','--nemo_freq_months',help='If type==-1')
parser.add_argument('-M','--mor',help='CMOR frequency/directory (e.g., 3hr or Amon) for atmospheric forcing.')
parser.add_argument('-v','--rvr',help='River file for remapping/scaling of freshwater inputs (if type==3)')
parser.add_argument('-c','--his2cmor',help='If 1, look for history file rather than CMORized file.')
parser.add_argument('-V','--flip_vel',help='Flip velocities in the "north" of the domain if regional fold')

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

# keep track of how long it takes to remap
start=time.time()

# get necessary info about run
if args.type==-1:
    calc_nemo_chunk_dates(args)
elif args.type==0:
    if args.hot_start==0:
        # restart from rest using history/CMORized files
        ic_remap(args)
    else:
        # hot start
        rs_remap(args)
    print(f'Done making IC files!\n{time.time()-start} s elapsed.')
elif args.type==1:
    fcount=bdy_remap(args)
    print(f'Done making the boundary files!\n{fcount} files created.\n{time.time()-start} s elapsed.')
elif args.type==2:
    fcount,fexpect=frc_slice(args)
    #fcount,fexpect=frc_remap(args)
    print(f'Done finding forcing files!\n({fcount}/{fexpect} found)\n{time.time()-start} s elapsed.')
elif args.type==3:
    print('Remapping/scaling rivers.')
    rvr_remap(args)
    print(f'Done river rivermapping.')
elif args.type==4:
    print('Remapping sos.')
    sos_remap(args)
    print('Done remapping sos.')
elif args.type==10:
    bdy_slc(args)
    print(f'Done extracting file as boundary input!')
else:
    sys.exit(f'ERROR: type must be -1, 0, 1 [/10], 2, 3, or 4...not {args.type}')
