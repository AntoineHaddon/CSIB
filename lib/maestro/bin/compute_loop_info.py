#!/usr/bin/python
"""Sets resources requirements in resources.def
 
   Calculates the run date information, and
   the number of maestro loops required given
   nemo_freq. 

   Information is read from experiment.cfg,
   processed and then written to resources.def

   	   NOT FOR PRODUCTION

   NCS, 05/2017. 
"""
import math
# dictionary mapper to get days in month
days_in_month = {1 : 31,
                 2 : 28,
                 3 : 31,
                 4 : 30,
                 5 : 31,
                 6 : 30,
                 7 : 31,
                 8 : 31,
                 9 : 30,
                 10 : 31,
                 11 : 30,
                 12 : 31
                 }

def calc_chunk_nloops(start, stop, nemo_freq):
    """Compute the number of loops required to
       run from start to end at nemo_freq
       
       Parameters:
       -----------
         start : str
           The start date with format "yyyy:mm"
         end : str
           The end date with format "yyyy:mm"
         nemo_freq : str
           The chunk size that NEMO is run in, with
           format 'NNm' for months or 'nny' for years
           where 'nn' is an integer.
           
       Returns:
       ---------
         nloops : int
           The number of loop iterations required
         nemo_freq_months : int
           The nemo frequency in months
    """

    ss = start.split(':')
    ee = end.split(':')
   
    sy = int(ss[0])
    sm = int(ss[1])
    ey = int(ee[0])
    em = int(ee[1])

    if sm not in range(1,13):
        raise ValueError('Start month must be 1-12')
    if em not in range(1,13):
        raise ValueError('End month must be 1-12')

    print 'Start year, month', sy, sm
    print 'End year, month', ey, em
 
    if 'y' in nemo_freq:
        nf = int(nemo_freq.replace('y',''))*12 
    elif 'm' in nemo_freq:
        nf = int(nemo_freq.replace('m','')) 
    else:
        raise ValueError('nemo_freq must have a format such as 3m or 2y')

    print 'nemo_freq in months:', nf

    # Number of months in the run
    # run_nmonth0 = (sy-1)*12 + sm
    nmonths = (ey - sy)*12 + (em -sm) + 1
    dyears = nmonths / 12.0

    print 'Months in run:', nmonths , '(which is {dyears} years)'.format(dyears=dyears)

    # Number of steps/loops to complete the run
    nloops = nmonths / nf
    print 'Number of loop iterations is:', nloops
    if nmonths%nf != 0:
      print 'Should be exiting here with sensible error'
      raise ValueError('\n\n ERROR: nemo_freq must divide evenly into the number of months in the run')

    return '{:04d}'.format(sy), '{:02d}'.format(sm), nloops, nf


def setdef(ifile, value):
    """Does a sed-like replace of values in file
    """
    import os
    import re
    ifile_bak = ifile + '.bak'
    os.rename(ifile, ifile_bak)
    fin = open(ifile_bak,'r')
    fout = open(ifile,'w')
    
    val = value.split('=')
    vleft = val[0]
    vright = val[1]

    #print vleft, vright
    for line in fin:
        str1=vleft + '=.*\n'
        line = re.sub(str1, value + '\n', line)
        fout.write(line)
        
    fin.close()
    fout.close()

if __name__ == '__main__':
    import subprocess
    import os

    # Read information from experiment.cdf
    if not os.path.exists('resources/resources.def'):
        raise SystemExit('No resources/resources.def file found. Configuration MUST be executed from $SEQ_EXP_HOME)')

#    try:
    if True:
        start = subprocess.check_output(['getdef', 'experiment.cfg', 'start']).strip()
        end = subprocess.check_output(['getdef', 'experiment.cfg', 'end']).strip()
        nemo_freq = subprocess.check_output(['getdef', 'experiment.cfg', 'nemo_freq']).strip()
        nemo_rdt = int(subprocess.check_output(['getdef', 'experiment.cfg', 'nemo_rdt']).strip())
        #backend_machine = subprocess.check_output(['getdef', 'experiment.cfg', 'machine']).strip()
        #frontend_machine = subprocess.check_output(['getdef', 'experiment.cfg', 'machine']).strip()
        hall = subprocess.check_output(['getdef', 'experiment.cfg', 'hall']).strip()
        nemo_wallclock = subprocess.check_output(['getdef', 'experiment.cfg', 'nemo_wallclock']).strip()
  
        # Compute how many loops are required in total
        run_start_year, run_start_month, nloops, nemo_freq_months = calc_chunk_nloops(start, end, nemo_freq)
        
        if hall == "hall1":
            frontend='eccc-ppp1'
            backend='hare'
        elif hall == "hall2":
            frontend='eccc-ppp2'
            backend='brooks'
        else:
            raise ValueError('Hall be must hall1 or hall2 in experiment.cfg')

        modstr = "NEMO_LOOP_END={0}".format(nloops)
        setdef('resources/resources.def', modstr)
        modstr = "NEMO_FREQ_MONTHS={0}".format(nemo_freq_months)
        setdef('resources/resources.def', modstr)
        modstr = "NEMO_RUN_START_YEAR={0}".format(run_start_year)
        setdef('resources/resources.def', modstr)
        modstr = "NEMO_RUN_START_MONTH={0}".format(run_start_month)
        setdef('resources/resources.def', modstr)
        modstr = "NEMO_RDT={0}".format(nemo_rdt)
        setdef('resources/resources.def', modstr)
        modstr = "BACKEND={0}".format(backend)
        setdef('resources/resources.def', modstr)
        modstr = "FRONTEND={0}".format(frontend)
        setdef('resources/resources.def', modstr)
        modstr = "NEMO_WALLCLOCK={0}".format(nemo_wallclock)
        setdef('resources/resources.def', modstr)

#    except:
#        print("FAILURE: Configuration of resources.def \n Make sure that experiment.cfg is valid.")
 
