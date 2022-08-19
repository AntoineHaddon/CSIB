#!/usr/bin/python
""" Compute essential counter information within an maestro loop.

    Computes information such as the chunk start and stop dates,
    the NEMO first and last iterations. This information is written
    to the file nemo_counter_info.cfg, which can then be sourced.

              NOT FOR PROD

    NCS, 05/2017
"""

# Look at how loop counter will work
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

def calc_nemo_chunk_dates(run_start_year, run_start_month, loop, nemo_freq_months):
    """ Compute the start and end date of the current chunk
    
    Parameters:
    ------------
         run_start_year : int
           The start year (of the run) 
         run_start_month : int
           The start month (of the run) 
         loop: int
           The current iteration of the loop counter
         nemo_freq_months : int
           The chunk size that NEMO is run in units of months.
           
       Returns:
       ---------
         chunk_start_day : str
           The start day of the chunk
         chunk_start_month : str
           The start month of the chunk
         chunk_start_year : str
           The start year of the chunk
         chunk_end_day : str
           The end day of the chunk
         chunk_end_month : str
           The end month of the chunk
         chunk_end_year : str
           The end year of the chunk
    """
    nf = int(nemo_freq_months)

    if run_start_month not in range(1,13):
        raise ValueError('Start month must be 1-12')
 
    sm = run_start_month
    sy = run_start_year

    #print "inputs"
    #print (run_start_year, run_start_month, loop, nemo_freq_months)
 
    ll = loop -1
  
    cl_start_nmonth = sm + ll*nf
    cl_end_nmonth = sm + (ll+1)*nf - 1

    # Get the calendar start/end month for each loop segment
    cl_start_cal_month = (cl_start_nmonth)%12 if cl_start_nmonth%12 > 0 else 12
    cl_end_cal_month =  cl_end_nmonth%12 if cl_end_nmonth%12 > 0 else 12
    # Get the calendar start/end year for each loop segment
    cl_start_cal_year = int((sm + ll*nf-1)/12 + sy)
    cl_end_cal_year = int(math.ceil((sm + (ll+1)*nf -1)/12.0)) + sy -1

    cl_start_cal_day = 0o1
    cl_end_cal_day = days_in_month[cl_end_cal_month]
    
    # List of all years in this chunk
    chunk_years = range(cl_start_cal_year, cl_end_cal_year+1,1)
    chunk_years_str = " ".join('%04d' % year for year in chunk_years)
    #print 'chunk_years ', chunk_years_str

    #print "outputs"
    #print cl_start_cal_day, cl_start_cal_month, cl_start_cal_year, cl_end_cal_day, cl_end_cal_month, cl_end_cal_year

    return ('%02d' % (cl_start_cal_day),'%02d' % (cl_start_cal_month), '%04d' % (cl_start_cal_year), 
           '% 02d' % (cl_end_cal_day), '%02d' % (cl_end_cal_month), '%04d' % (cl_end_cal_year), chunk_years_str)

def calc_nemo_iters(chunk_start_month, chunk_start_year, chunk_end_month, chunk_end_year, nemo_freq_months, nemo_rdt):
    """Calculate the number of nemo iterations/steps in this chunk
    
       This is used to set the end iteration in the namelist.
    """

    sec_per_day = 24 * 60 * 60

    chunk_start_month = int(chunk_start_month)
    chunk_end_month = int(chunk_end_month)
    chunk_start_year = int(chunk_start_year)
    chunk_end_year = int(chunk_end_year)

    # Figure out how many days are in each loop chunk
    months_in_chunk = [chunk_start_month + i%12  if chunk_start_month + i%12 <= 12 else (chunk_start_month + i)%12for i in range(nemo_freq_months)]
    days_permonth_in_chunk = [ days_in_month[m] for m in months_in_chunk ]
    total_days_in_chunk = sum(days_permonth_in_chunk)
    
    nemo_n_iters = total_days_in_chunk * sec_per_day / nemo_rdt

    # Compute how many days have been run up until now.
    print ( [1 + i%12 for i in range(chunk_start_month-1)] )
    months_before_chunk = [1 + i%12  for i in range(chunk_start_month-1)]
    print ( months_before_chunk )
    days_permonth_before_chunk = [ days_in_month[m] for m in months_before_chunk ]
    total_days_before_chunk = sum(days_permonth_before_chunk) + 365*(chunk_start_year-1)
    
    nemo_nn_it000 = total_days_before_chunk * sec_per_day / nemo_rdt + 1
    nemo_nn_itend = nemo_nn_it000 + nemo_n_iters -1
 
    #print total_days_in_chunk, total_days_in_chunk/365. 
    #print "number of nemo iterations:", nemo_n_iters
    return nemo_n_iters, total_days_in_chunk, nemo_nn_it000, nemo_nn_itend

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
    import sys

    # The loop index number should be fed in as the only command line arg
    if len(sys.argv) != 2:
        raise SystemExit("compute_iteration_info FAILURE: must pass exactly 1 command line arg (loop index)")
    
    loop_index = int(sys.argv[1])
    #print loop_index

    if type(loop_index) != int:
        raise TypeError("compute_iteration_info FAILURE: Loop index must be an integer")

    # Get some run information from resources.def
    if not "SEQ_EXP_HOME" in os.environ:
        raise SystemExit("compute_iteration_info FAILURE:: SEQ_EXP_HOME must be defined")     

#    try:
    if True:
        run_start_year = subprocess.Popen(['getdef', 'resources', 'NEMO_RUN_START_YEAR'],stdout=subprocess.PIPE).communicate()[0].strip()
        run_start_month = subprocess.Popen(['getdef', 'resources','NEMO_RUN_START_MONTH'],stdout=subprocess.PIPE).communicate()[0].strip()
        nemo_freq_months = subprocess.Popen(['getdef', 'resources','NEMO_FREQ_MONTHS'],stdout=subprocess.PIPE).communicate()[0].strip()
        nemo_rdt = subprocess.Popen(['getdef', 'resources', 'NEMO_RDT'],stdout=subprocess.PIPE).communicate()[0].strip()
#    except:
#        raise RuntimeError("Could not get defitions from $SEQ_EXP_HOME/resources/resources.def")


    # start / end dates of this loop
    chunk_start_day, chunk_start_month, chunk_start_year, chunk_end_day,chunk_end_month, chunk_end_year, chunk_years_str = calc_nemo_chunk_dates(int(run_start_year), 
                                int(run_start_month), loop_index, int(nemo_freq_months))
    
    # Nemo iterations this loop
    nemo_n_iters, days_in_chunk, nemo_nn_it000, nemo_nn_itend = calc_nemo_iters(chunk_start_month, chunk_start_year, chunk_end_month, chunk_end_year, int(nemo_freq_months), int(nemo_rdt))


    # Write to nemo_counter_info.cfg
    with open('nemo_counter_info.cfg', 'w') as ff:
        ff.write('NEMO_CHUNK_START_DAY=%s\n' % (chunk_start_day))
        ff.write('NEMO_CHUNK_START_MONTH=%s\n' % (chunk_start_month))
        ff.write('NEMO_CHUNK_START_YEAR=%s\n' % (chunk_start_year))
        ff.write('NEMO_CHUNK_END_DAY=%s\n' % (chunk_end_day.replace(" ","")))
        ff.write('NEMO_CHUNK_END_MONTH=%s\n' % (chunk_end_month))
        ff.write('NEMO_CHUNK_END_YEAR=%s\n' % (chunk_end_year))
        ff.write('NEMO_CHUNK_START_DATE=%s\n' % ( chunk_start_year + chunk_start_month + chunk_start_day))
        ff.write('NEMO_CHUNK_END_DATE=%s\n' % ( (chunk_end_year + chunk_end_month + chunk_end_day).replace(" ","")))
        ff.write('NEMO_CHUNK_NITERS=%08d \n' % (nemo_n_iters))
        ff.write('NEMO_CHUNK_NN_IT000=%08d \n' % (nemo_nn_it000))
        ff.write('NEMO_CHUNK_NN_ITEND=%08d \n' % (nemo_nn_itend))
        ff.write("NEMO_CHUNK_YEARS='%s'\n" % (chunk_years_str))
        



