#####################################################
# Author : Simona Flavoni for NEMO
# Contact : sflod@locean-ipsl.upmc.fr
#
# Some scripts called by sette.sh 
# fcm_job.sh   : simple job to run NEMO with fcm 
######################################################
#set -x
set -o posix
#set -u
#set -e
#+
#
# ================
# fcm_job.sh
# ================
#
# --------------------------
# Simple job for NEMO tests 
# --------------------------
#
# SYNOPSIS
# ========
#
# ::
#
#  $ ./fcm_job.sh INPUT_FILE_CONFIG_NAME NUMBER_PROC TEST_NAME MPI_INTERACT MPI_FLAG
#
#
# DESCRIPTION
# ===========
#
# Simple job for SET TESTS for NEMO (SETTE)
# 
#   get input files (if needed) : tar file  
#  (note this job needs to have an input_CONFIG.cfg in which can be found input tar file name)
#
#   runs job in interactive or batch mode : all jobs using 1 process are run interactive, and all MPP jobs are
#
#   run in batch (MPI_INTERACT="no") or interactive (MPI_INTERACT="yes") see sette.sh and BATCH_TEMPLATE directory
#
#   and call post_test_tidyup function (that moves in NEMO_VALIDATION_DIR solver.stat, tracer.stat (for LOBSTER & PISCES) & ocean.output)
#
# EXAMPLES
# ========
#
# ::
#
#  $ ./fcm_job.sh INPUT_FILE_CONFIG_NAME NUMBER_PROC TEST_NAME MPI_INTERACT MPI_FLAG
#
#  run a job of config GYRE with 1 processor SHORT test ( 5 days ) using an interactive run without mpirun
#  $ ./fcm_job.sh input_GYRE.cfg 1 SHORT yes no
#
#  run a job of config ORCA2_LIM_PISCES	with 8 processors test RESTARTABILITY submitting the job to the batch queue system and using mpirun
#  $ ./fcm_job.sh input_ORCA2_LIM_PISCES.cfg 8 LONG no yes
#
#
# TODO
# ====
#
# option debug
#
#
# EVOLUTIONS
# ==========
#
# $Id: fcm_job.sh 3294 2012-01-28 16:44:18Z rblod $
#
#
#
#   * creation
#
#-
#

usage=" Usage : ./fcm_job.sh input_CONFIG_NAME.cfg  NUMBER_OF_PROCS TEST_NAME INTERACT MPI_FLAG"
usage=" example : ./fcm_job.sh input_ORCA2_LIM_PISCES.cfg 8 SHORT no/yes no/yes"


minargcount=2
        if [ ${#} -lt ${minargcount} ]
        then
                echo "not enought arguments for fcm_job.sh script"
                echo "control number of argument of fcm_job.sh in sette.sh"
                echo "${usage}"
        exit 1
        fi
        unset minargcount
	if [ ! -f ${SETTE_DIR}/output.sette ] ; then
	        touch ${SETTE_DIR}/output.sette
	fi
       

export NB_PROCS=$1
export JOB_FILE=$2
################################################################
# RUN OPA
cd ${EXE_DIR}

# submit job to batch system 
        if [ ${NB_PROC} -eq 1 ]; then
	   eval ${BATCH_COMMAND_SEQ} $JOB_FILE 
        else
	   eval ${BATCH_COMMAND_PAR} $JOB_FILE
        fi

