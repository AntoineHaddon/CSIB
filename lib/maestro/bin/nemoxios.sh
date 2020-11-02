#!/bin/bash
# ==============================================================
# Script to partition cores on a same node between XIOS and NEMO
# used with aprun in nemo_run.tsk
# xios_jpnij: total number of cores used for XIOS
# D. Yang, 2/Nov/2020
# ==============================================================

xios_jpnij=$1
xios_jpnijm1=$(($xios_jpnij - 1))
if [ $ALPS_APP_PE -le $xios_jpnijm1 ] ; then
  exec ./xios_server.exe
else
  exec ./nemo.exe
fi
exit 0
