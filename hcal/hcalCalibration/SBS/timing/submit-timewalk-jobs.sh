#!/bin/bash

# swif2 job submit for timewalk.C for hcal calibration
# note that diagnostics should be done first in order to have good elastic cuts by kine/target/field-setting located in configuration files
#Usage
#./submit-timewalk-jobs.sh <string: experiment> <int: configuration> <bool: quasi-replay> <int: replay-pass>

exper=$1
conf=$2
iter=$3
pass=$4

work_flow='timewalk_'$exper'_sbs'$conf'_i'$iter'_pass'$pass

swif2 create $work_flow

jobname=$conf-$iter

script='/w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/hcalCalibration/SBS/timing/run_timewalk.sh'

cd $SWIF_WORK_DIRECTORY

swif2 add-job -workflow $work_flow -partition production -name timewalk-$jobname -cores 1 -disk 25GB -ram 1500MB $script $exper $conf $iter $pass

# run the workflow and then print status
swif2 run $work_flow
echo -e "\n Getting workflow status.. \n"
swif2 status $work_flow
