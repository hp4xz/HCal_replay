#!/bin/bash

#sseeds - 12.12.23

exper=$1
kine=$2
pass=$3
all_targets=$4

work_flow='timing_diagnostic'

swif2 create $work_flow

jobname=$exper'-'$kine'-'$pass'-'$all_targets

script='/w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/hcalCalibration/SBS/timing/run_timing_diagnostic.sh'

#cd $SWIF_WORK_DIRECTORY
cd /w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/hcalCalibration/SBS/timing

echo -e "\n Submitting " $script $exper $kine $pass $all_targets "\n"

swif2 add-job -workflow $work_flow -partition production -name timingDiagnostic-$jobname -cores 1 -disk 20GB -ram 2500MB $script $exper $kine $pass $all_targets

# run the workflow and then print status
swif2 run $work_flow
echo -e "\n Getting workflow status.. \n"
swif2 status $work_flow
