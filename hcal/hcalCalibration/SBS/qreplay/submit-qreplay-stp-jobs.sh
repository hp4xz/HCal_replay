#!/bin/bash

# swif2 job submit for timewalk.C for hcal calibration
# note that diagnostics should be done first in order to have good elastic cuts by kine/target/field-setting located in configuration files
#Usage
#./submit-timewalk-jobs.sh <string: experiment> <int: configuration> <bool: quasi-replay> <int: replay-pass>

exper=$1
config=$2
pass=$3
type=$4
pp1=$5
pp2=$6
pp3=$7
pp4=$8

work_flow='timepatch_pass2'

swif2 create $work_flow

jobname=$exper'-'$config'-'$type'-'$pp1'-'$pp2'-'$pp3'-'$pp4

script='/w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/hcalCalibration/SBS/qreplay/run_qreplay_settimepatch.sh'

#cd $SWIF_WORK_DIRECTORY
cd /w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/hcalCalibration/SBS/qreplay

echo -e "\n Submitting " $script $exper $config $pass $type $pp1 $pp2 $pp3 $pp4 "\n"

swif2 add-job -workflow $work_flow -partition production -name timepatch-$jobname -cores 1 -disk 5GB -ram 1500MB $script $exper $config $pass $type $pp1 $pp2 $pp3 $pp4

# run the workflow and then print status
swif2 run $work_flow
echo -e "\n Getting workflow status.. \n"
swif2 status $work_flow
