#!/bin/sh

# SSeeds - 9.1.23 qreplay batch farm submission script. path-to-n-patch should contain offset data in database format. See qreplay_settimepatch.C for more details.

## Usage
#./run_qreplay.sh <experiment> <configuration> <replay-pass> <type(tdc/adct)> <path-to-first-patch> <path-to-second-patch> <path-to-third-patch> <path-to-fourth-patch>
## for last four params, pass 0 if not used.

exper=$1
config=$2
pass=$3
type=$4
pp1=$5
pp2=$6
pp3=$7
pp4=$8

cd /w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/hcalCalibration/SBS/qreplay

root -l -b -q 'qreplay_settimepatch.C("'$exper'",'$config','$pass',"'$type'","'$pp1'","'$pp2'","'$pp3'","'$pp4'")'
