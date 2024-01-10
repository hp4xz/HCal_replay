#!/bin/sh

# SSeeds - 9.23.23

## Usage
#./run_timewalk.sh <kine> <iter>

exper=$1
conf=$2
iter=$3
pass=$4

cd /w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/hcalCalibration/SBS/timing

root -l -b -q 'timewalk.C("'$exper'",'$conf','$iter','$pass')'
