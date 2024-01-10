#!/bin/sh

# SSeeds -12.12.23

exp=$1
kine=$2
pass=$3
alltarg=$4

cd /w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/hcalCalibration/SBS/timing

root -l -b -q 'timing_diagnostic.C("'$exp'",'$kine','$pass','$alltarg')'
