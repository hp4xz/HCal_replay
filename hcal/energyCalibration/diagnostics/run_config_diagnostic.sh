#!/bin/sh

# SSeeds - 3.3.23 - simple shell script to run ecal_diagnostic.C for use with batch farm scripts

## Usage
#./run_config_diagnostic <kine> <target> <field> <iter>

kine=$1
tar=$2
field=$3
iter=$4

cd /w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/energyCalibration/diagnostics

root -l -b -q 'ecal_config_diagnostic.C('$kine',"'$tar'",'$field','$iter')'

