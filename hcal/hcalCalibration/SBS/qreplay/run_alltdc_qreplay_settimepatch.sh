#!/bin/sh

# SSeeds - 9.1.23 - cross fingers and hope this doesn't get shut down for no reason

cd /w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/hcalCalibration/SBS/qreplay

#run all tdc
root -l -b -q 'qreplay_settimepatch.C("gmn",11,1,"tdc","../timing/parameters/tdcsetoffsets_gmn_conf11_pass1_12314_to_12995_exclude_0_to_0.txt","../timing/parameters/tdcsetoffsets_gmn_conf11_pass1_12996_to_13063_exclude_0_to_0.txt")'

wait

root -l -b -q 'qreplay_settimepatch.C("gmn",9,1,"tdc","../timing/parameters/tdcsetoffsets_gmn_conf9_pass1_13656_to_13796_exclude_13682_to_13682.txt","../timing/parameters/tdcsetoffsets_gmn_conf9_pass1_13682_to_13682_exclude_0_to_0.txt")'

wait

root -l -b -q 'qreplay_settimepatch.C("gmn",8,1,"tdc","../timing/parameters/tdcsetoffsets_gmn_conf8_pass1_0_to_0_exclude_0_to_0.txt")'

wait

root -l -b -q 'qreplay_settimepatch.C("gmn",14,1,"tdc","../timing/parameters/tdcsetoffsets_gmn_conf14_pass1_0_to_0_exclude_0_to_0.txt")'

wait

root -l -b -q 'qreplay_settimepatch.C("gmn",4,0,"tdc","../timing/parameters/tdcsetoffsets_gmn_conf4_pass0_0_to_0_exclude_0_to_0.txt")'

wait

root -l -b -q 'qreplay_settimepatch.C("gmn",7,0,"tdc","../timing/parameters/tdcsetoffsets_gmn_conf7_pass0_0_to_0_exclude_0_to_0.txt")'

wait


