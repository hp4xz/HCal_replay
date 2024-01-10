#!/bin/sh

# SSeeds - 12.12.23
# Contains all GMn timing sets accounting for timing shifts where they are manifest.

root -l -b -q 'qreplay_settimepatch.C("gmn",4,0,false,"adct","../timing/parameters/simple_adctoffsets_gmn_conf4_pass1_0_to_0_exclude_11589_to_11595.txt","../timing/parameters/simple_adctoffsets_gmn_conf4_pass1_11589_to_11595_exclude_0_to_0.txt")'

wait

root -l -b -q 'qreplay_settimepatch.C("gmn",7,0,false,"adct","../timing/parameters/simple_adctoffsets_gmn_conf7_pass1_0_to_0_exclude_0_to_0.txt")'

wait

root -l -b -q 'qreplay_settimepatch.C("gmn",11,1,false,"adct","../timing/parameters/simple_adctoffsets_gmn_conf11_pass1_0_to_0_exclude_12314_to_12995.txt","../timing/parameters/simple_adctoffsets_gmn_conf11_pass1_12314_to_12995_exclude_0_to_0.txt")'

wait

root -l -b -q 'qreplay_settimepatch.C("gmn",14,1,false,"adct","../timing/parameters/simple_adctoffsets_gmn_conf14_pass1_13239_to_13260_exclude_0_to_0.txt","../timing/parameters/simple_adctoffsets_gmn_conf14_pass1_13261_to_13407_exclude_0_to_0.txt")'

wait

root -l -b -q 'qreplay_settimepatch.C("gmn",8,1,false,"adct","../timing/parameters/simple_adctoffsets_gmn_conf8_pass1_0_to_0_exclude_0_to_0.txt")'

wait

root -l -b -q 'qreplay_settimepatch.C("gmn",9,1,false,"adct","../timing/parameters/simple_adctoffsets_gmn_conf9_pass1_0_to_0_exclude_13682_to_13682.txt","../timing/parameters/simple_adctoffsets_gmn_conf9_pass1_13682_to_13682_exclude_0_to_0.txt")'

wait
