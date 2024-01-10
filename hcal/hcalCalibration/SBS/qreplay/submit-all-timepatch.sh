#!/bin/sh

# SSeeds - shell script to submit all qreplay_settimepatch.C jobs to farm for adct

cd /w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/hcalCalibration/SBS/qreplay

echo -e "\n Submitting GMn SBS4 with one patch to batch farm..\n"
./submit-qreplay-stp-jobs.sh gmn 4 0 adct ../timing/parameters/adctsetoffsets_gmn_conf4_pass0_0_to_0_exclude_0_to_0.txt 0 0 0

wait

echo -e "\n Submitting GMn SBS7 with one patch to batch farm..\n"
./submit-qreplay-stp-jobs.sh gmn 7 0 adct ../timing/parameters/adctsetoffsets_gmn_conf7_pass0_0_to_0_exclude_0_to_0.txt 0 0 0

wait

echo -e "\n Submitting GMn SBS11 with two patches to batch farm..\n"
./submit-qreplay-stp-jobs.sh gmn 11 1 adct ../timing/parameters/adctsetoffsets_gmn_conf11_pass1_0_to_0_exclude_12450_to_12860.txt ../timing/parameters/adctsetoffsets_gmn_conf11_pass1_12450_to_12860_exclude_0_to_0.txt 0 0

wait

echo -e "\n Submitting GMn SBS14 with two patches to batch farm..\n"
./submit-qreplay-stp-jobs.sh gmn 14 1 adct ../timing/parameters/adctsetoffsets_gmn_conf14_pass1_13239_to_13260_exclude_0_to_0.txt ../timing/parameters/adctsetoffsets_gmn_conf14_pass1_13261_to_13407_exclude_0_to_0.txt 0 0

wait

echo -e "\n Submitting GMn SBS8 with one patch to batch farm..\n"
./submit-qreplay-stp-jobs.sh gmn 8 1 adct ../timing/parameters/adctsetoffsets_gmn_conf8_pass1_0_to_0_exclude_0_to_0.txt 0 0 0

wait

echo -e "\n Submitting GMn SBS9 with one patch to batch farm..\n"
./submit-qreplay-stp-jobs.sh gmn 9 1 adct ../timing/parameters/adctsetoffsets_gmn_conf9_pass1_0_to_0_exclude_0_to_0.txt 0 0 0


