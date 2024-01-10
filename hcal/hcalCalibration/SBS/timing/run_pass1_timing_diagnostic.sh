#!/bin/sh

# SSeeds - 12.12.23
# Contains all GMn timing sets accounting for timing shifts where they are manifest

root -l -b -q 'timing_diagnostic.C("gmn",4,0,11589,11595,0,0,true)'

wait

root -l -b -q 'timing_diagnostic.C("gmn",4,0,0,0,11589,11595,true)'

wait

root -l -b -q 'timing_diagnostic.C("gmn",7,0,0,0,0,0,true)'

wait

root -l -b -q 'timing_diagnostic.C("gmn",11,1,12314,12995,0,0,true)'

wait

root -l -b -q 'timing_diagnostic.C("gmn",11,1,0,0,12314,12995,true)'

wait

root -l -b -q 'timing_diagnostic.C("gmn",14,1,13239,13260,0,0,true)'

wait

root -l -b -q 'timing_diagnostic.C("gmn",14,1,13261,13407,0,0,true)'

wait

root -l -b -q 'timing_diagnostic.C("gmn",8,1,0,0,0,0,true)'

wait

root -l -b -q 'timing_diagnostic.C("gmn",9,1,13682,13682,0,0,true)'

wait

root -l -b -q 'timing_diagnostic.C("gmn",9,1,0,0,13682,13682,true)'

wait
