#!/bin/sh

# SSeeds - 9.23.23 - run all timing calibration sets tmax.C

maxrun=$1
nsig=$2

echo "Setting up tmax run all with a " $maxrun " run maximum"

echo "Running tmax.C on SBS4 set 1"

root -l -b -q 'tmax.C("gmn",4,0,1,'$maxrun',0,0,6,'$nsig')'

wait

echo "Running tmax.C on SBS7 set 1"

root -l -b -q 'tmax.C("gmn",7,0,1,'$maxrun',0,0,6,'$nsig')'

wait

echo "Running tmax.C on SBS14 set 1"

root -l -b -q 'tmax.C("gmn",14,1,1,'$maxrun',0,0,6,'$nsig')'

wait

echo "Running tmax.C on SBS14 set 2"

root -l -b -q 'tmax.C("gmn",14,1,2,'$maxrun',0,0,6,'$nsig')'

wait

echo "Running tmax.C on SBS11 set 1"

root -l -b -q 'tmax.C("gmn",11,1,1,'$maxrun',0,0,6,'$nsig')'

wait

echo "Running tmax.C on SBS11 set 2"

root -l -b -q 'tmax.C("gmn",11,1,2,'$maxrun',0,0,6,'$nsig')'

wait

echo "Running tmax.C on SBS8 set 1"

root -l -b -q 'tmax.C("gmn",8,1,1,'$maxrun',0,0,6,'$nsig')'

wait

echo "Running tmax.C on SBS9 set 1"

root -l -b -q 'tmax.C("gmn",9,1,1,'$maxrun',0,0,6,'$nsig')'

echo "Finished!"
