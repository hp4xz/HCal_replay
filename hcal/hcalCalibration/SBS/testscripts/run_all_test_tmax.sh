#!/bin/sh

# SSeeds - 9.23.23 - run all timing calibration sets test_tmax.C

echo "Running tmax.C on SBS4 set 1"

root -l -b -q 'test_tmax.C()'

wait

echo "Running test_tmax.C on SBS4"

root -l -b -q 'test_tmax.C(11587,15000)'

wait

echo "Running test_tmax.C on SBS7"

root -l -b -q 'test_tmax.C(12049,15000)'

wait

echo "Running test_tmax.C on SBS11"

root -l -b -q 'test_tmax.C(12369,15000)'

wait

echo "Running test_tmax.C on SBS14"

root -l -b -q 'test_tmax.C(13314,15000)'

wait

echo "Running test_tmax.C on SBS8"

root -l -b -q 'test_tmax.C(13486,15000)'

wait

echo "Running test_tmax.C on SBS9"

root -l -b -q 'test_tmax.C(13662,15000)'

echo "Finished!"
