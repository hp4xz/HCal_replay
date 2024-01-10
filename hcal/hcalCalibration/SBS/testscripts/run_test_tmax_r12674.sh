#!/bin/sh

# SSeeds - 9.23.23 - run all timing calibration sets test_tmax.C

echo "Running tmax.C on SBS11, tmax=1000"

root -l -b -q 'test_tmax.C(12674,10009)'

wait

echo "Running tmax.C on SBS11, tmax=30"

root -l -b -q 'test_tmax.C(12674,10010)'

wait

echo "Running tmax.C on SBS11, tmax=9.0"

root -l -b -q 'test_tmax.C(12674,10011)'

wait

echo "Finished!"
