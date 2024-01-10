#!/bin/sh

# SSeeds - 12.12.23
# Contains all GMn timing sets accounting for timing shifts where they are manifest

root -l -b -q 'simple_setalign.C("gmn",4,0,11589,11595,0,0,false,0,0,-249.448)'

wait

root -l -b -q 'simple_setalign.C("gmn",4,0,0,0,11589,11595,false,0,0,-249.948)'

wait

root -l -b -q 'simple_setalign.C("gmn",7,0,0,0,0,0,false,0,0,-239.234)'

wait

root -l -b -q 'simple_setalign.C("gmn",11,1,12314,12995,0,0,false,0,0,-271.841)'

wait

root -l -b -q 'simple_setalign.C("gmn",11,1,0,0,12314,12995,false,0,0,-267.471)'

wait

root -l -b -q 'simple_setalign.C("gmn",14,1,13239,13260,0,0,false,0,0,-243.765)'

wait

root -l -b -q 'simple_setalign.C("gmn",14,1,13261,13407,0,0,false,0,0,-243.765)'

wait

root -l -b -q 'simple_setalign.C("gmn",8,1,0,0,0,0,false,0,0,-253.499)'

wait

root -l -b -q 'simple_setalign.C("gmn",9,1,13682,13682,0,0,false,0,0,-250.951)'

wait

root -l -b -q 'simple_setalign.C("gmn",9,1,0,0,13682,13682,false,0,0,-250.951)'

wait
