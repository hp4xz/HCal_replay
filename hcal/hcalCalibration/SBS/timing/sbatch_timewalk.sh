#!/bin/bash

#SBATCH --partition=production
#SBATCH --account=halla
#SBATCH --mem-per-cpu=1500

# if [[ $(type -t module) != function && -r ${MODULES} ]]; then
# source ${MODULES}
# fi

# if [ -d /apps/modulefiles ]; then
# module use /apps/modulefiles
# fi

#module load gcc/9.2.0

#source /site/12gev_phys/softenv.sh 2.6

exper=$1
kine=$2
iter=$3
pass=$4

cd $SWIF_WORK_DIRECTORY

root -b -q 'timewalk.C("'$exper'",'$kine','$iter','$pass')'
