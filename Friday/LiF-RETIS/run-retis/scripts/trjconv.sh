#!/bin/bash

ARGC=$#
NARGS=5

if [ $ARGC -ne $NARGS ]; then
    echo
    echo "Usage: trjconv.sh [epath] [grofilepath] [harvestfilepath] [harvesttime] [outfilepath]"
    echo "This must be run from the main level of each interface ensemble"
    echo
    return
fi

epath=$1
grofilepath=$2
harvestfilepath=$3
harvesttime=$4
outfilepath=$5

cd $epath

# Harvest config

export OMP_NUM_THREADS=1
gmx trjconv -s ${grofilepath} -f ${harvestfilepath} -vel -dump ${harvesttime} -o ${outfilepath} <<< "System"


# Run a gmx grompp and gmx mdrun on the npt-trjconv.mdp (which has 0 step) to generate a new gro with randomized velocities
export OMP_NUM_THREADS=1
gmx grompp -f /ocean/projects/see220002p/$(whoami)/LiF-RETIS/run-retis/masters/npt-trjconv.mdp \
    -c ${outfilepath} \
    -p /ocean/projects/see220002p/$(whoami)/LiF-RETIS/run-retis/masters/LiF.top \
    -maxwarn 1 \
    -o tpr-dummy.tpr
gmx mdrun -deffnm dummy -ntmpi 12 -ntomp 1 -s tpr-dummy.tpr -c dummy.gro
mv dummy.gro ${outfilepath} 
rm -f tpr-dummy.tpr mdout.mdp dummy*

oldgro=`sed -e 's/\.gro/_old.gro/g' <<<"$outfilepath"`
cp $outfilepath $oldgro


