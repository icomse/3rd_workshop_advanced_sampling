#!/bin/bash

ARGC=$#
NARGS=14

if [ $ARGC -ne $NARGS ]; then
    echo
    echo "Usage: simulate.sh [epath] [grofilepath] [inputmasterpath] [randn1] \ "
    echo "       [nsteps] [moveid] [inputfilepath] [xtcfilepath] [opfilepath] [topfilepath] [tprpath] [edrpath] [cptpath]"
    echo
    return
fi

epath=${1}
grofilepath=${2}
inputmasterpath=${3}
randn1=${4}
nsteps=${5}
moveid=${6}
inputfilepath=${7}
xtcfilepath=${8}
opfilepath=${9}
topfilepath=${10}
tprpath=${11}
edrpath=${12}
cptpath=${13}
trrpath=${14}


cd $epath

# Generate input file
sed "s/RRRR/${randn1}/" ${inputmasterpath} | sed "s/NNNN/${nsteps}/" > ${inputfilepath}

export OMP_NUM_THREADS=1

# Run simulation
gmx grompp -f ${inputfilepath} -c ${grofilepath} -p ${topfilepath} -o ${tprpath}
gmx mdrun -ntmpi 12 -ntomp 1 -s ${tprpath} -noappend -deffnm sim
mv sim*edr ${edrpath}
mv sim*xtc ${xtcfilepath}
mv sim*cpt ${cptpath}
mv sim*trr ${trrpath}
rm -f mdout.mdp sim.*

# Analyze simulation

python /ocean/projects/see220002p/$(whoami)/LiF-RETIS/op.py -g ${grofilepath} -t ${trrpath} -o ${opfilepath}

