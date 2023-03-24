#!/bin/bash

ARGC=$#
NARGS=12

if [ $ARGC -ne $NARGS ]; then
    echo
    echo "Usage: extend.sh [epath] [grofilepath] [rstfilepath] \ "
    echo "       [nsteps] [moveid] [xtcfilepath] [opfilepath] [tprpath] [oldtprpath] [cptpath] [edrpath] "
    echo
    return
fi

epath=${1}
grofilepath=${2}
rstfilepath=${3}
nsteps=${4}
moveid=${5}
xtcfilepath=${6}
opfilepath=${7}
tprpath=${8}
oldtprpath=${9}
cptpath=${10}
edrpath=${11}
trrpath=${12}

cd $epath

cp ${rstfilepath} rst/active_in.cpt

export OMP_NUM_THREADS=1

# Run simulation
gmx convert-tpr -s ${oldtprpath} -extend 0.1 -o ${tprpath}
gmx mdrun -s ${tprpath} -ntmpi 12 -ntomp 1 -cpi rst/active_in.cpt -noappend -deffnm sim
mv sim*edr log/${moveid}.edr
mv sim*xtc ${xtcfilepath}
mv sim*trr ${trrpath}
mv sim*cpt rst/${moveid}.cpt
rm -f mdout.mdp sim.*
rm -f rst/active_in.cpt

# Analyze simulation

python /ocean/projects/see220002p/$(whoami)/LiF-RETIS/op.py -g ${grofilepath} -t ${trrpath} -o ${opfilepath}

