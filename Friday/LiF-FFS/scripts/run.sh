#!/bin/bash
#****************************************************************
#*                             FFS                              *
#*                  Single Run with OP checking                 *
#*                                                              *
#*                     Author: Naomi Trampe                     *
#*                          SAMPEL Lab                          *
#*                   Last update: 03/20/2023                    *
#****************************************************************
ARGC=$#
NARGS=12

if [ $ARGC -ne $NARGS ] 
then
    echo
    echo "Usage: simulate.sh [programpath] [resultpath] [grofilename] [interface] \ "
    echo "       [traj_num] [basinA] [lambda] [basinB] [time] [steps] [exploringscouts (y/n)] [grow]"
    echo
    exit 0
fi

#if necessary, load gromacs here


ROOT=${1}
MPATH=${2}
FILENAME=${3}
INT=${4}
NUM=${5}
A=${6}
l=${7}
B=${8}
T_RUN=${9}
STEPS=${10}
if [ ${11} == n ]
then
  SIM=simulations/int_${INT}/sim_${NUM}/
elif [ ${11} == y ]
then
  SIM=simulations/int_${INT}/ex-sim_${NUM}/
fi
GROW=${12}

    


j=1
cd ${MPATH}${SIM}
sed 's/STEPS/'${STEPS}'/g' ${MPATH}data/mdp.mdp > t_run.mdp
#Run first timestep
gmx grompp -f t_run.mdp -c ${FILENAME}.gro -n ${MPATH}data/index.ndx -p ${MPATH}data/LiF.top -o ${FILENAME}.tpr -maxwarn 1
gmx mdrun -deffnm ${FILENAME} -noappend -ntmpi 1 -ntomp 1
#Analyze first timestep
FILEPATH=${MPATH}${SIM}${FILENAME}.part000${j}
SUCCESS=`python3 ${ROOT}scripts/op.py "run" "F" "Li" ${FILEPATH}".gro" ${FILEPATH}".xtc" ${A} ${l} ${GROW}`
echo "Result of first timeperiod: "${SUCCESS}
#Repeat by extending the tpr, running, and analyzing
#while (( $(bc <<<"$SUCCESS == 0 && $j < 5") ))
while [ "${SUCCESS}" -eq 0 ] && [ "${j}" -lt 3 ]
do
  echo "Starting again..."
  j=$(($j+1))
  gmx convert-tpr -s ${FILENAME}.tpr -extend ${T_RUN} -o ${FILENAME}.tpr
  gmx mdrun -deffnm ${FILENAME} -cpi ${FILENAME}.cpt -noappend -ntmpi 1 -ntomp 1
  if [ ${j} -lt 10 ]
  then
    FILEPATH=${MPATH}${SIM}${FILENAME}.part000${j}
  elif [ ${j} -lt 100 ]
  then
    FILEPATH=${MPATH}${SIM}${FILENAME}.part00${j}
  elif [ ${j} -lt 1000 ]
  then
    FILEPATH=${MPATH}${SIM}${FILENAME}.part0${j}
  else
    FILEPATH=${MPATH}${SIM}${FILENAME}.part${j}
  fi
  SUCCESS=`python3 ${ROOT}scripts/op.py "run" "F" "Li" ${FILEPATH}".gro" ${FILEPATH}".xtc" ${A} ${l} ${GROW}`
done
END=`python3 ${ROOT}scripts/op.py "end" "F" "Li" ${FILEPATH}".gro" ${FILEPATH}".xtc" ${A} ${l} ${GROW}`
echo "last file ends at: "${END}
#move configuration to the next folder
gmx trjconv -f ${FILEPATH}.xtc -e ${END} -o ${FILEPATH}.xtc
gmx trjcat -f *.xtc -o ${FILENAME}.xtc
gmx eneconv -f *.edr -o ${FILENAME}.edr
echo ${SUCCESS}
#Run next block if not an exploring scout
if [ ${11} == n ]
then
  echo "about to write"
  #Write output config
  if [ ${SUCCESS} -eq 2 ]
  then
    python3 ${ROOT}scripts/op.py "write" "F" "Li" ${FILEPATH}".gro" ${FILEPATH}".xtc" ${l} ${MPATH}"simulations/int_"$(($INT+1))"/"${FILENAME}".gro" ${GROW}
  fi
fi
#clean files
rm ${FILENAME}.part*
rm *.tpr
rm *.cpt
rm *.mdp
rm .${FILENAME}*
rm \#${FILENAME}*
