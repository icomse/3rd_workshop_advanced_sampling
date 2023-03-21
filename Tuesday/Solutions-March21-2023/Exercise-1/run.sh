GMX_BIN="gmx"
GMX_MDRUN_BIN="mpirun -np 4 mdrun_mpi"
RunFilename=NaCl_NPT-300K

${GMX_MDRUN_BIN}  -deffnm ${RunFilename}   -plumed plumed.dat 
echo "0" | ${GMX_BIN} trjconv -f ${RunFilename}.xtc -s ${RunFilename}.tpr -pbc whole -o ${RunFilename}.pbc-whole.xtc


