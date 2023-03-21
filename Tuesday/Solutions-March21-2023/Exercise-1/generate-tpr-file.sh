GMX_BIN=gmx
StartingGeometry=${1}
RunFilename=NaCl_NPT-300K
${GMX_BIN}  grompp -f MD-NPT.mdp -c ${StartingGeometry} -p NaCl.top -o ${RunFilename}.tpr  -maxwarn 1
