mpirun -np 1 gmx_mpi editconf -f sys_init.gro -o sys_box.gro -bt cubic -d 0.6
mpirun -np 1 gmx_mpi solvate -cp sys_box.gro -p sys.top -o sys_sol.gro -cs
mpirun -np 1 gmx_mpi grompp -f em.mdp -c sys_sol.gro -p *top -o em.tpr
mpirun -np 1 gmx_mpi mdrun -deffnm em -ntomp 1
mpirun -np 1 gmx_mpi grompp -f nvt_equil.mdp -c em.gro -p sys.top -o nvt_equil.tpr
mpirun -np 1 gmx_mpi mdrun -deffnm nvt_equil -ntomp 1
mpirun -np 1 gmx_mpi grompp -f npt_equil.mdp -c nvt_equil.gro -t nvt_equil.cpt -p sys.top -o npt_equil.tpr -maxwarn 1
mpirun -np 1 gmx_mpi mdrun -deffnm npt_equil -ntomp 1
mpirun -np 1 gmx_mpi grompp -f md.mdp -c npt_equil.gro -p sys.top -t npt_equil.cpt -o md.tpr
mpirun -np 1 gmx_mpi mdrun -deffnm md -ntomp 1
