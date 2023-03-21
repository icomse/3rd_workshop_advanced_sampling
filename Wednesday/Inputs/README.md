# Systems of interest
In this session, we consider two the following two systems.

## NaCl system
The NaCl system is used in Exercises 1, 2, 3, and 5, in which the general goal is to run advanced sampling methods to get the free energy surface of NaCl as a function of the ion-pair distance. In the folder `NaCl`, there are the following files
- `NaCl.gro`: A configuration of NaCl used in Exercise 1, and 2. This file is the same as `NaCl_StartingStructure-1.gro` provided by [this repo](https://github.com/valsson-group/masterclass-22-11/tree/main/SetupSystem).
- `NaCl.top`: A topology file of NaCl used in Exercises 1, 2, 3, and 5. This file is the same as `NaCl.top` used in [this repo](https://github.com/valsson-group/masterclass-22-11/tree/main/SetupSystem) except that in our topology file, the `[ moleculetype ]` and `[ atoms ]` directives are explicitly specified rather than included at the top of the file. This makes it easier to apply a position restraint to the NA atom in Exercises 2.
- `MD-NPT.mdp`: An `mdp` file for standard MD simulation in the NPT ensemble. This file is modified from `MD-NPT.mdp` provided by [this repo](https://github.com/valsson-group/masterclass-22-11/tree/main/SetupSystem), with the only differences being `nsteps` and ` nstxout-compressed`. This file is used in Exercise 1.
- `MD-NVT.mdp`: An `mdp` file for standard MD simulation in the NVT ensemble. This file is modified from `MD-NPT.mdp` by removing the settings for the barostat. This file serves as the base from which the `mdp` files in Exercises 2, 3, and 5 are modified from.
- `configs`: A folder of extra `gro` files of NaCl generated from an NVT MD simulation initialized with `NaCl.gro`, `NaCl.top`, and `MD-NVT.mdp`. Specifically, 8 configurations were extracted from timeframes from 430, 440, 450, ..., 500 ps of the simulation. NaCl_1.gro and NaCl_2.gro are used in Exercise 3 and all 8 configurations are used in Exercise 5. 

Importantly, smaller cutoffs (rlist, rvdw, and rcoloumb) are used in the `mdp` files to allow a smaller water box, hence a shorter time to finish the simulations in the tutorials. 

## 4-site system
The 4-site system is only used in Exercise 4. This is a toy model composed of 4 linearly connected interaction sites. As shown in the presentation, there are two torsional metastable states (cis and trans) for this system due to the diple from the first and last atoms (with charges of +/- 0.2) of the molecule. The files were prepared in the folder `Prep` by the script `prep.sh` (all inputs are provided), with the configuration having the average NPT volume extracted from the final NPT MD simulation as the `sys.gro` used in the tutorial. 

