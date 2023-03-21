# Systems of interest
In this session, we consider two systems: NaCl, and a 4-site system.

## NaCl system
The NaCl system is used in Exercises 1, 2, 3, and 5, in which the general goal is to run advanced sampling methods to get the free energy surface of NaCl as a function of the ion-pair distance. In the folder `NaCl`, there are the following files/folders:
- `NaCl.gro`: A configuration of NaCl used in Exercise 1 and 2. This file is the same as `NaCl_StartingStructure-1.gro` provided by [this repo](https://github.com/valsson-group/masterclass-22-11/tree/main/SetupSystem).
- `NaCl.top`: A topology file of NaCl used in Exercises 1, 2, 3, and 5. This file is the same as `NaCl.top` used in [this repo](https://github.com/valsson-group/masterclass-22-11/tree/main/SetupSystem) except that in our topology file, the `[ moleculetype ]` and `[ atoms ]` directives are explicitly specified rather than included at the top of the file. This makes it easier to add a position restraint to the NA atom in Exercises 2.
- `MD-NPT.mdp`: An `mdp` file for standard MD simulation in the NPT ensemble. This file is modified from `MD-NPT.mdp` provided by [this repo](https://github.com/valsson-group/masterclass-22-11/tree/main/SetupSystem), with the only differences being `nsteps` and `nstxout-compressed`. This file is only used in Exercise 1.
- `MD-NVT.mdp`: An `mdp` file for standard MD simulation in the NVT ensemble. This file is modified from `MD-NPT.mdp` by removing the settings for the barostat. This file serves as the base from which the `mdp` files in Exercises 2, 3, and 5 are modified.
- `configs`: A folder of extra `gro` files of NaCl generated from a pulling simulation as the one in Exercise 2. Specifically, 8 configurations with an ion-pair distance between 0.27 to 0.30 nm were extracted from the simulation. `NaCl_0.gro`, `NaCl_1.gro`, `NaCl_2.gro`, and `NaCl_3.gro` are used in Exercise 3 and all 8 configurations are used in Exercise 5. 

Importantly, smaller cutoffs (`rlist`, `rvdw`, and `rcoloumb`) are used in the `mdp` files to allow a smaller water box, hence a shorter time to finish the simulations in the tutorials. (An interaction cutoff should be smaller than half of the smallest dimension of the box.)

## 4-site system
The 4-site system is only used in Exercise 4. This is a toy model composed of 4 linearly connected interaction sites. As shown in the presentation, there are two torsional metastable states (cis and trans) for this system due to the diple from the first and last atoms (with charges of +/- 0.2) of the molecule. The files `sys.gro` and `sys.top` were prepared in the folder `Prep` by the script `prep.sh`, with all required inputs for the preparation process provided in the `Prep` folder. `sys.gro` is the configuration extracted from the NPT MD simulation that has a volume closet to the volume averaged overall the simulation trajectory.

