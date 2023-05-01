# Setting up software for the material from this i-CoMSE workshop:

### If you have your own account/allocation on PSC

Load our required modules
``` 
module load gromacs/2020.2-cpu
module load plumed
module load openmpi/4.0.5-gcc10.2.0
```
### If you don't have an account on PSC
You can get a copy of the the dockerfile we're using with singularity using:

```
singularity pull docker://ghcr.io/cmelab/icomse:latest
```
In principle this container can be used locally on your desktop/laptop, or on your singularity-enabled HPC cluster, but we unfortunately can't provide support getting that set up.

Most of the workshop functionality is provided through python libraries. If you wish to replicate the python environment from the container locally without gromacs use this [environment.yml](environment.yml):

```
conda env create -f environment.yml
conda activate icomseW23
```

And to get gromacs/plumed via conda:
```
	conda install --strict-channel-priority -c \
	    plumed/label/masterclass-2022 -c conda-forge plumed gromacs
```
