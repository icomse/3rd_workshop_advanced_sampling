# Setting up software for this i-CoMSE workshop:

## Official invitees

* Before we start the workshop, you should test your connection for Bridges2, which we will be using as the primary source of computing resources.  We will be accessing the system in two ways: running a Jupyter notebook, and directly running 

### Logging in via SSH

Use ssh to connect to Bridges-2 using ACCESS credentials and (optionally, if you have it set up) DUO MFA:

* Using your ssh client, use your ACCESS credentials and connect to hostname bridges2.psc.edu.

* ```ssh access-username@bridges2.psc.edu```
* Enter your ACCESS password when prompted.
* If you are registered with ACCESS DUO, you will receive a prompt on your device.  Once you have approved it, you will be logged in.

* If you have a previous Pittsburgh Supercomputing Center (PSC) username and password, you may login using this password, but you will be using your own computational resources. 

### Logging in as a Jupyter Notebook via on-demand.

* To connect to Bridges-2 via OnDemand, point your browser to https://ondemand.bridges2.psc.edu.

* You will be prompted for a username and password.  Enter your your access username prefixed with 'access-' and your ACCESS password.  i.e. if your access password is amy3, then use the username 'access-amy3'. 

*  The OnDemand Dashboard will open.  From this page, you can use the menus across the top of the page to manage files and submit jobs to Bridges-2.
To end your OnDemand session, choose Log Out at the top right of the Dashboard window and close your browser.

* For additional information, go to the "HOME BASE" link sent to you earlier. 

## Accessing the software during the workshop for additional participants watching

### If you have your own account/allocation on PSC

Load our required modules
``` 
module load gromacs/2020.2-cpu
module load plumed
module load openmpi/4.0.5-gcc10.2.0
```
You 

### If you don't have an account on PSC
You can check out this dockerfile. And to copy the dockerfile with singularity, use:

You can install it by running:
'''
singularity pull docker://ghcr.io/cmelab/icomse:latest
'''

If you wish to replicate the python environment locally without gromacs/plumed, you can use the following environment.yml:
name: pathsamplingenv
channels:
  - conda-forge
  - defaults
dependencies:
  - python=3.8
  - dask
  - mdanalysis
  - seaborn
  - jupyter
  - matplotlib
  - numpy

And to get gromacs/plumed via conda:
'''
	conda install -strict-channel-priority -c \
	    plumed/label/masterclass-2022 -c conda-forge plumed gromacs
'''      

## Post-workshop software setup

Instructions for asynchronously running the software after the workshop will be posted here afterwards.
