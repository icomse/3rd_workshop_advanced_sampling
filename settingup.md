# Setting up software for this i-CoMSE workshop:

## Official invitees

* Before we start the workshop, you should test your connection for Bridges2, which we will be using as the primary source of computing resources.  We will be accessing the system in two ways: running a Jupyter notebook, and directly running 

### Logging in via SSH

Set your Bridges-2 password(if you are not already aware of what your password is). To do this, go to https://apr.psc.edu and click Start. Enter in your Bridges-2 username and the email address that you created your ACCESS account with, then click Submit. (If you do not know your Bridges-2 Username, please visit this page in ACCESS: https://allocations.access-ci.org/profile). You will receive an email containing a security code. Enter that into the box listed "Security Code" and then again click Submit. Choose your PSC password, noting the password rules on the page, and click Submit again.

* Using your ssh client, use this Bridges-2 username and your new password.

* If you have a previous Pittsburgh Supercomputing Center (PSC) username and password, you may login using this password, but you will end up using your own computational resources. 

### Logging in as a Jupyter Notebook via on-demand.

* To connect to Bridges-2 via OnDemand, point your browser to https://ondemand.bridges2.psc.edu.

* You will be prompted for a username and password.  These are your Bridges-2 username and password from above.

*  The OnDemand Dashboard will open.  From this page, you can use the menus across the top of the page to manage files and submit jobs to Bridges-2.
See the "Home Base" document for specific settings. To end your OnDemand session, choose Log Out at the top right of the Dashboard window and close your browser.

* For additional information, go to the "Home Base" document link sent to you earlier. 

## Accessing the software during the workshop for additional participants watching

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

## Post-workshop software setup

Instructions for asynchronously running the software after the workshop will be posted in this section afterwards.
