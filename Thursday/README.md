# Path Sampling Tutorials

## Instructor/Authors

**Dr. Sapna Sarupria**
- Department of Chemistry, University of Minnesota
- https://sarupriagroup.github.io/
- sarupria@umn.edu

**Tutorials**
- PH (Porhouy Minh)
- Naomi Trampe

## Slides 
- Slides for for TPS/TIS/RETIS session can be found [here](https://github.com/icomse/3rd_workshop_advanced_sampling/blob/main/Thursday/2023-iCoMSE-PathSampling-P1-TPSTIS.pdf)

## Getting setup and downloading the tutorials ((NOTE: much of the material about accessing computing resources was for the workshop, and will not be available for people viewing the information after the workshop)

- To get this started, please login to Bridges2 OnDemand via: https://ondemand.bridges2.psc.edu/pun/sys/dashboard/

- Then navigate to Interactive Apps → Jupyter Lab 
  - Number of hours = 3
  - Number of nodes = 1
  - Account = see220002p
  - Partition = RM-shared**
  - Extra Slurm Args = -n 1

- After that, click on “Connect to Jupyter”. From here you should see a startup page which should have an option for you to open a “Terminal”.

- Then go into your iCoMSE directory by using the command: 
  - “$ cd [YOUR iCoMSE DIRECTORY]”

- Then using the command
  - “$ git fetch”: This will allow you to download the latest version of contents within the iCoMSE repository. (note: you will need to first clone this repository before fetching)
 
 In addition to the Monday - Wednesday directories, you should now also see Thursday and Friday directories 
 
## Thursday Tutorials (Toy models of path sampling methods):

This tutorial will be run through Jupyter notebooks (located in the directory at the left of your screen) via OnDemand. For each Jupyter notebook, please make sure to switch your kernel to “icomse-cpu” kernel.

Note: For participants without Bridges2 access, please refer to this link: [https://github.com/icomse/3rd_workshop_advanced_sampling/blob/main/settingup.md](y) for environment setup instructions

- In this tutorial, you will see Jupyter notebooks for each toy model of each path sampling method: TPS, TIS, RETIS. 
  - TPS = Transition path sampling
  - TIS = Transition interface sampling
  - RETIS = Replica exchange transition interface sampling 
 
- Each of these Jupyter notebooks will use langevin_dynamics.py, which uses stochastic dynamics to sample between two states (i.e., A and B) of the potential energy surfaces (PES) we provided.   

- Within each notebook, you will find short descriptions of what each notebook does and the learning objectives for that exercise.



