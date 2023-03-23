# Friday (FFS and cFFS Tutorials and Simulation of LiF Dissociation using RETIS and FFS):

# Getting setup and downloading the tutorials:

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
  - “$ git fetch”: This will allow you to download the latest version of contents within the iCoMSE repository. 
 
# Friday Tutorials (Toy models of path sampling methods):

This tutorial will be run through Jupyter notebooks (located in the directory at the left of your screen) via OnDemand. For each Jupyter notebook, please make sure to switch your kernel to “icomse-cpu” kernel.

Note: For participants without Bridges2 access, please refer to this link: [https://github.com/icomse/3rd_workshop_advanced_sampling/blob/main/settingup.md](y) for environment setup instructions

- In this tutorial, you will see Jupyter notebooks for each toy model of each path sampling method: FFS, cFFS. 
  - FFS = Forward flux sampling
  - cFFS = Contour forward flux sampling
 
- Each of these Jupyter notebooks will use langevin_dynamics.py, which uses stochastic dynamics to sample between two states (i.e., A and B) of the potential energy surfaces (PES) we provided.   

- Within each notebook, you will find short descriptions of what each notebook does and the learning objectives for that exercise.


# Simulation of LiF Dissociation using RETIS and FFS

To run this LiF dissociation simulation, with either RETIS or FFS, we will need to access a directory that allows for more storage than the home directory can provide. 

- If you are new to this workshop, you will need to create a directory in the “ocean” directory. To this please run the following command: 

  - “$ mkdir /ocean/projects/see220002p/$(whoami)”

If Running RETIS:

- Then navigate to your Ocean directory by using the following command:
  
  - “$ cd /ocean/projects/see220002p/$(whoami)” 

- Once there, please copy over the files you would need for the demonstration today by running the following command: 
  
  - “$ cp -r ~/[YOUR iCoMSE DIRECTORY]/Friday/* ."

You should see two directories:

- LiF-RETIS: [Readme.txt](https://docs.google.com/document/d/1ZmcKZ1IfSJoPowOvZbf-U9q_t0PGVA3hwY1i7q1FC8M/edit?usp=sharing)

- LiF-FFS: [Readme.txt](https://docs.google.com/document/d/1qZc1vysiBuHnXj2LjjFK-f2bSwkL7enEGSRacGvghh0/edit?usp=sharing)


Note: For participants without Bridges2 access, please refer to this link: [https://github.com/icomse/3rd_workshop_advanced_sampling/blob/main/settingup.md](y) for environment setup instructions.


