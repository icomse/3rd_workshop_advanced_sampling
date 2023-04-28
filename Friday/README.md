## Friday (FFS and cFFS Tutorials and Simulation of LiF Dissociation using RETIS and FFS):

## Slides:
- Slide for FFS/cFFS session can be found [here](https://github.com/icomse/3rd_workshop_advanced_sampling/blob/main/Friday/2023-iCoMSE-PathSampling-P2-FFScFFS.pdf)

## Videos: 
- [Video fpr FFS explanatipn](https://www.youtube.com/watch?v=xfzFO4RkGM8)
- [Video for cFFS explanantion](https://www.youtube.com/watch?v=xfzFO4RkGM8)
- [Video for FFS tutorial](https://www.youtube.com/watch?v=04Dxm57VLMI)

## Getting setup and downloading the tutorials:

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
  - “$ git pull”: This will allow you to download the latest version of contents within the iCoMSE repository.
    - If you are facing an issue due to the fact that you have previous commits, then you'd run do a "$ git revert" first before doing "$ git pull"
 
## Friday Tutorials (Toy models of path sampling methods):

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


Note: For participants without Bridges2 access, please refer to this link: [https://github.com/icomse/3rd_workshop_advanced_sampling/blob/main/settingup.md](y) for environment setup instructions.


