Enhanced sampling molecular dynamics on HPC systems
===================================================

## 1. Prerequisites
To run the Jupyter notebook `MD_HPC.ipynb` on Bridges-2, follow the steps below:
- Login to Bridges-2 and launch and interactive session using the command: `interact -N 1 --ntasks-per-node=4`.
- Find the hostname of the node you are running on using the `hostname` command in your terminal. For example:
  ```
  [wehs7661@r001 ~]$ hostname
  r001.ib.bridges2.psc.edu
  [wehs7661@r001 ~]$
  ```
- Start a Jupyter notebook with the command `jupyter notebook --no-browser --ip=0.0.0.0 &`. The output will provide a port number and token that you can use to access the notebook from your local machine. In the following example, the port number is 8889 and the token is `ec46ed9e04e9ba7b05d2be5d6526cbb3f8dfc23501027a4f`.  Note that this command runs on the background, but you can run it in the foreground by omitting the `&` at the end.
  ```
  (base) [wehs7661@r001 wehs7661]$ jupyter notebook --no-browser --ip=0.0.0.0 &
  [I 19:24:45.113 NotebookApp] Serving notebooks from local directory: /ocean/projects/cts160011p/wehs7661
  [I 19:24:45.113 NotebookApp] Jupyter Notebook 6.4.9 is running at:
  [I 19:24:45.113 NotebookApp] http://r001.ib.bridges2.psc.edu:8888/?token=ec46ed9e04e9ba7b05d2be5d6526cbb3f8dfc23501027a4f
  [I 19:24:45.113 NotebookApp]  or http://127.0.0.1:8889/?token=ec46ed9e04e9ba7b05d2be5d6526cbb3f8dfc23501027a4f
  [I 19:24:45.113 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
  [C 19:24:45.120 NotebookApp]

      To access the notebook, open this file in a browser:
          file:///jet/home/wehs7661/.local/share/jupyter/runtime/nbserver-39738-open.html
      Or copy and paste one of these URLs:
          http://r001.ib.bridges2.psc.edu:8888/?token=ec46ed9e04e9ba7b05d2be5d6526cbb3f8dfc23501027a4f
  ```
- On a new terminal and connect to Bridges-2 using the following command. (Remember to place `username` with your own username.)
  ```
  ssh -L 8888:r001.ib.bridges2.psc.edu:8889   bridges2.psc.edu -l username
  ```
   The 8888 is the local port number, while r001.ib.bridges2.psc.edu:8889 is the long name of the compute node and the port number where Jupyter is running. Note that even though the port on the compute node is 8889 in this case, we still use local port 8888  as usual.
- Finally, on your local machin, open a browser window and point it to http://localhost:8888. Enter the token to make the connection with the Bridges-2 compute node. When you are done, close your interactive session on Bridges-2. 

## 2. Outline
Here is the outline of this session:
- Introduction to High-Performing Computing (HPC) systems
- [Exercise 1](): Performing decoupled MD simulations in parallel
- The theory of umbrella sampling
- [Exercise 2](): Performing umbrella sampling in parallel
- [Exercise 3](): Analyzing umbrella sampling with MBAR
- The theory of multiple walkers metadynamics
- [Exercise 4](): Performing and analyzing multiple walkers metadynamics
- The theory of replica exchange molecular dynamics
- [Exercise 5](): Performing and analyzing Hamiltonian replica exchange
- The theory of replica exchange umbrella sampling
- [Exercise 6](): Performing and analyzing replica exchange umbrella sampling
