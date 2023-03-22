import numpy as np
import sys
import math
import os
import re
import glob
import subprocess
import MDAnalysis as mda
import dask
from dask.distributed import Client, LocalCluster
sys.path.insert(1, 'scripts')
import op

# job distribution settings
# number of cores available
n_cores_avail = 32
# number of cores per job
n_cores_per_job = 1
# number of jobs
n_jobs = int(n_cores_avail/n_cores_per_job)

with dask.config.set({"distributed.worker.resources.cores": n_cores_per_job}):
    cluster = LocalCluster(n_workers=n_jobs, threads_per_worker=n_cores_per_job, processes=False)
    
client = Client(cluster)


# DEFINE PATHS HERE ***** EDIT HERE *****
# result directory
MPATH = "/ocean/projects/see220002p/$(whoami)/lif-ffs/"
# program directory
PATH = "/jet/home/$(whoami)/[YOUR iCoMSE DIRECTORY]/Friday/LiF-FFS/"
# basin A trajectory directory
BPATH = "/ocean/projects/see220002p/shared/lif/"
# paths to bash scripts
SHPATH_SIMULATE = PATH + "scripts/run.sh"
SHPATH_SETUP = PATH + "scripts/setup.sh"

# standard filename
STDNAME = "LiF"

# MAKE SOME OTHER GLOBAL DEFS HERE
# Basin Simulation Parameters
FIRST_CROSS = 100

# Interface Simulation Parameters
# Does the OP increase from A-B
GROW = True
# Number of exploring scouts
M = 200
# Number of interface trajectories
N = 1000

# position of A and B basins
A = 1.923518185874568
B = 4.0
# interfaces to run
START = -1
END = 100
# position of starting interface
L = 2.0982887958982905
LAG = 500
N_CROSS = 100

# desired crossing probability
P_DES = 0.1
# minimimum interface displacement
D_MIN = 0.02

# simulation time between OP calculations
TIME = 0.5 # ps
STEPS = int(2000*TIME)  # with timestep = 0.0005



result = [None] * (N)


# run FFS
def main():
    global SHPATH_SETUP
    global PATH
    global MPATH
    global INT_PATH
    global START
    global END
    global L
    done = False
    basin = False
    if START == -1:
        basin = True
        command = SHPATH_SETUP + " " + PATH + " " + MPATH
        subprocess.call(['bash','-c',command])
        flux = basinA("F","Li")
        START = 0
    l_prev = L
    i = 0
    cumuprob = 1
    for i in range(START,END):                                                      # run from START interface to obtain END interface configs
        f = open(MPATH+"stdout.txt", "a")
        f.write(f"Starting interface: {i}\n")
        f.close()
        interface = i
        print(PATH)
        INT_PATH = MPATH + "simulations/int_" + str(interface) + "/"                # directory for current interface
        l = exploring_scouts(interface,l_prev)                                      # place next interface with exploring scouts
        f = open(MPATH+"interfaces.txt", "a")
        f.write("Interface {} is at: {}\n".format(i+1,l))
        f.close()
        cumuprob *= milestone(interface,l)                                                      # run the interface
        f = open(MPATH+"stdout.txt", "a")
        f.write("Cumulative probability from first interface through interface {}: {}\n".format(interface+1,cumuprob))
        f.close()
        l_prev = l
        if l == B:
            done = True
            break
    f = open(MPATH+"stdout.txt", "a")
    if done:
        if basin:
            calc_rate(flux,cumuprob)
        f.write(f"Calculation completed at interface: {i+1}!! Exiting normally...\n")
    else:
        f.write(f"Stopped calculation at interface: {i+1}. Exiting normally...\n")
    f.close()

# analyze basin A trajectory to place basin A and first intcess (PSC) today itself? Do you erface
def basinA(anion,cation):
    global MPATH
    global BPATH
    global STDNAME
    global GROW
    global LAG
    global FIRST_CROSS
    global A
    global L
    global D_MIN
    off = 20000
    gro = BPATH + STDNAME + "-basinA.gro"
    xtc = BPATH + STDNAME + "-basinA.xtc"
    lmda, time = op.op(anion,cation,gro,xtc,offset=off)
    op_list = [lmda]
    avg = np.mean(op_list)
    stdev = np.std(op_list)
    if GROW:
        A = avg + 0.5*stdev
    else:
        A = avg - 0.5*stdev
    f = open(MPATH+"interfaces.txt", "w")
    f.write("Basin A is at: {}\n".format(A))
    f.close()
    f = open(MPATH+"stdout.txt", "w")
    f.write("Basin A is at: {}\n".format(A))
    f.close()
    x_min = min([min(i) for i in op_list])
    x_max = max([max(j) for j in op_list])
    if GROW:
        inter = np.linspace(A,x_max,int((x_max-A)*100))
    else:
        inter = np.linspace(A,x_min,int((A-x_min)*100))
    done = False
    for i in inter:
        n_cross = 0
        for k in range(len(op_list)):
            last_cross = -LAG
            fromBasin = False
            for j in range(len(op_list[k])):
                if GROW:
                    if op_list[k][j] < A and j-last_cross >= LAG:
                        fromBasin = True
                    if fromBasin and op_list[k][j] >= i and j-last_cross >= LAG:
                        if op_list[k][j] < i + D_MIN:
                            n_cross += 1
                            last_cross = j
                        fromBasin = False
                else:
                    if op_list[k][j] > A and j-last_cross >= lag:
                        fromBasin = True
                    if fromBasin and op_list[k][j] <= i and j-last_cross >= LAG:
                        if op_list[k][j] > i - D_MIN:
                            n_cross += 1
                            last_cross = j
                        fromBasin = False                            
        if n_cross >= FIRST_CROSS:
            l0 = i
            done = True
    if not done:
        f = open(MPATH+"stdout.txt", "a")
        f.write("Not enough first crossings, exiting...\n")
        f.close()
        sys.exit("Not enough first crossings")
    else:
        L = l0
    f = open(MPATH+"interfaces.txt", "a")
    f.write("Interface 0 is at: {}\n".format(L))
    f.close()
    f = open(MPATH+"stdout.txt", "a")
    f.write("Interface 0 is at: {}\n".format(L))
    f.close()
    subprocess.call(['mkdir',MPATH+"simulations/int_0/"])
    univ = mda.Universe(gro,xtc)
    cnt = 0
    total_cross = 0
    for k in range(len(op_list)):
        last_cross = -LAG
        fromBasin = False
        for j in range(len(op_list[k])):
            if GROW:
                if op_list[k][j] < A and j-last_cross >= LAG:
                    fromBasin = True
                elif op_list[k][j] >= l0 and fromBasin:
                    if j-last_cross >= LAG:
                        if op_list[k][j] < l0 + D_MIN:
                            univ.trajectory[j+off]
                            with mda.Writer(MPATH+"simulations/int_0/"+STDNAME+f"_0_{cnt}_0.gro", univ.atoms.n_atoms) as W:
                                W.write(univ)
                            cnt += 1
                            last_cross = j
                    fromBasin = False
                    total_cross += 1
            else:
                if op_list[k][j] > A and j-last_cross >= LAG:
                    fromBasin = True
                elif op_list[k][j] <= l0 and fromBasin:
                    if j-last_cross >= LAG:
                        if op_list[k][j] < l0 + D_MIN:
                            univ.trajectory[j+off]
                            with mda.Writer(MPATH+"simulations/int_0/"+STDNAME+f"_0_{cnt}_0.gro", univ.atoms.n_atoms) as W:
                                W.write(univ)
                            cnt += 1
                            last_cross = j
                    fromBasin = False
                    total_cross += 1
    total_time = time[-1]-time[off]
    flux = total_cross/total_time
    f = open(MPATH+"stdout.txt", "a")
    f.write("Total flux through the first interface: {} ps^-1\n".format(flux))
    f.close()
    return flux
    
    
        
# perform exploring scouts to place next interface
def exploring_scouts(interface,l_prev):
    global INT_PATH
    global STDNAME
    global result
    global M
    global D_MIN
    global GROW
    global B
    configs = []
    for f in os.listdir(INT_PATH):
        name, ext = os.path.splitext(f)
        if ext == '.gro':
            configs += [f]                                                          # list of every config file
    length = len(configs)
    #f = open(MPATH+"stdout.txt", "a")
    #f.write(f"List of configurations: {configs}\n")
    #f.close()
    s = [None] * length                                                           
    sims = [int(M/length)] * length                                                 # distribute the simulations
    for i in range(length):
        if i < M%length:
            sims[i] += 1
    for i in range(length):
        s[i] = [int(float(n)) for n in re.findall(r'-?\d+\.?\d*', configs[i])]      # capture all the numbers in config filename
    for i in range(M):
        command = MPATH + "simulations/int_" + str(interface) + "/ex-sim_" + str(i) # make a directory for each sim
        subprocess.call(['mkdir',command])                                         
    cnt = 0
    for i in range(length):
        for j in range(sims[i]):
            module = STDNAME + "_" + str(interface) + "_" + str(s[i][-2]) + "_" \
            + str(s[i][-1]) + "_" + str(i) + "_" + str(j)
            command = INT_PATH + "ex-sim_" + str(cnt) + "/" + module + ".gro"
            subprocess.call(['cp',INT_PATH + configs[i],command])                   # copy configs ito sim directories and run
            result[cnt] = client.submit(run,interface,B,cnt,module,"y",pure=False,resources={'cores': 1})
            cnt += 1
    for i in range(M):
        if result[i].result():                                                      # ensure all runs finish before continuing
            val = True
    l = explore("F","Li",interface)                                                # calculate next interface location
    f = open(MPATH+"stdout.txt", "a")
    f.write(f"Exploring scouts interface at: {l}\n")
    f.close()
    if l_prev != 0:
        if GROW:
            if l > B:
                l = B
                f = open(MPATH+"stdout.txt", "a")
                f.write(f"Next interface is past B, moving interface to: {l}\n")
                f.close()
            elif l - l_prev < D_MIN:
                l = l_prev + D_MIN
                f = open(MPATH+"stdout.txt", "a")
                f.write(f"Next interface is not far enough away, moving interface to: {l}\n")
                f.close()
        else:
            if l < B:
                l = B
                f = open(MPATH+"stdout.txt", "a")
                f.write(f"Next interface is past B, moving interface to: {l}\n")
                f.close()
            if l_prev - l < D_MIN:                                                      # enforce minimum interface displacement
                l = l_prev - D_MIN
                f = open(MPATH+"stdout.txt", "a")
                f.write(f"Next interface is not far enough away, moving interface to: {l}\n")
                f.close()
    f = open(MPATH+"stdout.txt", "a")
    f.write(f"The next interface is at: {l}\n")
    f.close()
    return l
    
# calculate next interface location based on desired crossing probability
def explore(anion,cation,interface):
    global GROW
    global P_DES
    global M
    global STDNAME
    global INT_PATH
    global A
    l_min = np.zeros([2,M])
    l_max = np.zeros([2,M])
    for i in range(M):                                                             # calculate max OP value for each exploring scout
        lmda, time = op.op(anion,cation,glob.glob(INT_PATH+"ex-sim_"+str(i)+"/*.gro")[0],
                                         glob.glob(INT_PATH+"ex-sim_"+str(i)+"/*.xtc")[0])
        for j in lmda:
            if GROW:
                if j < A:
                    break
            else:
                if j > A:
                    break
            if l_min[0,i] == 0:
                l_min[0,i] = j
            elif j < l_min[0,i]:
                l_min[0,i] = j
            if l_max[0,i] == 0:
                l_max[0,i] = j
            elif j > l_max[0,i]:
                l_max[0,i] = j
        l_min[1,i] = i
        l_max[1,i] = i
    if GROW:
        l_sort = l_max[:, l_max[0,:].argsort()]
    else:
        l_sort = l_min[:, l_min[0,:].argsort()]                                        # sort OP values to select for desired crossing probability
    s = math.floor(M*P_DES)
    #f = open(MPATH+"stdout.txt", "a")
    #f.write(f"Sorted exploratory OP values: \n{l_sort[0]}\n")
    #f.close()
    if GROW:
        return l_sort[0,-s]
    else:
        return l_sort[0,s-1]
    
# simulate trajectories to collect configurations at the next interface
def milestone(interface,l):
    global INT_PATH
    global MPATH
    global STDNAME
    global result
    global N
    global M
    configs = []
    for f in os.listdir(INT_PATH):
        name, ext = os.path.splitext(f)
        if ext == '.gro':
            configs += [f]                                                         # list of every config file
    length = len(configs)
    s = [None] * length
    for i in range(length):
        s[i] = [int(float(n)) for n in re.findall(r'-?\d+\.?\d*', configs[i])]     # capture all the numbers in config filename
    sims = [int(N/(length))] * length
    for i in range(length):                                                        # distribute the simulations           
        if i < N%length:
            sims[i] += 1
    for i in range(N):
        command = MPATH + "simulations/int_" + str(interface) + "/sim_" + str(i)   # make a directory for each sim
        subprocess.call(['mkdir',command])
    cnt = 0
    conf_cnt = 0
    command = MPATH + "simulations/int_" + str(interface+1)
    subprocess.call(['mkdir',command])
    for i in range(length):
        for j in range(sims[i]):
            #naming convention: module_(int)_(prev_conf)_(prev_traj)_(conf)_(traj)
            module = STDNAME + "_" + str(interface) + "_" \
            + str(s[i][-2]) + "_" + str(s[i][-1]) + "_" + str(conf_cnt) + "_" + str(cnt)
            command = INT_PATH + "sim_" + str(cnt) + "/" + module + ".gro"
            subprocess.call(['cp',INT_PATH + configs[i],command])                  # copy configs ito sim directories and run
            result[cnt] = client.submit(run,interface,l,cnt,module,"n",pure=False,resources={'cores': 1})
            cnt += 1
        conf_cnt += 1
    for i in range(N):
        if result[i].result():                                                     # ensure all runs finish before continuing
            val = True
    new_configs = []
    for f in os.listdir(MPATH + "simulations/int_" + str(interface+1) + "/"):
        name, ext = os.path.splitext(f)
        if ext == '.gro':
            new_configs += [f]
    prob = len(new_configs)/N
    f = open(MPATH+"stdout.txt", "a")
    f.write("Crossing probability through interface {}: {}\n".format(interface+1,prob))
    f.close()
    return prob
    
    
# run simulation
def run(interface,l,traj_num,filename,explore):
    global GROW
    global MPATH
    global PATH
    global A
    global B
    global TIME
    global STEPS
    global SHPATH_SIMULATE
    #arguments [PATH] [MPATH] [grofilename] [interface] [traj_num] [basinA] [lambda] [basinB] [TIME] [STEPS] [exploringscouts (y/n)]
    command = SHPATH_SIMULATE + " " + PATH + " " + MPATH + " " + filename + " " \
              + str(interface) + " " + str(traj_num) + " " + str(A) + " " \
              + str(l) + " " + str(B) + " " + str(TIME) + " " + str(STEPS) + " " + explore + " " + str(GROW)
    subprocess.call(['bash','-c',command])
    #subprocess.call(['bash','-c',command])
    return True
    
def calc_rate(flux,cumuprob):
    rate = flux*cumuprob
    f = open(MPATH+"stdout.txt", "a")
    f.write("The rate of LiF dissociation: {} ps^-1\n".format(rate))
    f.close()
    
    
if __name__ == "__main__":
    main()


