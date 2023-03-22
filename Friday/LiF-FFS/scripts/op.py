import numpy as np
import MDAnalysis as mda
import sys
import math

# run through run.sh to calculate the OP to continue simulation or write out configuration 
def main():
    # [operation]
    if sys.argv[1] == "run":
        # [anion] [cation] [config_file] [traj_file]
        lmda, time = op(sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])       # calculate OP
        for i in lmda:
            # [basinA] [l] [grow]
            if sys.argv[8] == "True":
                if i < float(sys.argv[6]):
                    print(1)
                    return
                elif i >= float(sys.argv[7]):
                    print(2)
                    return
            elif sys.argv[8] == "False":
                if i > float(sys.argv[6]):
                    print(1)
                    return
                elif i <= float(sys.argv[7]):
                    print(2)
                    return   
            else:
                sys.exit("Error in GROW: "+sys.argv[8])
        print(0)
    elif sys.argv[1] == "end":
        lmda, time = op(sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])
        for i in range(len(time)):
            if sys.argv[8] == "True":
                if lmda[i] >= float(sys.argv[7]) or lmda[i] < float(sys.argv[6]):
                    print("%.4f" % time[i+1])
                    return
            else:
                if lmda[i] <= float(sys.arv[7]) or lmda[i] > float(sys.argv[6]):
                    print("%.4f" % time[i+1])
                    return
    elif sys.argv[1] == "write":                                                 # write next interface config
        # [anion] [cation] [config_file] [traj_file] [l] [next_int_path] [grow]
            op_w(sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],float(sys.argv[6]),sys.argv[7],sys.argv[8]) 

# calculate OP [anion] [cation] [config_file] [traj_file]
def op(anion,cation,gro,xtc,offset=0,cutoff=-1,skip=1):
    system = mda.Universe(gro,xtc)
    F = system.select_atoms("resname F")
    Li = system.select_atoms("resname Li")
    dis = []
    time = []
    for ts in system.trajectory[offset:cutoff:skip]:
        time.append(system.trajectory.time)
        F_pos = F.positions[0]
        Li_pos = Li.positions[0]
        box = system.dimensions[0:3]
        x = abs(Li_pos[0]-F_pos[0])
        y = abs(Li_pos[1]-F_pos[1])
        z = abs(Li_pos[2]-F_pos[2])
        dx = min(x,abs(box[0]-x))
        dy = min(y,abs(box[1]-y))
        dz = min(z,abs(box[2]-z))
        dis.append((dx**2+dy**2+dz**2)**0.5)
    return dis, time

# write next interface config [anion] [cation] [config_file] [traj_file] [next_int] [next_int_path] [grow]
def op_w(anion,cation,gro,xtc,l,write,grow):
    system = mda.Universe(gro,xtc)
    lmda, time = op(anion,cation,gro,xtc)
    print("writing to "+write)
    for i in range(len(time)):
        if grow == "True":
            if lmda[i] >= l:
                system.trajectory[i]
                break
        else:
            if lmda[i] <= l:
                system.trajectory[i]
                break
    with mda.Writer(write, system.atoms.n_atoms) as W:
        W.write(system)

if __name__ == "__main__":
    main()
