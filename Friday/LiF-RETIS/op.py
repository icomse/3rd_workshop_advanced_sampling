# Code Written by PH Minh
# Date Last Modified: 03-15-2023
# OP code for LiF Ion Dissociation
# SAMPEL Group
# UMN

import sys
import os.path
import numpy as np
import pandas as pd
import MDAnalysis as mda
import argparse

# Define parser 
def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
            description="OP Code for Ion Dissociation")
    parser.add_argument("-g", "--grofile", help="specify the .gro file", default="input.gro", metavar='')
    parser.add_argument("-x", "--xtcfile", help="specify the trajectory file", default="input.xtc", metavar='')
    parser.add_argument("-t", "--trrfile", help="specify the trajectory file", default="input.trr", metavar='')
    parser.add_argument("-o", "--output", help="specify the name of the output file", default="out.txt", metavar='')

    args = parser.parse_args()
    return args

# Call the function
args = parseargs()

def main():
    args = parseargs()
    gro = args.grofile
    xtc = args.xtcfile
    trr = args.trrfile
    out = args.output

    if os.path.isfile(gro) and os.path.isfile(trr):
    #if os.path.isfile(gro) and os.path.isfile(xtc):
        
        # Load in the gro and trajectory file
        system = mda.Universe(gro,trr)
        #system = mda.Universe(gro,xtc)
        
        # Define your ion
        F = system.select_atoms("resname F")
        Li = system.select_atoms("resname Li")
        
        # Get time and distance 
        dis = []
        time = []
        cutoff=-1
        skip=1
        #offset=0

        #for ts in system.trajectory[offset:cutoff:skip]:
        for ts in system.trajectory[:cutoff:skip]:
           
            time.append(system.trajectory.time)        # Time in ps
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

        np.savetxt(str(out), np.transpose((time,dis)), fmt='%10.5f')
    
    else:
        print("Either the gro file or the traj. file is incorrect")
        print("Current gro filename {}, current trj filename {}".format(gro,trr))
        exit()
        

# Boilerplate code to call main() fn. 
if __name__ == '__main__':
    main()
