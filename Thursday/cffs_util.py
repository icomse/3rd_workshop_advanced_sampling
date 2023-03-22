"""

This module is designed to provide functions which
are identical for many cffs scripts. 

"""

import sys
import math
import numpy as np
import langevin_dynamics as ld
import matplotlib.pyplot as plt
from matplotlib import colors, ticker, cm

# Code written for python 3.4
# Ensure correct version
min_version = (3,4)
if sys.version_info < min_version:
    sys.exit("Must be run with at least python 3.4")

# Function to address 'weird' IEEE 754 behavior
def normal_round(n):
    if n - math.floor(n) < 0.5:
        return math.floor(n)
    return math.ceil(n)

# Class that defines the grid we use to define
# interface sets in cffs. Defined to help pass
# around fewer vars.
class Grid:
    def __init__(self,gridsize,gridmin,gridmax):
        self.size_cv1 = gridsize[0]
        self.size_cv2 = gridsize[1]
        self.min_cv1 = gridmin[0]
        self.min_cv2 = gridmin[1]
        self.max_cv1 = gridmax[0]
        self.max_cv2 = gridmax[1]
        self.nbins_cv1 = int((gridmax[0]-gridmin[0])/gridsize[0])
        self.nbins_cv2 = int((gridmax[1]-gridmin[1])/gridsize[1])

# Function to plot set with contours
def plot_set(x_range,y_range,contours,sites,grid,pointsize,pes_type):
    N = 100
    x_vec = np.linspace(x_range[0], x_range[1], N)
    y_vec = np.linspace(y_range[0], y_range[1], N)
    energy = np.zeros((N, N))
    for i in range(len(x_vec)):
        for j in range(len(y_vec)):
            energy[j][i] = ld.potential(pes_type,x_vec[i],y_vec[j])
    locs = []
    for site in sites: 
        cords = [(site[0]+0.5)*grid.size_cv1 + grid.min_cv1,(site[1]+0.5)*grid.size_cv2 + grid.min_cv2] 
        locs.append(cords)
    plt.contour(x_vec,y_vec,energy,[-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3])
    locs = np.asarray(locs)
    plt.plot(locs[:,0],locs[:,1],'bs',markersize=pointsize)
    plt.xlim((x_range[0],x_range[1]))
    plt.ylim((y_range[0],y_range[1]))
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()    
    
def plot_edges(x_range,y_range,contours,sites,avg_edgecount,grid,pointsize,pes_type):
    N = 100
    x_vec = np.linspace(x_range[0], x_range[1], N)
    y_vec = np.linspace(y_range[0], y_range[1], N)
    energy = np.zeros((N, N))
    for i in range(len(x_vec)):
        for j in range(len(y_vec)):
            energy[j][i] = ld.potential(pes_type,x_vec[i],y_vec[j])
    locs = []
    for site in sites: 
        cords = [(site[0]+0.5)*grid.size_cv1 + grid.min_cv1,(site[1]+0.5)*grid.size_cv2 + grid.min_cv2, avg_edgecount[site]] 
        locs.append(cords)
    plt.contour(x_vec,y_vec,energy,[-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3])
    locs = np.asarray(locs)
    plt.scatter(locs[:,0],locs[:,1],c=locs[:,2],marker='s',s=pointsize)
    plt.xlim((x_range[0],x_range[1]))
    plt.ylim((y_range[0],y_range[1]))
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()

def plot_set_cross(x_range,y_range,contours,basin_A_sites,set_sites,first_crosses,grid,pointsize,pes_type):
    N = 100
    x_vec = np.linspace(x_range[0], x_range[1], N)
    y_vec = np.linspace(y_range[0], y_range[1], N)
    energy = np.zeros((N, N))
    for i in range(len(x_vec)):
        for j in range(len(y_vec)):
            energy[j][i] = ld.potential(pes_type,x_vec[i],y_vec[j])

    plt.contour(x_vec,y_vec,energy,[-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3])
    locs = []
    for site in set_sites: 
        cords = [(site[0]+0.5)*grid.size_cv1 + grid.min_cv1,(site[1]+0.5)*grid.size_cv2 + grid.min_cv2] 
        locs.append(cords)      
    locs = np.asarray(locs)
    plt.plot(locs[:,0],locs[:,1],'gs',markersize=pointsize)
    locs = []
    for site in basin_A_sites: 
        cords = [(site[0]+0.5)*grid.size_cv1 + grid.min_cv1,(site[1]+0.5)*grid.size_cv2 + grid.min_cv2] 
        locs.append(cords)      
    locs = np.asarray(locs)
    plt.plot(locs[:,0],locs[:,1],'bs',markersize=pointsize)
    first_crosses = np.asarray(first_crosses)
    plt.plot(first_crosses[:,0],first_crosses[:,1],'ro',markersize=pointsize)
    plt.xlim((x_range[0],x_range[1]))
    plt.ylim((y_range[0],y_range[1]))
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()

# Extend or 'grow' the set in all directions
def extend_set_everywhere(prev_set,basin_B_set,grid):
    newset = set()
    for site in prev_set:
        bin_cv1 = site[0]
        bin_cv2 = site[1]
        # First add site 
        newset = add_site(newset,site,basin_B_set,grid)
        # Then add neighbors
        newset = add_site(newset,(bin_cv1+1,bin_cv2),basin_B_set,grid)
        newset = add_site(newset,(bin_cv1-1,bin_cv2),basin_B_set,grid)
        newset = add_site(newset,(bin_cv1,bin_cv2+1),basin_B_set,grid)
        newset = add_site(newset,(bin_cv1,bin_cv2-1),basin_B_set,grid)

    return newset

# Extend or 'grow' the set based upon the edge crossing density
def extend_set_selective(newset,edges,avg_edgecount,density_threshold,basin_B_set,grid):
    for site in edges:
        bin_cv1 = site[0]
        bin_cv2 = site[1]
        # Grow edges of set with more crosses more quickly
        if avg_edgecount[site] > 3.*density_threshold:
            # Add first layer
            if bin_cv1-1 >= 0 and (bin_cv1-1,bin_cv2) not in basin_B_set:
                newset.add((bin_cv1-1,bin_cv2))
            if bin_cv1+1 < grid.nbins_cv1 and (bin_cv1+1,bin_cv2) not in basin_B_set:
                newset.add((bin_cv1+1,bin_cv2))
            if bin_cv2-1 >= 0 and (bin_cv1,bin_cv2-1) not in basin_B_set:
                newset.add((bin_cv1,bin_cv2-1))
            if bin_cv2+1 < grid.nbins_cv2 and (bin_cv1,bin_cv2+1) not in basin_B_set:
                newset.add((bin_cv1,bin_cv2+1))
            # Add corners
            if bin_cv1-1 >= 0 and bin_cv2-1 >= 0 and (bin_cv1-1,bin_cv2-1) not in basin_B_set:
                newset.add((bin_cv1-1,bin_cv2-1))
            if bin_cv1+1 < grid.nbins_cv1 and bin_cv2-1 >= 0 and (bin_cv1+1,bin_cv2-1) not in basin_B_set:
                newset.add((bin_cv1+1,bin_cv2-1))
            if bin_cv1-1 >= 0 and bin_cv2+1 < grid.nbins_cv2 and (bin_cv1-1,bin_cv2+1) not in basin_B_set:
                newset.add((bin_cv1-1,bin_cv2+1))
            if bin_cv1+1 < grid.nbins_cv1 and bin_cv2+1 < grid.nbins_cv2 and (bin_cv1+1,bin_cv2+1) not in basin_B_set:
                newset.add((bin_cv1+1,bin_cv2+1))
            # Add second layer
            if bin_cv1-2 >= 0 and (bin_cv1-2,bin_cv2) not in basin_B_set:
                newset.add((bin_cv1-2,bin_cv2))
            if bin_cv1+2 < grid.nbins_cv1 and (bin_cv1+2,bin_cv2) not in basin_B_set:
                newset.add((bin_cv1+2,bin_cv2))
            if bin_cv2-2 >= 0 and (bin_cv1,bin_cv2-2) not in basin_B_set:
                newset.add((bin_cv1,bin_cv2-2))
            if bin_cv2+2 < grid.nbins_cv2 and (bin_cv1,bin_cv2+2) not in basin_B_set:
                newset.add((bin_cv1,bin_cv2+2))
        # Grow edges of set with fewer crosses less quickly
        elif avg_edgecount[site] > density_threshold:
            if bin_cv1-1 >= 0 and (bin_cv1-1,bin_cv2) not in basin_B_set:
                newset.add((bin_cv1-1,bin_cv2))
            if bin_cv1+1 < grid.nbins_cv1 and (bin_cv1+1,bin_cv2) not in basin_B_set:
                newset.add((bin_cv1+1,bin_cv2))
            if bin_cv2-1 >= 0 and (bin_cv1,bin_cv2-1) not in basin_B_set:
                newset.add((bin_cv1,bin_cv2-1))
            if bin_cv2+1 < grid.nbins_cv2 and (bin_cv1,bin_cv2+1) not in basin_B_set:
                newset.add((bin_cv1,bin_cv2+1))
    
    return newset

# Add a new site if possible
# For non-periodic cvs 
def add_site(currentset,newsite,basin_B_set,grid):
    bin_cv1 = newsite[0]
    bin_cv2 = newsite[1]
    if bin_cv1 < 0 or bin_cv2 < 0:
        return currentset
    elif bin_cv1 >= grid.nbins_cv1 or bin_cv2 >= grid.nbins_cv2:
        return currentset
    else:
        if (bin_cv1,bin_cv2) not in basin_B_set:
            currentset.add((bin_cv1,bin_cv2))

    return currentset

# Find the edges of a set
# Again for non-periodic cvs
def loc_edges(sset):
    # Set object to store edges
    edges = set()
    for site in sset:
        bin_cv1 = site[0]
        bin_cv2 = site[1]
        countn = 0
        if (bin_cv1-1,bin_cv2) in sset:
            countn+=1
        if (bin_cv1+1,bin_cv2) in sset:
            countn+=1
        if (bin_cv1,bin_cv2-1) in sset:
            countn+=1
        if (bin_cv1,bin_cv2+1) in sset:
            countn+=1
        # If countn < 4 we have found an edge
        if countn < 4:
            edges.add(site)

    return edges

def idsetcls(sset,grid):
    # Find the largest touching set
    # of things that DONT meet the threshold
    cls = []
    clscount = 0
    incls = set()
    for site in sset:
        if site not in incls:
            # Start a new cluster and add
            # the current site to the cluster
            cls.append(set())
            incls.add(site)
            cls[clscount].add(site)
            clscheck = []
            clscheck.append(site)
            while clscheck:
                y = clscheck[0]
                del clscheck[0]
                if y[0]-1 >= 0:
                    trial = (y[0]-1,y[1])
                    if trial in sset and trial not in incls:
                        cls[clscount].add(trial)
                        incls.add(trial)
                        clscheck.append(trial)
                if y[0]+1 < grid.nbins_cv1:
                    trial = (y[0]+1,y[1])
                    if trial in sset and trial not in incls:
                        cls[clscount].add(trial)
                        incls.add(trial)
                        clscheck.append(trial)
                if y[1]-1 >= 0:
                    trial = (y[0],y[1]-1)
                    if trial in sset and trial not in incls:
                        cls[clscount].add(trial)
                        incls.add(trial)
                        clscheck.append(trial)
                if y[1]+1 < grid.nbins_cv2:
                    trial = (y[0],y[1]+1)
                    if trial in sset and trial not in incls:
                        cls[clscount].add(trial)
                        incls.add(trial)
                        clscheck.append(trial)
            # Now increment counter for clusters
            clscount+=1
    return cls

def idnosetcls(newset,grid):
    # Find the largest touching set
    # of sites that DONT meet the threshold
    # This function helps us remove 'holes'
    # in the interface set.
    cls = []
    clscount = 0
    incls = set()
    for bin_cv1 in range(grid.nbins_cv1):
        for bin_cv2 in range(grid.nbins_cv1):
            if (bin_cv1,bin_cv2) not in newset:
                if (bin_cv1,bin_cv2) not in incls:
                    # Start a new cluster and add
                    # the current site to the cluster
                    cls.append(set())
                    incls.add((bin_cv1,bin_cv2))
                    cls[clscount].add((bin_cv1,bin_cv2))
                    clscheck = []
                    clscheck.append((bin_cv1,bin_cv2))
                    while clscheck:
                        y = clscheck[0]
                        del clscheck[0]
                        if y[0]-1 >= 0:
                            trial = (y[0]-1,y[1])
                            if trial not in newset and trial not in incls:
                                cls[clscount].add(trial)
                                incls.add(trial)
                                clscheck.append(trial)
                        if y[0]+1 < grid.nbins_cv1:
                            trial = (y[0]+1,y[1])
                            if trial not in newset and trial not in incls:
                                cls[clscount].add(trial)
                                incls.add(trial)
                                clscheck.append(trial)
                        if y[1]-1 >= 0:
                            trial = (y[0],y[1]-1)
                            if trial not in newset and trial not in incls:
                                cls[clscount].add(trial)
                                incls.add(trial)
                                clscheck.append(trial)
                        if y[1]+1 < grid.nbins_cv1:
                            trial = (y[0],y[1]+1)
                            if trial not in newset and trial not in incls:
                                cls[clscount].add(trial)
                                incls.add(trial)
                                clscheck.append(trial)
                    # Now increment counter for clusters
                    clscount+=1
    return cls

def calc_exit_basin(basintraj,grid,basin_A_set,lambda_0_set,edges):
    # Initialize vars
    from_basin = False
    first_crosses = []
    edgecount = {}
    basin_count = 0
    cross_count = 0
    loc = (-10,-10)
    for timestep in basintraj:
        prevloc = loc
        bin_cv1 = int((timestep[0]-grid.min_cv1)/grid.size_cv1)
        bin_cv2 = int((timestep[1]-grid.min_cv2)/grid.size_cv2)
        loc = (bin_cv1,bin_cv2)
        if loc in basin_A_set:
            from_basin = True
            basin_count += 1
        else:
            if loc not in lambda_0_set and from_basin == True:
                # Found a crossing
                from_basin = False
                cross_count += 1
                first_crosses.append(timestep)
                # And track where it exited from
                if prevloc in edges:
                    if prevloc in edgecount:
                        edgecount[prevloc] += 1
                    else:
                        edgecount[prevloc] = 1

    # Now we average the edge counts to smooth data a bit
    # Average edge counts with edges +-3 sites 
    # in cv1/cv2 directions
    avg_edgecount = {}
    for site in edges:
        bin_cv1 = site[0]
        bin_cv2 = site[1]
        avg_edgecount[site] = 0
        count = 0
        for d_cv1 in range(-3,4):
            for d_cv2 in range(-3,4):
                loc = (bin_cv1+d_cv1,bin_cv2+d_cv2)
                if loc in edges:
                    if loc in edgecount:
                        avg_edgecount[site] += edgecount[loc]
                    count +=1 
        avg_edgecount[site] /= count

    return first_crosses, cross_count, basin_count, avg_edgecount

def calc_exit(inttrajs,grid,basin_A_set,basin_B_set,active_set,edges):
    first_crosses = []
    edgecount = {}
    basin_count = 0
    cross_count = 0
    success_count = 0
    for traj in inttrajs:
        loc = (-10,-10)
        for timestep in traj:
            prevloc = loc
            bin_cv1 = int((timestep[0]-grid.min_cv1)/grid.size_cv1)
            bin_cv2 = int((timestep[1]-grid.min_cv2)/grid.size_cv2)
            loc = (bin_cv1,bin_cv2)
            if loc in basin_A_set:
                basin_count += 1
                break
            if loc in basin_B_set:
                success_count += 1 
                break
            if loc not in active_set:
                cross_count += 1
                first_crosses.append(timestep)
                # And track where it exited from
                if prevloc in edges:
                    if prevloc in edgecount:
                        edgecount[prevloc] += 1
                    else:
                        edgecount[prevloc] = 1
                break

    # Now we average the edge counts to smooth data a bit
    # Average edge counts with edges +-3 sites 
    # in cv1/cv2 directions
    avg_edgecount = {}
    for site in edges:
        bin_cv1 = site[0]
        bin_cv2 = site[1]
        avg_edgecount[site] = 0
        count = 0
        for d_cv1 in range(-3,4):
            for d_cv2 in range(-3,4):
                loc = (bin_cv1+d_cv1,bin_cv2+d_cv2)
                if loc in edges:
                    if loc in edgecount:
                        avg_edgecount[site] += edgecount[loc]
                    count +=1 
        avg_edgecount[site] /= count

    return first_crosses, cross_count, basin_count, success_count, avg_edgecount

def calc_exit_final(inttrajs,grid,basin_A_set,basin_B_set):
    basin_count = 0
    success_count = 0
    for traj in inttrajs:
        loc = (-10,-10)
        for timestep in traj:
            prevloc = loc
            bin_cv1 = int((timestep[0]-grid.min_cv1)/grid.size_cv1)
            bin_cv2 = int((timestep[1]-grid.min_cv2)/grid.size_cv2)
            loc = (bin_cv1,bin_cv2)
            if loc in basin_A_set:
                basin_count += 1
                break
            if loc in basin_B_set:
                success_count += 1 
                break

    return basin_count, success_count

