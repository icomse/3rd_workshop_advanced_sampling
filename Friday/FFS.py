# FFS.py
# Created by PH Minh and Naomi Trampe
# Last Updated Date: 03-18-2023 
# SAMPEL Group

# Import necessary packages
import numpy as np
import sys
import math
import langevin_dynamics as ld
import random
import matplotlib.pyplot as plt
from matplotlib import colors, ticker, cm
import time
from matplotlib.colors import Normalize 
from scipy.interpolate import interpn

def performffs(pes_type,basin_pos,place,op_type_list,interfaces,basinlen,init_coords,init_p,dt,beta,gamma,lag,interface_trajs,n_confs=30,explore_frac=0.2,p_des=0.3,d_min=0.1):
    print("performing FFS...", end="\r")
    op_xtype = op_type_list[0]
    op_ytype = op_type_list[1]
    op_xcoef = op_type_list[2]
    op_ycoef = op_type_list[3]
    well=4
    basineqlen = 5000
    space = np.linspace(-3.5,3.5,1000)
    sp_lim = max(space)-min(space)
    if place:
        interfaces = []
        basinB = basin_pos[0]
        basinB_inter = [[],[]]
        for i in range(space.shape[0]):
            print(f"Placing B... {i+1}/{space.shape[0]}", end="\r")
            for j in range(space.shape[0]):
                if -(sp_lim/space.shape[0]/2) < (abs(ld.calc_op_f(op_xtype,op_ytype,op_xcoef,op_ycoef,space[i],space[j])))**(1/max(op_xtype,op_ytype))-((abs(basinB))**(1/max(op_xtype,op_ytype))) < (sp_lim/space.shape[0]/2):
                    if ld.calc_op_f(op_xtype,op_ytype,op_xcoef,op_ycoef,space[i],space[j])*basinB > 0:
                        basinB_inter[0].append(space[i])
                        basinB_inter[1].append(space[j])
    else:
        basinB = basin_pos[0]
        basinA = basin_pos[1]
    print("Equilibrating...          ", end="\r")
    # declare array to store basin trajectory
    basintraj = np.zeros((basinlen + 1, 6),dtype=float)
    # calculate initial forces
    fx,fy = ld.force(init_coords[0],init_coords[1],init_p[0],init_p[1],dt,beta,gamma,pes_type)
    # combine positions, momenta, and forces to make an initial phase point
    init_phasepoint = init_coords + init_p + [fx,fy]
    basintrajeq = ld.vv_step(init_phasepoint,dt,beta,gamma,pes_type)
    # equilibrate in basin
    eqtraj = []
    for i in range(1,basineqlen + 1):
        new_basintrajeq = ld.vv_step(basintrajeq,dt,beta,gamma,pes_type)
        basintrajeq = new_basintrajeq
        eqtraj.append(basintrajeq)
        op = ld.calc_op_f(op_xtype,op_ytype,op_xcoef,op_ycoef,basintrajeq[0],basintrajeq[1])
        # check if trajectory reaches basin B
        if op >= basinB:
            sys.exit("Basin trajectory reached B! Exiting...")
    # start at final point of equilibration
    basintraj[0] = basintrajeq
    fromBasin = False
    first_crosses = []
    harvest_crosses = []
    n_cross = 0
    n_harvest = 0
    cross_time = []
    harvest_time = []
    inter = []
    last = -lag
    print("Running basin simulation...       ", end='\r')
    # run basin A simulation
    ## check for first crossings if not placing interfaces
    for j in range(1,basinlen + 1):
        basintraj[j] = ld.vv_step(basintraj[j-1],dt,beta,gamma,pes_type)
        op = ld.calc_op_f(op_xtype,op_ytype,op_xcoef,op_ycoef,basintraj[j][0],basintraj[j][1])
        if not place:
            # collect first crossings|
            if op < basinA:
                fromBasin = True
            if fromBasin == True and op >= interfaces[0]:
                first_crosses.append(basintraj[j])
                n_cross += 1
                if j - last > lag:
                    harvest_crosses.append(basintraj[j])
                    n_harvest += 1
                    last = j
                fromBasin = False
            cross_time.append(n_cross)
            harvest_time.append(n_harvest)
        # check if trajectory reaches basin B
        if op >= basinB:
            sys.exit("Basin trajectory reached B! Exiting...")

    # if we are placing interfaces, we need to 
    # define basin A and the first interface
    # based on the basin trajectory
    if place:
        lmda = [ld.calc_op_f(op_xtype,op_ytype,op_xcoef,op_ycoef,basintraj[j][0],basintraj[j][1]) for j in range(len(basintraj))]
        avg = np.mean(lmda)
        stdev = np.std(lmda)
        # define basin A
        basinA = avg + 0.5*stdev
        #print("Basin A is at {}\n".format(basinA))
        basinA_inter = [[],[]]
        for i in range(space.shape[0]):
            print(f"Placing A... {i+1}/{space.shape[0]}               ", end="\r")
            for j in range(space.shape[0]):
                if -(sp_lim/space.shape[0]/2) < (abs(ld.calc_op_f(op_xtype,op_ytype,op_xcoef,op_ycoef,space[i],space[j])))**(1/max(op_xtype,op_ytype))-((abs(basinA))**(1/max(op_xtype,op_ytype))) < (sp_lim/space.shape[0]/2):
                    if ld.calc_op_f(op_xtype,op_ytype,op_xcoef,op_ycoef,space[i],space[j])*basinA > 0:
                        basinA_inter[0].append(space[i])
                        basinA_inter[1].append(space[j])
        interfaces = [0]
        done = False
        print("Finding interface 0...             ", end='\r')
        # define l0 to satisfy desired harvested crossings
        for i in np.linspace(basinA,max(lmda),int(200*abs(basinA-max(lmda)))):
            n_cross = 0
            n_harvest = 0
            fromBasin = False
            last = -lag
            for j in range(len(lmda)):
                if lmda[j] < basinA:
                    fromBasin = True
                if fromBasin == True and lmda[j] >= i:
                    n_cross += 1
                    fromBasin = False
                    if j - last > lag:
                        n_harvest += 1
                        last = j
            if n_harvest >= n_confs:
                interfaces[0] = i
                done = True
        # exit if not enough crossings are harvested
        if not done:
            sys.exit("Not enough harvested crossings. Exiting...")
        #print("Interface 0 is at {}\n".format(interfaces[0]))          
        inter_i = [[],[]]
        for i in range(space.shape[0]):
            print(f"Placing 0... {i+1}/{space.shape[0]}               ", end="\r")
            for j in range(space.shape[0]):
                if -(sp_lim/space.shape[0]/2) < (abs(ld.calc_op_f(op_xtype,op_ytype,op_xcoef,op_ycoef,space[i],space[j])))**(1/max(op_xtype,op_ytype))-((abs(interfaces[0]))**(1/max(op_xtype,op_ytype))) < (sp_lim/space.shape[0]/2):
                    if ld.calc_op_f(op_xtype,op_ytype,op_xcoef,op_ycoef,space[i],space[j])*interfaces[0] > 0:
                        inter_i[0].append(space[i])
                        inter_i[1].append(space[j])
        inter.append(inter_i)
        # collect first crossings
        fromBasin = False
        n_cross = 0
        n_harvest = 0
        last = -lag
        for j in range(len(lmda)):
            if lmda[j] < basinA:
                fromBasin = True
            if fromBasin == True and lmda[j] >= interfaces[0]:
                first_crosses.append(basintraj[j])
                n_cross += 1
                if j - last > lag:
                    harvest_crosses.append(basintraj[j])
                    n_harvest += 1
                    last = j
                fromBasin = False
            cross_time.append(n_cross)
            harvest_time.append(n_harvest)
    # exit if there are no first crossings
    if n_cross == 0:
        sys.exit("No first crossings obtained from basin A to interface 0. Exiting...")
    flux = n_cross/(basinlen*dt)
#    print("Flux through first interface: {}\n".format(flux))
    # evenly distribute configurations through simulations
#    print("Number of first crossings: {}\n".format(len(first_crosses)))
    configs = []
    configs_list = []
    configs_group = []
    configs_groups = []
    for i in range(len(first_crosses)):
        configs_group.append(0)
    configs_groups.append(configs_group)
    alltraj_group = []
    alltraj_groups = []
    for j in range(len(first_crosses)):
        for i in range(int(interface_trajs/len(first_crosses))):
            configs.append(first_crosses[j])
            alltraj_group.append(j)
    # randomly select the remainder without replacement
    for i in np.asarray(random.sample(range(len(first_crosses)),k=interface_trajs%len(first_crosses))):
        configs.append(first_crosses[i])
        alltraj_group.append(i)
    configs = np.asarray(configs)
    configs_list.append(configs)
    alltraj_groups.append(alltraj_group)
    # run interfaces sequentially
    starttime = time.monotonic()
    alltrajs = []
    cross_probs = []
    allyestrajs = []
    allnotrajs = []
    if not place:
        for i in range(len(interfaces) - 1):
            print(f"Starting interface {i}...               ", end="\r")
            inttrajs = []
            yestrajs = []
            notrajs = []
            configs_group = []
#            print("Starting interface {}...".format(i))
            first_crosses = []
            n_cross = 0
            # run simulations for each selected configuration
            for config in range(len(configs)):
                op = ld.calc_op_f(op_xtype,op_ytype,op_xcoef,op_ycoef,configs[config][0],configs[config][1])
                step = 0
                traj = []
                traj.append(configs[config])
                while op >= basinA and op < interfaces[i+1]:
                    traj.append(ld.vv_step(traj[step],dt,beta,gamma,pes_type,well))
                    step += 1
                    op = ld.calc_op_f(op_xtype,op_ytype,op_xcoef,op_ycoef,traj[step][0],traj[step][1])
                if op >= interfaces[i+1]:
                    n_cross += 1
                    first_crosses.append(traj[step])
                    configs_group.append(config)
                    yestrajs.append(np.asarray(traj))
                else:
                    notrajs.append(np.asarray(traj))
                inttrajs.append(np.asarray(traj))
            if n_cross == 0:
                sys.exit("No first crossings obtained from interface {} to interface {}. Exiting...".format(i,i+1))
            print(f"Finishing interface {i}...               ", end="\r")
            # determine the crossing probability
            cross_prob = n_cross/(len(inttrajs))
            cross_probs.append(cross_prob)
            alltrajs.append(inttrajs)
            allyestrajs.append(yestrajs)
            allnotrajs.append(notrajs)
#            print("Interface {} first crossings from {}: {}".format(i+1,i,len(first_crosses)))
#            print("{} to {} crossing prob: {}\n".format(i,i+1,cross_prob))
            configs = []
            alltraj_group = []
            for j in range(len(first_crosses)):
                for i in range(int(interface_trajs/len(first_crosses))):
                    configs.append(first_crosses[j])
                    alltraj_group.append(j)
            # randomly select the remainder without replacement
            for i in np.asarray(random.sample(range(len(first_crosses)),k=interface_trajs%len(first_crosses))):
                configs.append(first_crosses[i])
                alltraj_group.append(i)
            configs = np.asarray(configs)
            configs_list.append(configs)
            configs_groups.append(configs_group)
            alltraj_groups.append(alltraj_group)
    else:
        cnt = 0
        while interfaces[-1] < basinB:
            print(f"Exploring from interface {cnt}...               ", end="\r")
            ##EXPLORING SCOUTS##
#            print("Exploring from interface {}...".format(cnt))
            op_max = []
            for config in configs[::int(1/explore_frac)]:
                op = ld.calc_op_f(op_xtype,op_ytype,op_xcoef,op_ycoef,config[0],config[1])
                op_maxi = op
                step = 0
                traj = []
                traj.append(config)
                while op >= basinA and op < basinB and step < 5000:
                    traj.append(ld.vv_step(traj[step],dt,beta,gamma,pes_type,well))
                    step += 1
                    op = ld.calc_op_f(op_xtype,op_ytype,op_xcoef,op_ycoef,traj[step][0],traj[step][1])
                    if op > op_maxi:
                        op_maxi = op
                op_max.append(op_maxi)
            op_max.sort()
            if op_max[int((1-p_des)*len(op_max))] > basinB:
                interfaces.append(basinB)
            elif op_max[int((1-p_des)*len(op_max))]-interfaces[-1] < d_min:
                interfaces.append(interfaces[-1]+d_min)
            else:
                interfaces.append(op_max[int((1-p_des)*len(op_max))])
#            print("Next interface is at {}\n".format(interfaces[-1]))
            inter_i = [[],[]]
            for i in range(space.shape[0]):
                print(f"Placing {cnt+1}... {i+1}/{space.shape[0]}                               ", end="\r")
                for j in range(space.shape[0]):
                    if -(sp_lim/space.shape[0]/2) < (abs(ld.calc_op_f(op_xtype,op_ytype,op_xcoef,op_ycoef,space[i],space[j])))**(1/max(op_xtype,op_ytype))-((abs(interfaces[cnt+1]))**(1/max(op_xtype,op_ytype))) < (sp_lim/space.shape[0]/2):
                        if ld.calc_op_f(op_xtype,op_ytype,op_xcoef,op_ycoef,space[i],space[j])*interfaces[cnt+1] > 0:
                            inter_i[0].append(space[i])
                            inter_i[1].append(space[j])
            inter.append(inter_i)
            inttrajs = []
            yestrajs = []
            notrajs = []
            configs_group = []
            print(f"Starting interface {cnt}...               ", end="\r")
            ##INTERFACE SIMULATION##
#            print("Starting interface {}...".format(cnt))
            first_crosses = []
            n_cross = 0
            for config in range(len(configs)):
                op = ld.calc_op_f(op_xtype,op_ytype,op_xcoef,op_ycoef,configs[config][0],configs[config][1])
                step = 0
                traj = []
                traj.append(configs[config])
                while op >= basinA and op < interfaces[cnt+1]:
                    traj.append(ld.vv_step(traj[step],dt,beta,gamma,pes_type,well))
                    step += 1
                    op = ld.calc_op_f(op_xtype,op_ytype,op_xcoef,op_ycoef,traj[step][0],traj[step][1])
                if op >= interfaces[cnt+1]:
                    n_cross += 1
                    first_crosses.append(traj[step])
                    configs_group.append(config)
                    yestrajs.append(np.asarray(traj))
                else:
                    notrajs.append(np.asarray(traj))
                inttrajs.append(np.asarray(traj))
            if n_cross == 0:
                sys.exit("No first crossings obtained from interface {} to interface {}. Exiting...".format(cnt,cnt+1))
            cross_prob = n_cross/(len(inttrajs))
            cross_probs.append(cross_prob)
            alltrajs.append(inttrajs)
            allyestrajs.append(yestrajs)
            allnotrajs.append(notrajs)
            print(f"Finishing interface {cnt}...               ", end="\r")
#            print("Interface {} first crossings from {}: {}".format(cnt+1,cnt,len(first_crosses)))
#            print("{} to {} crossing prob: {}\n\n".format(cnt,cnt+1,cross_prob))
            configs = []
            alltraj_group = []
            for j in range(len(first_crosses)):
                for i in range(int(interface_trajs/len(first_crosses))):
                    configs.append(first_crosses[j])
                    alltraj_group.append(j)
            # randomly select the remainder without replacement
            for i in np.asarray(random.sample(range(len(first_crosses)),k=interface_trajs%len(first_crosses))):
                configs.append(first_crosses[i])
                alltraj_group.append(i)
            configs = np.asarray(configs)
            configs_list.append(configs)
            configs_groups.append(configs_group)
            alltraj_groups.append(alltraj_group)
            cnt += 1
            if cnt > 30:
                sys.exit("Too many interfaces placed. OP may not be able to sample from A to B. Exiting...")
    endtime = time.monotonic()
#    print("Total Simulation Time: %.2f s" % (endtime-starttime))
        
    rate = flux*np.prod(np.asarray(cross_probs))
#    print('Rate = %8.3e' % rate)
    print("Stitching the TPE...                ", end='\r')
    TPE = []
    for i in range(len(configs_groups[-1])): 
        complete_path = []
        prev_traj = configs_groups[-1][i]
        complete_path.append(alltrajs[-1][prev_traj])
        for j in range(2,len(interfaces)):
            conf = alltraj_groups[-j][prev_traj]
            prev_traj = configs_groups[-j][conf]
            complete_path.append(alltrajs[-j][prev_traj])
        TPE.append(complete_path)
    
    
    
    print("Plotting denisty...                     ", end='\r')
    skip = 1
    x = []
    y = []
    for i in range(len(TPE)):
        for j in range(len(TPE[i])):
            for k in TPE[i][j][:,0]:
                x.append(k)
            for k in TPE[i][j][:,1]:
                y.append(k)
    y = np.asarray(y)
    x = np.asarray(x)


    bins = 20
    sort = True
    # Plot potential energy surface contours
    N = 100
    x_vec = np.linspace(-3.5, 3.5, N)
    y_vec = np.linspace(-3.5, 3.5, N)
    X, Y = np.meshgrid(x_vec, y_vec)
    energy = np.zeros((N, N))

    # Plot contours
    for i in range(len(x_vec)):
        for j in range(len(y_vec)):
            energy[j][i] = ld.potential(pes_type,x_vec[i],y_vec[j])
    fig , ax = plt.subplots(figsize=(11.25,9))
    ax.contour(x_vec,y_vec,energy,np.linspace(-3,3,15), colors='silver')
    data , x_e, y_e = np.histogram2d( x, y, bins = bins, density = True )
    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , np.vstack([x,y]).T , method = "splinef2d", bounds_error = False)

    #To be sure to plot all data
    z[np.where(np.isnan(z))] = 0.0

    # Sort the points by density, so that the densest points are plotted last
    if sort :
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]

    m = 'jet'
    ax.scatter( x, y, c=z, s=2, cmap = m, alpha = 0.002)

    norm = Normalize(vmin = np.min(z), vmax = np.max(z))

    cbar = fig.colorbar(cm.ScalarMappable(norm = norm, cmap = m), ax=ax)
    cbar.set_ticks([])
    cbar.ax.set_ylabel('Density',size=12)
    ax.set_xlabel('x',fontsize=15)
    ax.set_ylabel('y',fontsize=15)
    ax.set_xlim(-3.5,3.5)
    ax.set_ylim(-3.5,3.5)
    ax.tick_params(axis='both',labelsize=12)
    if place:
        for i in range(len(inter)):
            ax.plot(inter[i][0],inter[i][1],'ko',markersize=0.25)
        ax.plot(basinA_inter[0],basinA_inter[1],'ro',markersize=0.25)
        ax.plot(basinB_inter[0],basinB_inter[1],'bo',markersize=0.25)
    else:
        if op_xtype == 1 and op_ycoef == [0]:
            for i in interfaces:
                ax.axvline(x=i,color='k',linewidth=1.5,alpha=0.5)
            ax.axvline(x=basinB,color='b',linewidth=1.5,alpha=0.5)
            ax.axvline(x=basinA,color='r',linewidth=1.5,alpha=0.5)
        elif op_ytype == 1 and op_xcoef == [0]:
            for i in interfaces:
                ax.axhline(y=i,color='k',linewidth=1.5,alpha=0.5)
            ax.axhline(y=basinB,color='b',linewidth=1.5,alpha=0.5)
            ax.axhline(y=basinA,color='r',linewidth=1.5,alpha=0.5)
        elif op_xtype == 1 and op_ytype == 1:
            x = np.linspace(-3.5,3.5,1000)
            for i in interfaces:
                y = i - x*op_xcoef/op_ycoef
                ax.plot(x,y,color='k',linewidth=1.5,alpha=0.5)
            y = basinB - x*op_xcoef/op_ycoef
            ax.plot(x,y,color='b',linewidth=1.5,alpha=0.5)
            y = basinA - x*op_xcoef/op_ycoef
            ax.plot(x,y,color='r',linewidth=1.5,alpha=0.5)
    plt.show()

          
          
    return print("""\
          
            🎉 CONGRATULATIONS 🎉
                    
                """)