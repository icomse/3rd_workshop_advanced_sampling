import numpy as np
import sys
import re
import subprocess
import MDAnalysis as mda
import copy
import dask
from dask.distributed import Client, LocalCluster

with dask.config.set({"distributed.worker.resources.cores": 12}):
    cluster = LocalCluster(n_workers=1, threads_per_worker=12, processes=False)
client = Client(cluster)


#NOTE: Swapping move is not enabled right now! -- See line 108 
      # Changed lines 174 if you'd like to ignore pathlength limit ( do --> alwaysAccept = True) -- we want to ignore pathlength because max length right now is determined by the length of the current path's shooting slice (and since we only have 1 config right now, we would have a really small numerator when calculating the max length). As such, we would have a high rejection chance if we limit pathlength because we would keep having short paths and they will get rejected before we get an AA or AB path.  
      # TIME REVERSALS CURRENTLY DISABLED and PATH LENGTH CRITERIA IS TURNED OFF



# DEFINE PATHS HERE
MPATH = "/ocean/projects/see220002p/$(whoami)/LiF-RETIS/run-retis/"
SHPATH_TRJCONV = MPATH + "scripts/trjconv.sh"
SHPATH_SIMULATE = MPATH + "scripts/simulate.sh"
SHPATH_EXTEND = MPATH + "scripts/extend.sh"
GROPATH = MPATH + "masters/npteq_LiF.gro"
INMASTERPATH = MPATH + "masters/npt.mdp"
TOPPATH = MPATH + "masters/LiF.top"
LOG_PATH = "RETIS.log"
FLOG = open(LOG_PATH,'w')

# MAKE SOME OTHER GLOBAL DEFS HERE
# current code only works with commit length 1
commit_length = 1 # number of frames in basin before being classified as returned to basin


interfaces = [ 1.9199, 2.1319, 2.1999, 2.2599,
               2.2799, 2.2999, 2.3099, 2.3399,
               2.3799, 2.3999, 2.4199, 2.4299, 2.4999,
               2.5999, 4.0000]

# number of steps to run before checking for crossing
nsteps = [ 50, 50, 50, 50, 50, 50, 50, 2000, 2000, 2000, 2000, 2000, 2000, 2000 ]

# Quick sanity check
ninterfaces = len(interfaces)-1
assert len(nsteps) == ninterfaces, \
        "Error, nsteps not specified for each interface" 

# Global vars to store paths
paths = [None] * (ninterfaces)
pathlogs = [None] * (ninterfaces)
basin_path = None

def main():
    global interfaces
    global ninterfaces
    global paths
   
    # Read in paths from starting move
    read_paths('initf')

    # Now we actually start RETIS
    for glmoveid in range(0,51):
        print('Move {}'.format(glmoveid))
        global_move_generator(glmoveid)

def read_paths(glmoveid):
    global MPATH
    global paths
    global basin_path
    # Read in interface paths
    for interface in range(ninterfaces):
        filen = MPATH + str(interface) + "/paths/path_gl-" + str(glmoveid) + \
                "_ens-" + str(interface) + ".txt"
        path = np.genfromtxt(filen,comments=None,dtype=[('pathid','U50'),('time','f8'),('op','f8'),('rev','i8')],encoding='utf8')
        paths[interface] = path
    # And read in basin
    filen = MPATH + "basin/paths/path_gl-" + str(glmoveid) + \
            "_ens-basin" + ".txt"
    path = np.genfromtxt(filen,comments=None,dtype=[('pathid','U50'),('time','f8'),('op','f8'),('rev','i8')],encoding='utf8')
    basin_path = path

def write_paths(glmoveid):
    global MPATH
    global paths
    global basin_path
    for interface in range(ninterfaces):
        filen = MPATH + str(interface) + "/paths/path_gl-" + str(glmoveid) + \
                "_ens-" + str(interface) + ".txt"
        np.savetxt(filen,paths[interface],fmt='%-25s%20f%20f%10d')
    # And write basin
    filen = MPATH + "basin/paths/path_gl-" + str(glmoveid) + \
            "_ens-basin" + ".txt"
    np.savetxt(filen,basin_path,fmt='%-25s%20f%20f%10d')

def global_move_generator(glmoveid):
    global interfaces
    global ninterfaces
    global paths
    global pathlogs
    global basin_path
    FLOG.write("Global move {}: ".format(glmoveid)) 
    # Shooting or swapping?
    swap = True if np.random.ranf() < -0.5 else False # if -0.5, this will always be false, which means, swap is turned off
    if swap:
        FLOG.write("Swap: ")
        # Which ensembles to swap?
        zero_minus_ensembles = True if np.random.ranf() < 0.5 else False
        if zero_minus_ensembles:
            print('performing zero minus')
            FLOG.write("-0<-->0,1<-->2,3<-->4...\n")
            swap_status = client.submit(perform_zerominus_move, glmoveid, resources={"cores": 12})
            tmpstat = swap_status.result()
            if tmpstat[0] == True:
                paths[0] = tmpstat[1]
                basin_path = tmpstat[2]
                FLOG.write("\t[-0]<-->[0]: Success\n")
            elif tmpstat[0] == False:
                print("Error: zero minus move failed")
                exit(1)
                FLOG.write("\t[-0]<-->[0]: Fail\n")
            for i in range(1,ninterfaces-1,2):
                FLOG.write("\t[{}]<-->[{}] Attempt\n".format(i,i+1))
                swap_status = attempt_swap(i)
                if swap_status:
                    FLOG.write("\t[{}]<-->[{}]: Success\n".format(i,i+1))
                else:
                    FLOG.write("\t[{}]<-->[{}]: Fail\n".format(i,i+1))
        else:
            FLOG.write("0<-->1,2<-->3,4<-->5...\n")
            for i in range(0,ninterfaces-1,2):
                FLOG.write("\t[{}]<-->[{}] Attempt\n".format(i,i+1))
                swap_status = attempt_swap(i)
                if swap_status:
                    FLOG.write("\t[{}]<-->[{}]: Success\n".format(i,i+1))
                else:
                    FLOG.write("\t[{}]<-->[{}]: Fail\n".format(i,i+1))
    else:
        print('starting shooting moves')
        FLOG.write("\n")
        for i in range(ninterfaces):
            path = paths[i]
            paths[i] = client.submit(ensemble_move_generator, i, path, glmoveid, pure=False, resources={"cores": 12})
        tmppaths = [paths[j].result() for j in range(len(paths))]
        paths = [tmppaths[j][0] for j in range(len(tmppaths))]
        logs = [tmppaths[j][1] for j in range(len(tmppaths))]
        for i in range(ninterfaces):
            FLOG.write(logs[i])
    # After each global move record the
    # paths for each ensemble
    FLOG.flush()
    write_paths(glmoveid)

def attempt_swap(interface):
    global paths 
    # Need to know structure of paths
    # If i path meets requirements of i+1
    # then accept. Else reject and re-count
    nslice = check_path(paths[interface],interface+1)
    if nslice is not None:
        tmp_i = paths[interface]
        tmp_ni = paths[interface+1]
        paths[interface] = tmp_ni
        paths[interface+1] = tmp_i
        return True
    else:
        return False

def ensemble_move_generator(interface,path,glmoveid):
    moveid = "gl-" + str(glmoveid) + "_ens-" + str(interface)
    # Shooting or time reversal? 
    if np.random.ranf() < 1.5: # we're not doing a reversal move, because this if is always going to be true
        path = perform_shooting_move(interface,path,moveid,alwaysAccept=True)
    else:
        path = perform_reversal_move(path)
    return path

def perform_reversal_move(path):
    global commit_length
    global interfaces
    trev = np.flip(path,axis=0)
    trev['rev'] = -1*trev['rev']
    # Verify that traj does not have either end in B
    basinB_frames = np.where(trev['op'] > interfaces[-1])[0]
    if basinB_frames.shape[0] == 0:
        #FLOG.write("Success\n")
        return trev
    else:
        #FLOG.write("Reject\n")
        return path

def perform_zerominus_move(glmoveid):
    global paths
    global basin_path
    global commit_length
    global nsteps
    zero_path = paths[0]
    # First run minus section of path
    interface="basin"
    moveid = "gl-" + str(glmoveid) + "_ens-" + interface
    # Harvest initial frame of 0+ path
    ntime = zero_path['time'][0]
    trajname = zero_path['pathid'][0]
    if zero_path['rev'][0] == -1:
        reverse = True
    else:
        reverse = False

    # Harvest gro file for shooting point
    call_trjconv(interface,trajname,ntime,moveid,reverse,zeroswap=True)
    # Run simulation/analysis
    nstep = nsteps[0]
    call_backward_simulation(interface,moveid,nstep)
    # Extract output from op file
    path_seg = extract_path_minus(interface,moveid,0,None)
    # Move rejected if path is too short
    if path_seg is None:
        print("Error: empty zerominus path.\n")
        exit(1)
        return False

    # If not none, we have a new basin path
    new_basin_path = np.hstack((path_seg,zero_path[1]))
   
    # harvest .gro from end of old 
    # 0- path and use to start simulation in the 0+ ensemble
    interface=0
    moveid = "gl-" + str(glmoveid) + "_ens-" + str(interface)

    ntime = basin_path['time'][-1]
    trajname = basin_path['pathid'][-1]
    if zero_path['rev'][0] == -1:
        reverse = True
    else:
        reverse = False

    # Harvest gro file for shooting point
    call_trjconv(interface,trajname,ntime,moveid,reverse,zeroswap=True)
    # Run simulation/analysis
    call_simulation(interface,moveid,nstep)
    # Extract output from op file
    path_seg = extract_path_zero(interface,moveid,0,None)

    # Combine segment with old path portion (from end of basin path)
    new_zero_path = np.hstack((basin_path[-2:],path_seg))

    return True, new_zero_path, new_basin_path
 
def call_trjconv(interface,trajname,ntime,moveid,reverse,zeroswap=False):
    global SHPATH_TRJCONV
    global GROPATH
    global MPATH
    #global beta
    #global temperature
    #sigma = np.sqrt(temperature)*10
    epath = MPATH + str(interface) + "/"
    traj_loc = re.sub("_rev","",re.sub("_extend-[0-9]*","",re.sub("gl-[\w]*ens-","",trajname)))

    #trjpath = epath + "/xtc/" + trajname + ".xtc"
    trjpath = MPATH + traj_loc + "/xtc/" + trajname + ".xtc"
    trrpath = MPATH + traj_loc + "/trr/" + trajname + ".trr"
    outpath = epath + "gro/" + moveid + ".gro"

    command = "source " + SHPATH_TRJCONV + " " + epath + " " + \
              GROPATH + " " + trrpath + " " + str(ntime) + " " + outpath
    subprocess.call(['bash','-c', command])
    # The command above calls trjconv script, which harvest a gro file and copy the gro file to an _old.gro file so we can call it in the if statement below

    if zeroswap == False:
        oldgro = epath + "gro/" + moveid + "_old.gro"
        u = mda.Universe(oldgro)
        velocityshape = np.shape(u.atoms.velocities)
        
        if reverse == True:
            if 'gl-0' in trrpath or 'gl-0' in oldgro or 'gl-initf' in oldgro or 'gl-initf' in trrpath:   
                u.atoms.velocities = 1*u.atoms.velocities
            else:
                u.atoms.velocities = -1*u.atoms.velocities
        #new_velocities = np.random.normal(0,sigma,velocityshape)
        #u.atoms.velocities = copy.deepcopy(new_velocities)
        
        # write gro file with modified velocities
        outputgro = u.select_atoms('all')
        outputgro.write(outpath)
        return False
    else:
        oldgro = epath + "gro/" + moveid + "_old.gro"
        u = mda.Universe(oldgro)
        if reverse == True:
            if 'gl-0' in trrpath or 'gl-0' in oldgro or 'gl-initf' in oldgro or 'gl-initf' in trrpath:   
                u.atoms.velocities = 1*u.atoms.velocities
            else:
                u.atoms.velocities = -1*u.atoms.velocities
        outputgro = u.select_atoms('all')
        outputgro.write(outpath)

def call_backward_simulation(interface,moveid,nstep):
    global SHPATH_SIMULATE
    global MPATH
    global INMASTERPATH
    global TOPPATH
    epath = MPATH + str(interface) + "/"

    ingropath = epath + "gro/" + moveid + ".gro"
    infilepath = epath + "input/" + moveid + ".mdp"
    xtcfilepath = epath + "xtc/" + moveid + "_rev.xtc"
    trrfilepath = epath + "trr/" + moveid + "_rev.trr"
    opfilepath = epath + "op/" + moveid + "_rev.txt"
    tprpath = epath + "tpr/" + moveid + "_rev.tpr"
    edrpath = epath + "log/" + moveid + "_rev.edr"
    cptpath = epath + "rst/" + moveid + "_rev.cpt"

    # generate gro file with reversed velocities
    revgropath = epath + "gro/" + moveid + "_rev.gro"
    u = mda.Universe(ingropath)
    u.atoms.velocities = -1*u.atoms.velocities
    reversedgro = u.select_atoms('all')
    reversedgro.write(revgropath)
    randn1 = -1
 
    command = "source " + SHPATH_SIMULATE + " " + epath + " " + revgropath + " " + \
              INMASTERPATH + " " + str(randn1) + " " + \
              str(nstep) + " " + moveid + " " + infilepath + " " + xtcfilepath + " " + opfilepath + \
              " " + TOPPATH + " " + tprpath + " " + edrpath + " " + cptpath + " " + trrfilepath
    subprocess.call(['bash','-c', command])

def call_simulation(interface,moveid,nstep):
    global SHPATH_SIMULATE
    global MPATH
    global INMASTERPATH
    global TOPPATH
    epath = MPATH + str(interface) + "/"

    ingropath = epath + "gro/" + moveid + ".gro"
    infilepath = epath + "input/" + moveid + ".mdp"
    xtcfilepath = epath + "xtc/" + moveid + ".xtc"
    trrfilepath = epath + "trr/" + moveid + ".trr"
    opfilepath = epath + "op/" + moveid + ".txt"
    tprpath = epath + "tpr/" + moveid + ".tpr"
    edrpath = epath + "log/" + moveid + ".edr"
    cptpath = epath + "rst/" + moveid + ".cpt"
    randn1 = -1
 
    command = "source " + SHPATH_SIMULATE + " " + epath + " " + ingropath + " " + \
              INMASTERPATH + " " + str(randn1) + " " + \
              str(nstep) + " " + moveid + " " + infilepath + " " + xtcfilepath + " " + opfilepath + \
              " " + TOPPATH + " " + tprpath + " " + edrpath + " " + cptpath + " " + trrfilepath
    subprocess.call(['bash','-c', command])

def extend_simulation(interface,moveid,nstep,itr,backward=False):
    global SHPATH_EXTEND
    global MPATH
    epath = MPATH + str(interface) + "/"

    ingropath = epath + "gro/" + moveid + ".gro"
    if backward == True:
        if itr == 1:
            inrstpath = epath + "rst/" + moveid + "_rev.cpt"
            oldtprpath = epath + "tpr/" + moveid + "_rev.tpr"
        else:
            prev_itr = itr-1
            inrstpath = epath + "rst/" + moveid + "_extend-" + str(prev_itr) + "_rev.cpt"
            oldtprpath = epath + "tpr/" + moveid + "_extend-" + str(prev_itr) + "_rev.tpr"
        tprpath = epath + "tpr/" + moveid + "_extend-" + str(itr) + "_rev.tpr"
        xtcfilepath = epath + "xtc/" + moveid + "_extend-" + str(itr) + "_rev.xtc"
        trrfilepath = epath + "trr/" + moveid + "_extend-" + str(itr) + "_rev.trr"
        opfilepath = epath + "op/" + moveid + "_extend-" + str(itr) + "_rev.txt"
        edrpath = epath + "log/" + moveid + "_extend-" + str(itr) + "_rev.edr"
        cptpath = epath + "rst/" + moveid + "_extend-" + str(itr) + "_rev.cpt"
        fullmoveid = moveid + "_extend-" + str(itr) + "_rev"
    else:
        if itr == 1:
            inrstpath = epath + "rst/" + moveid + ".cpt"
            oldtprpath = epath + "tpr/" + moveid + ".tpr"
        else:
            prev_itr = itr-1
            inrstpath = epath + "rst/" + moveid + "_extend-" + str(prev_itr) + ".cpt"
            oldtprpath = epath + "tpr/" + moveid + "_extend-" + str(prev_itr) + ".tpr"
        tprpath = epath + "tpr/" + moveid + "_extend-" + str(itr) + ".tpr"
        xtcfilepath = epath + "xtc/" + moveid + "_extend-" + str(itr) + ".xtc"
        trrfilepath = epath + "trr/" + moveid + "_extend-" + str(itr) + ".trr"
        opfilepath = epath + "op/" + moveid + "_extend-" + str(itr) + ".txt"
        edrpath = epath + "log/" + moveid + "_extend-" + str(itr) + ".edr"
        cptpath = epath + "rst/" + moveid + "_extend-" + str(itr) + ".cpt"
        fullmoveid = moveid + "_extend-" + str(itr)
 
    command = "source " + SHPATH_EXTEND + " " + epath + " " + ingropath + " " + \
              inrstpath + " " + str(nstep) + " " + \
              fullmoveid + " " + xtcfilepath + " " + opfilepath + \
              " " + tprpath + " " + oldtprpath + " " + cptpath + " " + edrpath + " " + trrfilepath
    subprocess.call(['bash','-c', command])

def extract_path(interface,moveid,itr,partial_path,max_length,backward=False):
    global MPATH
    global commit_length

    if backward == True:
        if partial_path is None:
            opfilepath = MPATH + str(interface) + "/op/" + moveid + "_rev.txt"
    else:
        if partial_path is None:
            opfilepath = MPATH + str(interface) + "/op/" + moveid + ".txt"
    
    # Save new path section into temporary structure
    tmp_seg = np.genfromtxt(opfilepath,comments=None,dtype=[('time','f8'),('op','f8')])
    # Create datatype to include moveid/extend iteration
    new_dt = np.dtype([('pathid','U50')] + tmp_seg.dtype.descr + [('rev','i8')])
    new_seg = np.zeros(tmp_seg.shape,dtype=new_dt)
    # Copy info into new structure
    new_seg['time'] = tmp_seg['time']
    new_seg['op'] = tmp_seg['op']
    if backward == True:
        new_seg['rev'] = -1
    else:
        new_seg['rev'] = 1
    # Include moveid and append old partial path if applicable 
    if partial_path is None:
        if backward == True:
            new_seg['pathid'] = moveid + "_rev"
        else:
            new_seg['pathid'] = moveid
        path_seg = new_seg
    else:
        print("Error, simulation " + moveid + " had partial path where it should not.")
        exit(1)
        new_seg['pathid'] = moveid + "_extend-" + str(itr)
        path_seg = np.hstack((partial_path,new_seg))

    basinA_frames = np.where(path_seg['op'] < interfaces[0])[0]
    basinA_fr = scan_arr(basinA_frames,commit_length)
    basinB_frames = np.where(path_seg['op'] > interfaces[-1])[0]
    basinB_fr = scan_arr(basinB_frames,commit_length)

    # Verify that move completed without weirdness (this should not happen)
    if basinA_fr is not None and basinB_fr is not None:
        print("Error, simulation " + moveid + " in both basins A and B.")
        exit(1)

    # Reject if backward path does not reach basin A
    if basinA_fr is None and backward == True:
        reject = True
        return path_seg, len(path_seg), reject

    # If backward segmant reaches A, reverse and return path segment
    if basinA_fr is not None and backward == True:
        trev = np.flip(path_seg[:basinA_fr+1],axis=0)
        if len(trev) > max_length:
            reject = True
            return trev, len(trev), reject
        else:
            reject = False
            return trev, len(trev), reject

    # Check forward segment
    if basinA_fr is None and basinB_fr is None and backward == False:
        reject = True
        return path_seg, len(path_seg), reject

    if basinA_fr is not None and backward == False:
        if len(path_seg[:basinA_fr+1]) > max_length:
            reject = True
            return path_seg, len(path_seg), reject
        else:
            reject = False
            return path_seg[:basinA_fr+1], len(path_seg[:basinA_fr+1]), reject

    if basinB_fr is not None and backward == False:
        if len(path_seg[:basinB_fr+1]) > max_length:
            reject = True
            return path_seg, len(path_seg), reject
        else:
            reject = False
            return path_seg[:basinB_fr+1], len(path_seg[:basinB_fr+1]), reject

def extract_path_zero(interface,moveid,itr,partial_path):
    global MPATH
    global commit_length
    global nsteps

    nstep = nsteps[interface]
 
    if partial_path is None:
        opfilepath = MPATH + str(interface) + "/op/" + moveid + ".txt"
    else:
        opfilepath = MPATH + str(interface) + "/op/" + moveid + "_extend-" + str(itr) + ".txt"

    # Save new path section into temporary structure
    tmp_seg = np.genfromtxt(opfilepath,comments=None,dtype=[('time','f8'),('op','f8')])
    # Create datatype to include moveid/extend iteration
    new_dt = np.dtype([('pathid','U50')] + tmp_seg.dtype.descr + [('rev','i8')])
    new_seg = np.zeros(tmp_seg.shape,dtype=new_dt)
    # Copy info into new structure
    new_seg['time'] = tmp_seg['time']
    new_seg['op'] = tmp_seg['op']
    new_seg['rev'] = 1
    # Include moveid and append old partial path if applicable 
    # Use [:-1] to avoid duplicate frames when extending.
    if partial_path is None:
        new_seg['pathid'] = moveid
        path_seg = new_seg[:-1]
    else:
        new_seg['pathid'] = moveid + "_extend-" + str(itr)
        path_seg = np.hstack((partial_path,new_seg[:-1]))

    basinA_frames = np.where(path_seg['op'] < interfaces[0])[0]
    basinA_fr = scan_arr(basinA_frames,commit_length)
    basinB_frames = np.where(path_seg['op'] > interfaces[-1])[0]
    basinB_fr = scan_arr(basinB_frames,commit_length)

    # Verify that move completed without weirdness (this should not happen)
    if basinA_fr is not None and basinB_fr is not None:
        print("Error, simulation " + moveid + " in both basins A and B.")
        exit(1)

    # Here is where we extend simulations if req'd (YAY, recursive coding)
    if basinA_fr is None and basinB_fr is None:
        itr += 1
        extend_simulation(interface,moveid,nstep,itr,backward=False)
        path_seg = extract_path_zero(interface,moveid,itr,path_seg)
        # This is how we cascade out of the recursive section 
        return path_seg

    # Return new trajectory segment
    if basinB_fr is None:
        return path_seg[:basinA_fr+1]
    else:
        return path_seg[:basinB_fr+1]

def extract_path_minus(interface,moveid,itr,partial_path):
    global MPATH
    global commit_length
    global nsteps

    nstep = nsteps[0]
    
    if partial_path is None:
        opfilepath = MPATH + str(interface) + "/op/" + moveid + "_rev.txt"
    else:
        opfilepath = MPATH + str(interface) + "/op/" + moveid + "_extend-" + str(itr) + "_rev.txt"

    # Save new path section into temporary structure
    tmp_seg = np.genfromtxt(opfilepath,comments=None,dtype=[('time','f8'),('op','f8')])
    # Create datatype to include moveid/extend iteration
    new_dt = np.dtype([('pathid','U50')] + tmp_seg.dtype.descr + [('rev','i8')])
    new_seg = np.zeros(tmp_seg.shape,dtype=new_dt)
    # Copy info into new structure
    new_seg['time'] = tmp_seg['time']
    new_seg['op'] = tmp_seg['op']
    new_seg['rev'] = -1
    # Include moveid and append old partial path if applicable 
    if partial_path is None:
        new_seg['pathid'] = moveid + "_rev"
        path_seg = new_seg[:-1]
    else:
        new_seg['pathid'] = moveid + "_extend-" + str(itr) + "_rev"
        path_seg = np.hstack((partial_path,new_seg[:-1]))

    cross_zero = np.where(path_seg['op'] > interfaces[0])
    if cross_zero[0].shape[0] == 0:
        itr += 1
        extend_simulation(interface,moveid,nstep,itr,backward=True)
        path_seg = extract_path_minus(interface,moveid,itr,path_seg)
        # This is how we cascade out of the recursive section 
        return path_seg
    else:
        cross_zero_fr = cross_zero[0][0]
        if cross_zero_fr < commit_length:
            return None
        else:
            return np.flip(path_seg[:cross_zero_fr+1],axis=0)

def check_path(path,interface):
    slices = np.where(path['op'] > interfaces[interface])
    if slices[0].shape[0] == 0:
        return None
    else:
        nslice = slices[0][0]
    return nslice

def find_shoot_region(path,interface):
    if interface > 0 and interface < len(interfaces)-2:
        slices = np.where((path['op'] > interfaces[interface-1]) & (path['op'] < interfaces[interface+1]))
    elif interface == 0:
        slices = np.where((path['op'] < interfaces[interface+1]) & (path['op'] > interfaces[0]))
    elif interface == len(interfaces)-2:
        slices = np.where((path['op'] > interfaces[interface-1]) & (path['op'] < interfaces[-1]))


    if slices[0].shape[0] == 0:
        return None
    else:
        return slices[0]

def scan_arr(arr,n):
    for i in range(len(arr)-n+1):
        seq = arr[i:i+n] 
        if seq[n-1] == (seq[0] + (n-1)): 
            return seq[n-1] 
    return None

def perform_shooting_move(interface,path,moveid,alwaysAccept):
    global nsteps
    reject = False
    # select shooting point from between i-1 and i+1
    slices = find_shoot_region(path,interface)
    if slices is None:
        print("Error, no valid shooting point found for move " + moveid)
        #exit(1)
        #reject = True
        nslice = check_path(path,interface)
        pathlength = 1
    else:
        nslice = np.random.choice(slices)
        pathlength = len(slices) # length of current path
    if reject == True:
        #FLOG.write("Shoot {} {} Reject no ShootingPoint\n".format(interface,moveid))
        shootlog = "Shoot {} {} Reject no ShootingPoint\n".format(interface,moveid)
        return path, shootlog

    pathlength_randn = 0.0
    while pathlength_randn == 0.0:
        pathlength_randn = np.random.uniform()
    max_length = round((pathlength + 1)/pathlength_randn) # max length allowed for a path is current pathlength(length of path between i-1 and i+1) divide by the random number

    ntime = path[nslice][1]
    trajname = path[nslice][0]
    if path[nslice][3] == -1:
        reverse = True
    else:
        reverse = False
    
    # Harvest gro file for shooting point
    reject = call_trjconv(interface,trajname,ntime,moveid,reverse,zeroswap=False)
    if reject == True:
        #FLOG.write("Shoot {} {} Reject momenta\n".format(interface,moveid))
        shootlog = "Shoot {} {} Reject momenta\n".format(interface,moveid)
        return path, shootlog

    # Run reverse simulation/analysis
    nstep = nsteps[interface]
    call_backward_simulation(interface,moveid,nstep)

    backward_path_seg, reject, bpathlength = extract_path_shooting(interface,moveid,0,None,max_length,alwaysAccept,backward=True)
    if bpathlength > max_length and alwaysAccept == False:
        #FLOG.write("Shoot {} {} Reject pathlength\n".format(interface,moveid))
        shootlog = "Shoot {} {} Reject pathlength\n".format(interface,moveid)
        return path, shootlog
    if reject == True:
        #FLOG.write("Shoot {} {} Reject backward\n".format(interface,moveid))
        shootlog = "Shoot {} {} Reject backward\n".format(interface,moveid)
        return path, shootlog
    max_flength = max_length - bpathlength

    # Run forward sims
    call_simulation(interface,moveid,nstep)

    forward_path_seg, reject, fpathlength = extract_path_shooting(interface,moveid,0,None,max_flength,alwaysAccept,backward=False)
    if reject == True:
        #FLOG.write("Shoot {} {} Reject pathlength\n".format(interface,moveid))
        shootlog = "Shoot {} {} Reject pathlength\n".format(interface,moveid)
        return path, shootlog

    # Combine segment with old path portion
    new_path = np.hstack((backward_path_seg[:-1],forward_path_seg))

    # Check for crossing of interface i
    crossing = check_path(new_path,interface)
    if crossing is not None:
        # Check pathlength
        new_slices = find_shoot_region(new_path,interface)
        if new_slices is None:
            new_slices = [0]
        if len(new_slices) > max_length and alwaysAccept == False:
            #FLOG.write("Shoot {} {} Reject pathlength\n".format(interface,moveid))
            shootlog = "Shoot {} {} Reject pathlength\n".format(interface,moveid)
            return path, shootlog
        else:
            #FLOG.write("Shoot {} {} Accept crossing\n".format(interface,moveid))
            shootlog = "Shoot {} {} Accept crossing\n".format(interface,moveid)
            return new_path, shootlog
    else:
        #FLOG.write("Shoot {} {} Reject crossing\n".format(interface,moveid))
        shootlog = "Shoot {} {} Reject crossing\n".format(interface,moveid)
        return path, shootlog

def extract_path_shooting(interface,moveid,itr,partial_path,max_length,alwaysAccept,backward=False):
    global MPATH
    global commit_length
    global nsteps

    nstep = nsteps[interface]
 
    if backward == True:
        if partial_path is None:
            opfilepath = MPATH + str(interface) + "/op/" + moveid + "_rev.txt"
        else:
            opfilepath = MPATH + str(interface) + "/op/" + moveid + "_extend-" + str(itr) + "_rev.txt"
    else:
        if partial_path is None:
            opfilepath = MPATH + str(interface) + "/op/" + moveid + ".txt"
        else:
            opfilepath = MPATH + str(interface) + "/op/" + moveid + "_extend-" + str(itr) + ".txt"
    
    # Save new path section into temporary structure
    tmp_seg = np.genfromtxt(opfilepath,comments=None,dtype=[('time','f8'),('op','f8')])
    # Create datatype to include moveid/extend iteration
    new_dt = np.dtype([('pathid','U50')] + tmp_seg.dtype.descr + [('rev','i8')])
    new_seg = np.zeros(tmp_seg.shape,dtype=new_dt)
    # Copy info into new structure
    new_seg['time'] = tmp_seg['time']
    new_seg['op'] = tmp_seg['op']
    if backward == True:
        new_seg['rev'] = -1
    else:
        new_seg['rev'] = 1
    # Include moveid and append old partial path if applicable 
    if partial_path is None:
        if backward == True:
            new_seg['pathid'] = moveid + "_rev"
        else:
            new_seg['pathid'] = moveid
        path_seg = new_seg[:-1]
    else:
        if backward == True:
            new_seg['pathid'] = moveid + "_extend-" + str(itr) + "_rev"
        else:
            new_seg['pathid'] = moveid + "_extend-" + str(itr)
        path_seg = np.hstack((partial_path,new_seg[:-1]))

    basinA_frames = np.where(path_seg['op'] < interfaces[0])[0]
    basinA_fr = scan_arr(basinA_frames,commit_length)
    basinB_frames = np.where(path_seg['op'] > interfaces[-1])[0]
    basinB_fr = scan_arr(basinB_frames,commit_length)

    # Cut the path when it first reaches A or B
    if basinA_fr is not None:
        tmp_path_seg = path_seg[:basinA_fr+1]
    elif basinB_fr is not None:
        tmp_path_seg = path_seg[:basinB_fr+1]
    else:
        tmp_path_seg = path_seg

    # Find pathlength
    currentslices = find_shoot_region(tmp_path_seg,interface)
    if currentslices is not None and len(currentslices) > max_length and alwaysAccept == False:
        reject = True
        return tmp_path_seg, reject, len(currentslices)
    if currentslices is None:
        currentslices = [0]


    # Verify that move completed without weirdness (this should not happen)
    if basinA_fr is not None and basinB_fr is not None:
        print("Error, simulation " + moveid + " in both basins A and B.")
        exit(1)

    reject = False
    # Here is where we extend simulations if req'd (YAY, recursive coding)
    if basinA_fr is None and basinB_fr is None:
        itr += 1
        extend_simulation(interface,moveid,nstep,itr,backward)
        path_seg, reject, currentpathlength = extract_path_shooting(interface,moveid,itr,path_seg,max_length,alwaysAccept,backward)
        return path_seg, reject, len(currentslices)

    # Reject if backward segment goes to B, always accept forward segment
    if basinB_fr is None:
        if backward == True:
            trev = np.flip(path_seg[:basinA_fr+1],axis=0)
            return trev, reject, len(currentslices)
        else:
            return path_seg[:basinA_fr+1], reject, len(currentslices)
    else:
        if backward == True:
            reject = True
            return path_seg, reject, len(currentslices)
        else:
            return path_seg[:basinB_fr+1], reject, len(currentslices)


if __name__ == "__main__":
    main()

