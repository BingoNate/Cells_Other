
""" Read hdf5 data of cells and write some data
        IMPORTANT NOTE:
        All trajectory information is assumed UNWRAPPED!"""

### example command line arguments: 
###

##############################################################################

import argparse
import numpy as np
import os
import h5py
import misc_tools

##############################################################################
    
def read_h5_file(folder):
    """ read data from hdf5 file"""
    
    ### file path
    
    fpath = folder + 'out.h5'
    assert os.path.exists(fpath), "The out.h5 file does NOT exist for " + fpath
    fl = h5py.File(fpath, 'r')
    
    ### positions of beads
    
    xu = np.array(fl['/beads/xu'], dtype=np.float32)
    cid = np.array(fl['/beads/cid'], dtype=np.float32)
    
    ### cell information
    
    comu = np.array(fl['/cells/comu'], dtype=np.float32)
    pol = np.array(fl['/cells/pol'], dtype=np.float32)
    nbpc = np.array(fl['/cells/nbpc'], dtype=np.float32)
    
    ### simulation information
    
    lx = fl['/info/box/x'][...]
    ly = fl['/info/box/y'][...]
    dt = fl['/info/dt'][...]
    nsteps = fl['/info/nsteps'][...]
    ncells = fl['/info/ncells'][...]
    nbeads = fl['/info/nbeads'][...]
    nsamp = fl['/info/nsamp'][...]
    
    ### simulation parameters
    
    eps = fl['/param/eps'][...]
    rho = fl['/param/rho'][...]
    fp = fl['/param/fp'][...]
    areak = fl['/param/areak'][...]
    bl = fl['/param/bl'][...]
    sigma = fl['/param/sigma'][...]
    
    fl.close()
    
    ### generate classes to submerge data
    
    sim = misc_tools.Simulation(lx, ly, dt, nsteps, ncells, nbeads, nsamp, nbpc, \
                     eps, rho, fp, areak, bl, sigma)
    cells = misc_tools.Cells(comu, pol, nbpc, sim)
    beads = misc_tools.Beads(xu, cid)
    
    return sim, cells, beads

##############################################################################
    
def read_sim_info(folder):
    """ read simulation info from hdf5 file"""
    
    ### file path
    
    fpath = folder + 'out.h5'
    assert os.path.exists(fpath), "The out.h5 file does NOT exist for " + fpath
    fl = h5py.File(fpath, 'r')    
    
    ### cell information
    
    nbpc = np.array(fl['/cells/nbpc'], dtype=np.float32)
    
    ### simulation information
    
    lx = fl['/info/box/x'][...]
    ly = fl['/info/box/y'][...]
    dt = fl['/info/dt'][...]
    nsteps = fl['/info/nsteps'][...]
    ncells = fl['/info/ncells'][...]
    nbeads = fl['/info/nbeads'][...]
    nsamp = fl['/info/nsamp'][...]
    
    ### simulation parameters
    
    eps = fl['/param/eps'][...]
    rho = fl['/param/rho'][...]
    fp = fl['/param/fp'][...]
    areak = fl['/param/areak'][...]
    bl = fl['/param/bl'][...]
    sigma = fl['/param/sigma'][...]
    
    fl.close()
    
    ### generate classes to submerge data
    
    sim = misc_tools.Simulation(lx, ly, dt, nsteps, ncells, nbeads, nsamp, nbpc, \
                     eps, rho, fp, areak, bl, sigma)
    
    return sim
    
##############################################################################
   
def write_single_analysis_data(data, f):
    """ write single analysis data to the corresponding file"""

    fl = open(f, 'w')
    fl.write(str(data) + "\n")
    fl.close()
    
    return
 
##############################################################################
    
def write_2d_analysis_data(x, y, savebase, savefolder, sim):
    """ write 2d analysis data to the corresponding file
    EXAMPLE:
    savebase: /usr/users/iff_th2/duman/Cells_in_LAMMPS/DATA/
    savefolder: MSD"""
    
    ### create the path
    
    base = savebase + savefolder + '/'
    os.system("mkdir -p " + base)
    fpath = base + savefolder + "_eps_" + \
        str(sim.eps) + "_fp_" + str(sim.fp) + "_areak_" + str(sim.areak) + ".txt"  
  
    ### write the data
    
    fl = open(fpath, 'w')
    
    N = len(x)
    for j in range(N):
        fl.write(str(x[j]) + '\t\t' + str(y[j]) + '\n')
    
    fl.close()

    return
    
##############################################################################

def read_single_analysis_data(f):
    """ read single analysis data"""
    
    data = np.loadtxt(f, dtype=np.float64)

    return data  

##############################################################################    
    
def read_2d_analysis_data(f):
    """ read 2d analysis data"""
    
    data = np.transpose(np.loadtxt(f, dtype=np.float64))
    x = data[0]
    y = data[1]

    return x, y 

##############################################################################
    
def read_multid_analysis_data(f):
    """ read multid analysis data"""
    
    data = np.transpose(np.loadtxt(f, dtype=np.float64))
    
    return data

##############################################################################

def gen_folders(eps, fp, areak, analysis, dbase, analysisdbase):
    """ data and analysis folders generator"""
    
    path1 = 'eps_' + str(eps) + '/fp_' + str(fp) + '/areak_' + str(areak)   
    path2 = analysis + '_eps_' + str(eps) + '_fp_' + str(fp) + '_areak_' + str(areak) + '.txt'  
    datafolder = dbase + path1 + '/'
    analysisfile = analysisdbase + path2  

    return datafolder, analysisfile
    
##############################################################################
    
def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-fl", "--folder", help="Folder containing data")        
    args = parser.parse_args()    
    read_h5_file(args.folder)
    
    return
    
##############################################################################

if __name__ == '__main__':
    main()    
    
##############################################################################
