
""" Turn bead information to cell information and save the data to make things faster"""

##############################################################################

import argparse
import numpy as np
import os
import h5py
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import math

##############################################################################

def neigh_min(dx, lx):
    """ calculate minimum image distance between neighboring beads"""
    
    dx1 = dx + lx
    dx2 = dx - lx
    if dx**2 < dx1**2 and dx**2 < dx2**2:
        return dx
    if dx1**2 < dx2**2:
        return dx1
    return dx2 

##############################################################################
    
def neigh_min_array(dx, lx):
    """ calculate minimum image distance between neighboring beads for an array of displacements"""
    
    for j in range(len(dx)):
        dx[j] = neigh_min(dx[j], lx)
    
    return dx
    
##############################################################################        
   
class Simulation:
    """ data structure for storing general simulation information"""

    def __init__(self, lx, ly, dt, nsteps, nbeads, nsamp, ncells, nbpc, eps, rho, \
               fp, areak, bl, sigma):
        
        self.lx = lx
        self.ly = ly
        self.dt = dt
        self.nsteps = nsteps
        self.nbeads = nbeads
        self.nsamp = nsamp
        self.ncells = ncells
        self.nbpc = nbpc
        self.eps = eps
        self.rho = rho
        self.fp = fp
        self.areak = areak
        self.bl = bl
        self.sigma = sigma
        
        ### normalize certain variables
        
#        self.lx /= self.bl
#        self.ly /= self.bl   
#        self.dt *= self.nsamp
        
        ### define more simulation parameters
        
#        self.kT = 1
#        self.gamma_n = 1
#        self.N_avg = np.average(self.nbpc)
#        self.r_avg = self.bl*self.N_avg/2/np.pi
#        self.tau_D = self.r_avg**2 * self.gamma_n * self.N_avg / self.kT
#        self.tau_A = 2 * self.r_avg * self.gamma_n / self.fp
        
        return
        
##############################################################################

class Beads:
    """ data structure for storing particle/bead information"""
    
    def __init__(self, x, sim):
        
        ### assign bead positions
        
        self.x = x
        #self.x = x/sim.bl  
        
        ### assign cell indices to beads
        
        self.cid = np.zeros((sim.nbeads), dtype=np.int32)-1
        k = 0
        for j in range(sim.ncells):
            for n in range(sim.nbpc[j]):
                self.cid[k] = j
                k += 1
        
        return

##############################################################################

class Cells:
    """ data structure for storing cell information"""

    def __init__(self, beads, sim):
        
        ### calculate the centre of mass positions of cells 
        
        self.com = np.zeros((sim.nsteps, 2, sim.ncells), dtype=np.float32)
        for tstep in range(sim.nsteps):
            print tstep, sim.nsteps
            k = 0
            for j in range(sim.ncells):
                self.com[tstep, 0, j], self.com[tstep, 1, j] = \
                    self.compute_com_single(beads.x[tstep, 0, k:(k+sim.nbpc[j])], \
                        beads.x[tstep, 1, k:(k+sim.nbpc[j])], \
                            sim.lx, sim.ly, sim.nbpc[j])
                k += sim.nbpc[j]
        
        return  
 
    def compute_com_single(self, x, y, lx, ly, nbeads):
        """ compute the centre of mass of cells per timestep, correcting pbc on the fly"""
            
        ### correct pbcs
        
        for i in range(1, nbeads):
            dx = x[i] - x[i-1]
            dy = y[i] - y[i-1]
            x[i] = x[i-1] + neigh_min(dx,lx)
            y[i] = y[i-1] + neigh_min(dy,ly)
            
        ### compute comx and comy
        
        comx = np.average(x)
        comy = np.average(y)
        
        ### put comx and comy into simulation box
        
        comx /= lx
        comx = comx - math.floor(comx)
        comx *= lx
        comy /= ly
        comy = comy - math.floor(comy)
        comy *= ly
    
        return comx, comy
        
    def correct_for_pbc(self, sim):
        """ correct the positions for periodic boundary conditions"""
        
        self.comi = np.zeros((sim.nsteps, 2, sim.ncells), dtype=np.float32)
        self.comi[0, :, :] = self.com[0, :, :]
        for tstep in range(1, sim.nsteps):
            dx = self.com[tstep, 0, :] - self.com[tstep-1, 0, :]
            dy = self.com[tstep, 1, :] - self.com[tstep-1, 1, :]
            self.comi[tstep, 0, :] = self.comi[tstep-1, 0, :] + neigh_min_array(dx, sim.lx)
            self.comi[tstep, 1, :] = self.comi[tstep-1, 1, :] + neigh_min_array(dy, sim.ly)            

        return        
        
##############################################################################

def read_data(folder):
    """ read simulation data through hdf5 file"""
    
    ### access the file
    
    fpath = folder + '/out.h5'
    assert os.path.exists(fpath), "out.h5 does NOT exist for " + fpath
    fl = h5py.File(fpath, 'r')
    
    ### read in the positions of beads
    
    x = np.array(fl['/positions/x'], dtype=np.float32)

    ### read in the box info

    lx = fl['/info/box/x'][...]
    ly = fl['/info/box/y'][...]

    ### read in the general simulation info
    
    dt = fl['/info/dt'][...]
    nsteps = fl['/info/nsteps'][...]
    nbeads = fl['/info/nbeads'][...]
    nsamp = fl['/info/nsamp'][...]

    ### read in the cell information
    
    ncells = fl['/cell/ncells'][...]
    nbpc = np.array(fl['/cell/nbpc'], dtype=np.int32)

    ### read in the simulation parameters
    
    eps = fl['/param/eps'][...]
    rho = fl['/param/rho'][...]
    fp = fl['/param/fp'][...]
    areak = fl['/param/areak'][...]
    bl = fl['/param/bl'][...]
    sigma = fl['/param/sigma'][...]

    ### close the file

    fl.close()
    
    sim = Simulation(lx, ly, dt, nsteps, nbeads, nsamp, ncells, nbpc, eps, rho, \
               fp, areak, bl, sigma)
    beads = Beads(x, sim)
    cells = Cells(beads, sim)
    
    return beads, cells, sim
    
##############################################################################

def write_cell_data(cells, sim, folder):
    """ write cell data to hdf5 files in corresponding folders"""
    
    ### data file to write to
    
    fpath = folder + '/out_cell.h5'
    fl = h5py.File(fpath, 'w')
    
    ### positions of cells
    
    pos = fl.create_group('positions')
    pos.create_dataset('x', (sim.nsteps, 2, sim.ncells), data=cells.com, dtype=np.float64, compression='gzip') 
    pos.create_dataset('xi', (sim.nsteps, 2, sim.ncells), data=cells.comi, dtype=np.float64, compression='gzip') 
    
    
    ### simulation information
    
    info = fl.create_group('info')
    box = info.create_group('box')
    box.create_dataset('x', data=sim.lx)
    box.create_dataset('y', data=sim.ly)
    info.create_dataset('dt', data=sim.dt)
    info.create_dataset('nsteps', data=sim.nsteps)
    info.create_dataset('nbeads', data=sim.nbeads)
    info.create_dataset('nsamp', data=sim.nsamp)
    
    ### simulation parameters
    
    param = fl.create_group('param')
    param.create_dataset('eps', data=sim.eps)
    param.create_dataset('rho', data=sim.rho)
    param.create_dataset('fp', data=sim.fp)
    param.create_dataset('areak', data=sim.areak)
    param.create_dataset('bl', data=sim.bl)
    param.create_dataset('sigma', data=sim.sigma)
    
    ### cell parameters
    
    cell = fl.create_group('cell')
    cell.create_dataset('ncells', data=sim.ncells)
    cell.create_dataset('nbpc', data=sim.nbpc)
    
    fl.close()    
        
    return
        
##############################################################################

def main():

    ### get the data folder
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-fl", "--folder", help="Folder containing data")
    args = parser.parse_args()
    
    ### read the data and general informaion from the folder
    
    beads, cells, sim = read_data(args.folder)
        
    ### write the cell data to hdf5 files
    
    cells.correct_for_pbc(sim)
    write_cell_data(cells, sim, args.folder)
    
    return
    
##############################################################################

if __name__ == '__main__':
    main()    
    
##############################################################################
