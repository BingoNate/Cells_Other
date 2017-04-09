
""" Calculate the self part of the intermediate scattering function of cells per parameter"""

##############################################################################

import argparse
import numpy as np
import os
import h5py
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import math
import numpy.ma as ma

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
        
        #self.lx /= self.bl
        #self.ly /= self.bl   
        self.dt *= self.nsamp
        
        ### define more simulation parameters
        
        self.kT = 1
        self.gamma_n = 1
        self.N_avg = np.average(self.nbpc)
        self.r_avg = self.bl*self.N_avg/2./np.pi
        self.tau_D = self.r_avg**2 * self.gamma_n * self.N_avg / self.kT
        self.tau_A = 2 * self.r_avg * self.gamma_n / self.fp
        
        return
        
##############################################################################

class Cells:
    """ data structure for storing cell information"""

    def __init__(self, x, sim):
        
       ### assign centre of mass cell positions
        
        self.x = x
        
        return  
        
    def correct_for_pbc(self, sim):
        """ correct the positions for periodic boundary conditions"""
        
        self.xi = np.zeros((sim.nsteps, 2, sim.ncells), dtype=np.float32)
        self.xi[0, :, :] = self.x[0, :, :]
        for tstep in range(1, sim.nsteps):
            dx = self.x[tstep, 0, :] - self.x[tstep-1, 0, :]
            dy = self.x[tstep, 1, :] - self.x[tstep-1, 1, :]
            self.xi[tstep, 0, :] = self.xi[tstep-1, 0, :] + neigh_min_array(dx, sim.lx)
            self.xi[tstep, 1, :] = self.xi[tstep-1, 1, :] + neigh_min_array(dy, sim.ly)            

        return
        
    def calculate_inter_scattering(self, sim, qvector):
        """ calculate and average the intermediate scattering function"""
        
        ### populate circular orientation of the q vector
        
        nks = 24
        kxs = np.zeros((nks), dtype=np.float64)
        kys = np.zeros((nks), dtype=np.float64)
        for j in range(nks):
            kxs[j] = np.cos(2*np.pi*j/nks) * qvector
            kys[j] = np.sin(2*np.pi*j/nks) * qvector
            
        ndelay = int(sim.nsteps/2)
        delay = np.zeros((ndelay), dtype=np.int64)
        Fs = np.zeros((ndelay), dtype=np.float64)
        Fs[0] = 1.0
        
        ### calculate the intermediate scattering function
        
        for d in range(1, ndelay):
            
            delay[d] = d*sim.dt
            dx = self.xi[d:, 0, :] - self.xi[:-d, 0, :]
            dy = self.xi[d:, 1, :] - self.xi[:-d, 1, :]
            
            ### take an average over the circular orientation of the q vector
            
            avg_over_qvector = 0.
            for j in range(nks):
                
                ### sum cosine and sum terms over the ensemble of cells
                
                cost = np.sum(np.cos(kxs[j]*dx), axis=1)
                sint = np.sum(np.sin(kys[j]*dy), axis=1)
                
                ### take a time average of the summed up terms and add it up to the average
                
                avg_over_qvector += np.mean(cost**2 + sint**2)/sim.ncells

            Fs[d] /= avg_over_qvector/nks

        return delay, Fs        
        
##############################################################################

def read_data(folder):
    """ read simulation data through hdf5 file"""
    
    ### access the file
    
    fpath = folder + '/out_cell.h5'
    assert os.path.exists(fpath), "out_cell.h5 does NOT exist for " + fpath
    fl = h5py.File(fpath, 'r')
    
    ### read in the positions of centre of mass of cells
    
    x = np.asarray(fl['/positions/x'], dtype=np.float64)
    #xi = np.asarray(fl['/positions/xi'], dtype=np.float64)

    ### read in the box info

    lx = fl['/info/box/x'][...]
    ly = fl['/info/box/y'][...]

    ### read in the general simulation info
    
    dt = fl['/info/dt'][...]
    nsteps = fl['/info/nsteps'][...]
    nbeads = fl['/info/nbeads'][...]
    nsamp = fl['/info/nsamp'][...]

    ### read in the cell information
    
    #ncells = fl['/cell/ncells'][...]
    #nbpc = np.asarray(fl['/cell/nbpc'], dtype=np.int32)
    ncells = 5000
    nbpc = np.zeros((ncells), dtype=np.int32)

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
    cells = Cells(x, sim)
    
    return cells, sim

##############################################################################
        
class Subplots:
    """ plot structure"""
    
    totcnt = -1             # Total number of subplots 
    
    def __init__(self, f, l, s, b, t):
        self.fig = f        # Figure axes handle
        self.length = l     # Length of the subplot box 
        self.sep = s        # Separation distance between subplots 
        self.beg = b        # Beginning (offset) in the figure box
        self.tot = t        # Total number of subplots in the x direction
        
        return
        
    def addSubplot(self):
        """ add a subplot in the grid structure"""
        
        ### increase the number of subplots in the figure
        
        self.totcnt += 1
        
        ### get indices of the subplot in the figure
        
        self.nx = self.totcnt%(self.tot)
        self.ny = self.totcnt/(self.tot)
        
        self.xbeg = self.beg + self.nx*self.length + self.nx*self.sep
        self.ybeg = self.beg + self.ny*self.length + self.ny*self.sep
        
        return self.fig.add_axes([self.xbeg,self.ybeg,self.length,self.length])

##############################################################################
        
def write_analysis_data(x, y, sim):
    """ write analysis data to the corresponding file"""
    
    base = "/usr/users/iff_th2/duman/Cells_in_LAMMPS/DATA/"
    base += "Inter_scattering/"
    os.system("mkdir -p " + base)
    fpath = base + "inter_scattering_eps_" + \
        str(sim.eps) + "_fp_" + str(sim.fp) + "_areak_" + str(sim.areak) + ".txt"
    
    fl = open(fpath, 'w')
    
    N = len(x)
    for j in range(N):
        fl.write(str(x[j]) + '\t\t' + str(y[j]) + '\n')
    
    fl.close()

    return
        
##############################################################################

def main():

    ### get the data folder
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-fl", "--folder", help="Folder containing data")     
    parser.add_argument("-s","--save_eps", action="store_true", help="Decide whether to save in eps or not")            
    args = parser.parse_args()
    
    ### read the data and general information from the folder
    
    cells, sim = read_data(args.folder)
        
    ### calculate the intermediate scattering fnc of the centre of mass of cells

    cells.correct_for_pbc(sim)
    qvector = 1./4.
    delay, Fs = cells.calculate_inter_scattering(sim, qvector)
    
    ### write the overlap fnc data to the corresponding file
    
    write_analysis_data(delay, Fs, sim)
    
    return
    
##############################################################################

if __name__ == '__main__':
    main()    
    
##############################################################################
