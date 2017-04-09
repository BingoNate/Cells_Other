
""" Calculate the MSD of the centre of mass of cells per parameter"""

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

    def __init__(self, lx, ly, dt, nsteps, nbeads, nsamp, nfils, nbpf, density, kappa, km, pa, bl, sigma):
        
        self.lx = lx
        self.ly = ly
        self.dt = dt
        self.nsteps = nsteps
        self.nbeads = nbeads
        self.nsamp = nsamp
        self.nfils = nfils
        self.nbpf = nbpf
        self.density = density
        self.kappa = kappa
        self.km = km
        self.pa = pa
        self.bl = bl
        self.sigma = sigma
        self.density = "{:.1f}".format(float(density))
        self.kappa = "{:.1f}".format(float(kappa))
        self.km = "{:.1f}".format(float(km))
        self.pa = "{:.1f}".format(float(pa))
        
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

class Cells:
    """ data structure for storing cell information"""

    def __init__(self, x, xi, sim):
        
       ### assign centre of mass cell positions
        
        self.x = x 
        self.xi = xi
        
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
        
    def calculate_MSD(self, sim):
        """ calculate and average the MSD (time and ensemble avgs)"""
            
        ndelay = int(sim.nsteps/2)
        delay = np.zeros((ndelay), dtype=np.int64)
        msd = np.zeros((ndelay), dtype=np.float64)
        
        for d in range(1, ndelay):
            delay[d] = d*sim.dt
            msd[d] = np.mean(np.mean( \
                (self.xi[d:,0,:]-self.xi[:-d,0,:])**2 + \
                    (self.xi[d:,1,:]-self.xi[:-d,1,:])**2, axis=1) )
            
        return delay, msd        
        
##############################################################################

def read_data(folder):
    """ read simulation data through hdf5 file"""
  
    ### access the file
    
    fpath = folder + '/out_fil.h5'
    assert os.path.exists(fpath), "out_fil.h5 does NOT exist for " + fpath
    fl = h5py.File(fpath, 'r')
    
    ### read in the positions of filaments
    
    x = np.array(fl['/positions/x'], dtype=np.float32)
    xi = np.array(fl['/positions/xi'], dtype=np.float32)

    ### read in the box info

    lx = fl['/info/box/x'][...]
    ly = fl['/info/box/y'][...]

    ### read in the general simulation info
    
    dt = fl['/info/dt'][...]
    nsteps = fl['/info/nsteps'][...]
    nbeads = fl['/info/nbeads'][...]
    nsamp = fl['/info/nsamp'][...]

    ### read in the filament information
    
    nfils = fl['/info/nfils'][...]
    nbpf = fl['/info/nbpf'][...]

    ### read in the simulation parameters
    
    density = fl['/param/density'][...]
    kappa = fl['/param/kappa'][...]
    km = fl['/param/km'][...]
    pa = fl['/param/pa'][...]
    bl = fl['/param/bl'][...]
    sigma = fl['/param/sigma'][...]

    ### close the file

    fl.close()
    
    sim = Simulation(lx, ly, dt, nsteps, nbeads, nsamp, nfils, nbpf, density, kappa, km, pa, bl, sigma)
    fils = Cells(x, xi, sim)
    
    return fils, sim

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
    
    base = "/usr/users/iff_th2/ravich/Cells_in_LAMMPS/DATA/"
    base += "MSD/"
    os.system("mkdir -p " + base)
    fpath = base + "msd_eps_" + \
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
    print cells.xi
        
    ### calculate MSD of the centre of mass of cells
    
    #cells.correct_for_pbc(sim)
    delay, msd = cells.calculate_MSD(sim)
    
    ### write the MSD data to the corresponding file
    
    write_analysis_data(delay, msd, sim)
    
    return
    
##############################################################################

if __name__ == '__main__':
    main()    
    
##############################################################################
