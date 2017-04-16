
""" A collection of helper functions and classes to do post-processing analysis
        IMPORTANT NOTE:
        All trajectory information is assumed UNWRAPPED!"""

##############################################################################

import numpy as np
import math

############################################################################

def nearbyint(x):
    """ Round to the nearby integer"""
    
    if x >= 0:
        return math.floor(x+0.5)
    else:
        return math.floor(x-0.5)

############################################################################
    
def min_img_dist(x1, x2, lx):
    """ compute the minimum image distance between two positions"""
    
    dx = x2 - x1 
    return dx-nearbyint(dx/lx)*lx

############################################################################

def get_img(x, lx):
    """ get the image position in the central box 
    -- can be numpy array or single pos"""
    
    return x-np.floor(x/lx)*lx

##############################################################################

class Beads:
    """ data structure for storing particle/bead information"""
    
    def __init__(self, x, cid):
        
        ### assign bead positions
        
        self.xu = x
        
        ### assign cell indices to beads
        
        self.cid = cid
        
        return
        
    ##############
    
    def get_img_pos(self, lx):
        """ get the image positions of beads in the central box"""
        
        self.xi = get_img(self.xu, lx)
        
        return
        
##############################################################################

class Cells:
    """ data structure for storing cell information"""
    
    def __init__(self, x, d, nbpc, sim):
        
        ### assign bead positions
        
        self.xu = x
        
        ### assign polarity
        
        self.pol = d
        
        ### assign number of beads info
        
        self.nbpc = nbpc
                
        return

    ##############
    
    def get_img_pos(self, lx):
        """ get the image positions of cells in the central box"""
        
        self.xi = get_img(self.xu, lx)
        
        return
        
##############################################################################

class Simulation:
    """ data structure for storing general simulation information"""

    def __init__(self, lx, ly, dt, nsteps, ncells, nbeads, nsamp, nbpc, \
                     eps, rho, fp, areak, bl, sigma, kappa=100.0):
        
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
        self.kappa = kappa
        
        ### normalize certain variables
          
        self.dt *= self.nsamp
        
        ### define more simulation parameters
        
        self.kT = 1
        self.gamma_n = 1
        self.N_avg = np.average(self.nbpc)
        self.r_avg = self.bl*self.N_avg/2/np.pi
        self.tau_D = self.r_avg**2 * self.gamma_n * self.N_avg / 4. / self.kT
        if self.fp == 0.:
            self.tau_A = 0.0
        else:
            self.tau_A = 2 * self.r_avg * self.gamma_n / self.fp
        
        return
        
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
        
    def addSubplot(self, multi=False):
        """ add a subplot in the grid structure"""
        
        ### increase the number of subplots in the figure
        
        self.totcnt += 1
        
        ### get indices of the subplot in the figure
        
        self.nx = self.totcnt % self.tot
        self.ny = self.totcnt / self.tot
        
        self.xbeg = self.beg + self.nx*self.length + self.nx*self.sep
        self.ybeg = self.beg + self.ny*self.length + self.ny*self.sep
            
        if (multi):
            if (self.nx > 0):
                self.xbeg -= 0.6
                
        return self.fig.add_axes([self.xbeg,self.ybeg,self.length,self.length])

##############################################################################
        
