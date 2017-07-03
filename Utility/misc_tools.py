
""" A collection of helper functions and classes to do post-processing analysis
        IMPORTANT NOTE:
        All trajectory information is assumed UNWRAPPED!"""

##############################################################################

import numpy as np
import math
import read_write

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

def coords_per_cells(x, y, nbpc):
    """ get the coordinates of beads in terms of cells"""
    
    splitter = np.cumsum(nbpc)[:-1]
    x_per_cell = np.split(x, splitter)
    y_per_cell = np.split(y, splitter)
    
    return x_per_cell, y_per_cell

##############################################################################
 
def gen_labels(params, param_choice, sim):
    """ generate labels according to the parameter choice for helping in plots"""
    
    eps_sym = r'$\epsilon$'
    f_sym = r'$\f_m$'
    ka_sym = r'$\kappa_A$'
    kb_sym = r'$\kappa_B$'
    
    if param_choice == 'eps':
        xlab = eps_sym
        label = f_sym + "=" + str(sim.fp) + "," + \
            ka_sym + "=" + str(sim.areak) + "," + \
            kb_sym + "=" + str(sim.kappa)
    elif param_choice == 'fp':
        xlab = f_sym
        label = eps_sym + "=" + str(sim.eps) + "," + \
            ka_sym + "=" + str(sim.areak) + "," + \
            kb_sym + "=" + str(sim.kappa)        
    elif param_choice == 'areak':
        xlab = ka_sym
        label = eps_sym + "=" + str(sim.eps) + "," + \
            f_sym + "=" + str(sim.fp) + "," + \
            kb_sym + "=" + str(sim.kappa)        
    elif param_choice == 'kappa':
        xlab = kb_sym
        label = eps_sym + "=" + str(sim.eps) + "," + \
            f_sym + "=" + str(sim.fp) + "," + \
            ka_sym + "=" + str(sim.areak)        
    else:
        raise ValueError("Parameter choice is non-existing!")
    
    return xlab, label
 
##############################################################################
    
def gen_save_props(param_choice, sim):
    """ generate saving properties for the figure"""
    
    name = ''
    pname = ''
    if param_choice == 'areak': 
        name = 'AREAK'
        pname = name + '_eps_' + str(sim.eps) + '_fp_' + str(sim.fp) + \
            '_kappa_' + str(sim.kappa)
        xlab = '$\kappa_A$'
        lab = '$\epsilon=$' + str(sim.eps) + ',$f_m=$' + str(sim.fp) + \
            ',$\kappa_B=$' + str(sim.kappa)
    elif param_choice == 'eps':
        name = 'EPS'
        pname = name + '_fp_' + str(sim.fp) + '_areak_' + str(sim.areak) + \
            '_kappa_' + str(sim.kappa)
        xlab = '$\epsilon$'
        lab = '$f_m=$' + str(sim.fp) + ',$\kappa_A=$' + str(sim.areak) + \
            ',$\kappa_B=$' + str(sim.kappa)        
    elif param_choice == 'fp':
        name = 'FP'
        pname = name + '_eps_' + str(sim.eps) + '_areak_' + str(sim.areak) + \
            '_kappa_' + str(sim.kappa)
        xlab = '$f_{m}$'
        lab = '$\epsilon=$' + str(sim.eps) + ',$\kappa_A=$' + str(sim.areak) + \
            ',$\kappa_B=$' + str(sim.kappa)          
    elif param_choice == 'kappa':
        name = 'KAPPA'
        pname = name + '_eps_' + str(sim.eps) + '_fp_' + str(sim.fp) + \
            '_areak_' + str(sim.areak)
        xlab = '$\kappa_B$'
        lab = '$\epsilon=$' + str(sim.eps) + ',$f_m=$' + str(sim.fp) + \
            ',$\kappa_A=$' + str(sim.areak) 
    else:
        raise ValueError("Parameter choice is non-existing!")            
            
    return name, pname, xlab, lab

##############################################################################
    
def collect_data(args, analysisdatabase, read_analysis_func):
    """ collect analysis data per chosen parameter,
    where the chosen param is an array and the rest are scalars"""

    param = []
    param_choice = ''
    if args.eps == -1:
        param_choice = 'eps'
        param = [0.05, 1.0, 5.0, 20.0]
    if args.fp == -1:
        param_choice = 'fp'
        param = [0.5, 1.0, 5.0]
    if args.areak == -1:
        param_choice = 'areak'
        param = [1.0, 10.0, 100.0]
    if args.kappa == -1:
        param_choice = 'kappa'
        param = [1.0, 10.0, 100.0, 1000.0]

    data = {}       # carries the data per parameter set
    sims = {}       # carries the simulation information per parameter set

    for p in param:
        
        if param_choice == 'areak':
            datafolder, analysisfile = read_write.gen_folders(args.eps, args.fp, p, args.kappa,
                                                              args.savefolder, 
                                                              args.folder, analysisdatabase)
        elif param_choice == 'eps':
            datafolder, analysisfile = read_write.gen_folders(p, args.fp, args.areak, args.kappa,
                                                              args.savefolder, 
                                                              args.folder, analysisdatabase)
        elif param_choice == 'fp':            
            datafolder, analysisfile = read_write.gen_folders(args.eps, p, args.areak, args.kappa,
                                                              args.savefolder, 
                                                              args.folder, analysisdatabase)  
        elif param_choice == 'kappa':            
            datafolder, analysisfile = read_write.gen_folders(args.eps, args.fp, args.areak, p,
                                                              args.savefolder, 
                                                              args.folder, analysisdatabase)  
            
        sims[p] = read_write.read_sim_info(datafolder, True)
        out = read_analysis_func(analysisfile)
        if type(out) == tuple:
            if len(out) == 2:
                x, y = out
        else:
            y = out
        data[p] = y

    if type(out) == tuple:
        return x, data, param_choice, sims   
    else:
        return param, data, param_choice, sims
    
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
        
        self.kT = 1.
        self.gamma_n = 1.
        self.N_avg = np.average(self.nbpc)
        self.r_avg = self.bl*self.N_avg/2/np.pi
        self.area_avg = np.pi*self.r_avg**2*0.9
        self.Dt = self.kT/(self.gamma_n*self.N_avg)
        if self.fp == 0.:
            self.tau_A = 0.0
        else:
            self.tau_A = 2 * self.r_avg * self.gamma_n / self.fp
        self.Dr = 0.001
        self.tau_D = self.r_avg**2/4./self.Dt
        self.v0 = self.fp/self.gamma_n
        
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
        
