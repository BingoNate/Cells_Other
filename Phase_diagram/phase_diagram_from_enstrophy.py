
""" Build a phase diagram based on the enstrophy"""

### example command line arguments: 
###    -fl=/local/duman/SIMULATIONS/Cells_in_LAMMPS/.../ 
###         -sb=/usr/users/iff_th2/duman/Cells_in_LAMMPS/MOVIES/
###             -ti=10 -tf=2000

##############################################################################

import sys
sys.path.append('../Utility')

import argparse
import numpy as np
import os
import matplotlib as mpl
mpl.use('Agg', warn=False)
import matplotlib.pyplot as plt
import read_write
import misc_tools
import glob
import pandas as pd
from string import atof
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from matplotlib._png import read_png
import colorsys 
import seaborn as sns
from scipy.stats.kde import gaussian_kde
#sns.set()
sns.set(style="white",context='paper',
        font_scale=1.2,font="Open Sans",
        rc={'mathtext.default': 'regular','font.size': 30, 
            'font.family': 'sans',"figure.dpi":300,
            "xtick.major.size": 8, "ytick.major.size": 8,
            'grid.linestyle': '--'})

##############################################################################

class Phase:
    """ properties of the point in phase space"""        
        
    def __init__(self, e, f, a, x, threshold):
        
        self.eps = e                  # adhesion
        self.fm = f                   # motility 
        self.areak = a                # deformability
        self.data = x                 # enstrophy
        
        ### determine the type of phase point based on the enstrophy
        
        if x < threshold:
            self.type = "solidlike"
        elif x >= threshold:
            self.type = "fluidlike"

        if a == 1.0:
            if e == 0.05:
                if f >= 1.0:
                    self.type = "none"
            elif e == 0.5 or e == 1.0:
                if f >= 3.0:
                    self.type = "none"
            else:
                if f == 10.0:
                    self.type = "none"
                    
        self.set_plot_props()
        
        return

    #############            
        
    def set_plot_props(self):
        """ determine plot properties based on the type"""
        
        if self.type == "fluidlike":
            self.marker = "s"
            self.color = "maroon"
            
        elif self.type == "turbulent":
            self.marker = "D"
            self.color = "red"
            
        elif self.type == "solidlike":
            self.marker = "o"
            self.color = "cyan"
            
        elif self.type == "none":
            self.marker = "x"
            self.color = "black"
            
        
        return

##############################################################################

def plot_slices(phases, slicer):
    """ plot phase diagram by slices in 2D"""

    npoints = len(phases)       # total number of data points
    
    ### load the data
    
    idx = []
    x = []
    y = []
    z = []
    markers = []
    colors = []
    for j in xrange(npoints):
        idx.append(j)
        x.append(phases[j].eps)
        y.append(phases[j].fm)
        z.append(phases[j].data)
        markers.append(phases[j].marker)
        colors.append(phases[j].color)
    
    ### convert the lists into numpy arrays
    
    x = np.array(x, dtype=np.float32)
    y = np.array(y, dtype=np.float32)
    z = np.array(z, dtype=np.float32)
    idx = np.array(idx, dtype=np.int32)
    
    ### set plot properties
 
    savefolder = "/usr/users/iff_th2/duman/Cells_in_LAMMPS/PLOTS/PHASE_DIAGRAM/"
    os.system("mkdir -p " + savefolder)
    sfile = savefolder + "phase_diagram_from_enstrophy_per_AREAK_" + \
        str(slicer)
    fig = plt.figure()
    ax_len = 0.9                          # Length of one subplot square box
    ax_b = 0.0                            # Beginning/offset of the subplot in the box
    ax_sep = 0.0                          # Separation length between two subplots
    total_subplots_in_x = 1               # Total number of subplots in the x direction  
    tick_num = 5                          # Number of ticks
    subp = misc_tools.Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x)  
    downlim_x = 0.03
    uplim_x = 35.0
    downlim_y = 0.3
    uplim_y = 15.0
    quant_steps = 2056
    norm = mpl.colors.Normalize(vmin=1e-7, vmax=1e-3)  

    ### make these smaller to increase the resolution
         
    ax = subp.addSubplot()
    
    k = 0
    for xi, yi, zi, mi, ci in zip(x, y, z, markers, colors):
        if phases[k].type == "None":
            xi = 0.0
            yi = 0.0        
        line0 = ax.scatter(xi, yi, marker=mi, c=zi, s=30,
                   cmap=plt.cm.get_cmap('jet',quant_steps),
                    edgecolors='None', alpha=1.0, norm=norm, 
                    vmin=1e-7, vmax=1e-3)         
        k += 1
        
    #ax.grid("on")
    
    cax = plt.colorbar(line0, ax=ax)
    cax.ax.set_title(r'$\Omega$', fontsize=30)
    
    ax.set_title(r'$\kappa_{A} = $' + str(slicer), fontsize=30)
        
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel(r'$\epsilon$', fontsize=30)        
    ax.set_ylabel(r'$f_{m}$', fontsize=30)
    ax.set_xlim((downlim_x,uplim_x))
    ax.set_ylim((downlim_y,uplim_y))
    #ax.xaxis.set_ticks( np.linspace(downlim, uplim, num=tick_num, endpoint=True) )
    #ax.yaxis.set_ticks( np.linspace(downlim, uplim, num=tick_num, endpoint=True) )
    ax.tick_params(axis='both', which='major', labelsize=30)   
    
    ### save the figure
    
    plt.savefig(sfile + '.eps', dpi=300, bbox_inches='tight', pad_inches=0.08)
    plt.savefig(sfile + '.png', dpi=300, bbox_inches='tight', pad_inches=0.08)    
    plt.clf()    
    
    return
                
##############################################################################
                
def load_data(filepath, totalData):
    """ load the data into a dictionary"""

    fl = open(filepath, 'r')
    
    data = {}
    for j in xrange(totalData):
        
        line = fl.readline()
        e, f, a, x = line.split()
        key = (float(e), float(f), float(a))
        data[key] = float(x)
        
    fl.close()
        
    return data

##############################################################################
    
def construct_phase_diagram(data, idx, threshold):
    """ build the phase diagram slicing by idx"""

    phases = []
    for e, f, a in data.keys():
        if a == idx:
            phase = Phase(e, f, a, data[(e,f,a)], threshold)
            phases.append(phase)
                
    return phases
              
##############################################################################

def main():

    ### parse the command line arguments
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--sliceby", type=float, \
                        help="Slice by argument in areak")
    parser.add_argument("-t", "--threshold", type=float, \
                        help="Threshold value of correlation length to distinguish the phases")
    args = parser.parse_args()
    
    ### index the phase space

    eps = [0.05, 0.5, 1.0, 5.0, 10.0, 20.0]
    fp = [0.0, 0.5, 1.0, 3.0, 5.0, 10.0]
    areak = [1.0, 10.0, 100.0]
    totalData = len(eps)*len(fp)*len(areak)
    
    ### load the data into a dictionary

    filepath = \
        '/usr/users/iff_th2/duman/Cells_in_LAMMPS/DATA/PHASE_DIAGRAM/PHASE_DIAGRAM_from_enstrophy.txt'    
    data = load_data(filepath, totalData)

    ### build the phase diagram

    phases = construct_phase_diagram(data, args.sliceby, args.threshold)
    
    ### plot the phase diagram
    
    plot_slices(phases, args.sliceby)
    
    
##############################################################################
    
if __name__ == '__main__':
    main()
    
##############################################################################
    
