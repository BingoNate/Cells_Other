
""" Plot energy spectrum per transition 
    in the phase diagram"""

### example command line arguments: 
###    -e=1.0 -f=0.5
###    -fl=/local/duman/SIMULATIONS/Cells_in_LAMMPS/density_0.8/ 
###         -sb=/usr/users/iff_th2/duman/Cells_in_LAMMPS/PLOTS/
###             -sf=Energy_spectrum

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
import seaborn as sns
from scipy.optimize import curve_fit
from scipy.stats.kde import gaussian_kde
#sns.set()
sns.set(style="white",context='paper',
        font_scale=1.2,font="Open Sans",
        rc={'mathtext.default': 'regular','font.size': 30, 
            'font.family': 'sans',"figure.dpi":300,
            "xtick.major.size": 8, "ytick.major.size": 8,
            'grid.linestyle': '--'})   


##########################################################################

def up_power_law(x, a, b):
    return a * x**b

##########################################################################

def down_power_law(x, a, b):
    return a * x**(-b)

##########################################################################

def exp_law(x, a, b):
    return a * np.exp(-x/b)

##########################################################################

def power_law_with_exp_tail(x, a, b):
    return (x**(-a))*np.exp(-x/b)
    
##########################################################################

def power_law_with_up_exp_tail(x, a, b):
    return (x**(-a))*np.exp(x/b)    

##########################################################################

def mixed_power_law_with_exp_tail(x, a, b, c, d, f):
    return (x**(-a))*np.exp(-x/b) + c*(x**(d))*np.exp(-x/f)

##########################################################################
    
def mixed_power_law_with_exp_tail_2(x, a, b, c, d):
    return (x**(b))*np.exp(-x/c) + a*(x**(-d))  

##########################################################################

def stretch_exp_law(x, a, b, c):
    return a * np.exp(-(x/c)**b)
    
##############################################################################
     
def plot_data(x, data, param_choice, sims, savebase, savefolder):
    """ plot the data as a function of the chosen parameter"""

    ### set normalization parameter 
    
    knorm = 2.*np.pi/sims[data.keys()[0]].r_avg
                          
    ### set general plot properties

    #downlim = -1
    #uplim = sim.lx/4.
    num_ticks = 5
    ax_len = 1.0                          # Length of one subplot square box
    ax_b = 0.0                            # Beginning/offset of the subplot in the box
    ax_sep = 0.0                          # Separation length between two subplots
    total_subplots_in_x = 1               # Total number of subplots    
    fig = plt.figure()
    subp = misc_tools.Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x) 
    ax0 = subp.addSubplot()
                          
    ### plot 

    subp = misc_tools.Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x) 
    ax0 = subp.addSubplot()
    
    for p in data.keys():
        
        sim = sims[p]       # simulation information
        y = data[p]         # analysis data
        
        popt, pcov = curve_fit(up_power_law, x, y)
        yfit = up_power_law(x, popt[0], popt[1])
        print "key = ", p,", fit param = ", popt[1]    

        label = r'$\epsilon=$' + str(sim.eps) + '$,f_{m}=$' + str(sim.fp) + \
            '$,\kappa_{A}=$' + str(sim.areak)
        line0 = ax0.loglog(x/knorm, y, \
                         linewidth=2.0, label=label)
        line1 = ax0.loglog(x/knorm, yfit, '--', \
                         linewidth=1.0, label='_nolegend_')    
    ### labels
        
    ax0.set_xlabel(r"$q/q_{R}$", fontsize=30)
    ax0.set_ylabel(r"$E(q)$", fontsize=30)

    ### limits

#    ax0.set_xlim((0, 22))
#    ax0.set_ylim((0.0, 1.9))
    
    ### ticks
    
#    ax0.xaxis.set_ticks([0, 10, 20])
#    ax0.yaxis.set_ticks([0, 0.5, 1.0, 1.5])
    ax0.tick_params(axis='both', which='major', labelsize=30)
    
    ### legend

    ax0.legend(bbox_to_anchor=(0.35, 0.,0.65, 1.), loc=2, borderaxespad=0., \
        prop={'size': 20}, mode="expand", frameon=False)
    
    ### save

    base = savebase + savefolder + '/'
    os.system("mkdir -p " + base)  
    
    savepath1 = base + savefolder + "_transition_per_" + param_choice + ".png"
    savepath2 = base + savefolder + "_transition_per_" + param_choice + ".eps"
    plt.savefig(savepath1, dpi=300, bbox_inches='tight', pad_inches=0.08)
    plt.savefig(savepath2, dpi=300, bbox_inches='tight', pad_inches=0.08)        
    fig.clf()                          
        
    return
        
##############################################################################

def main():

    ### get the folder structure

    analysistype = "Energy_spectrum"
    analysisdatabase = '/usr/users/iff_th2/duman/Cells_in_LAMMPS/DATA/' 
    analysisdatabase += analysistype + '/' 
    
    datafolderbase = '/local/duman/SIMULATIONS/Cells_in_LAMMPS/density_0.8/'
    
    savebase = "/usr/users/iff_th2/duman/Cells_in_LAMMPS/PLOTS/"
    savefolder = "Energy_spectrum"
    
    ### make the parameter choice
    
    # motility
#    fp = [1.0, 3.0, 10.0]
#    eps = 1.0
#    areak = 10.0
#    param = fp
#    param_choice = 'fp'
    
    # compressibility
#    fp = 5.0
#    eps = 5.0
#    areak = [1.0, 10.0, 100.0]
#    param = areak
#    param_choice = 'areak'
    
    # adhesion
    fp = 3.0
    areak = 10.0
    eps = [0.05, 1.0, 10.0]
    param = eps
    param_choice = 'eps'

    data = {}       # carries the data per parameter set
    sims = {}       # carries the simulation information per parameter set

    for p in param:
        
        if param_choice == 'areak':
            datafolder, analysisfile = read_write.gen_folders(eps, fp, p, savefolder, 
                                                   datafolderbase, analysisdatabase)
        elif param_choice == 'eps':
            datafolder, analysisfile = read_write.gen_folders(p, fp, areak, savefolder, 
                                                   datafolderbase, analysisdatabase)
        elif param_choice == 'fp':            
            datafolder, analysisfile = read_write.gen_folders(eps, p, areak, savefolder, 
                                                   datafolderbase, analysisdatabase)  
            
        sims[p] = read_write.read_sim_info(datafolder)
        x, y = read_write.read_2d_analysis_data(analysisfile)
        data[p] = y
        
    ### plot the data as a function of the parameter
    
    plot_data(x, data, param_choice, sims, savebase, savefolder)

    return
    
##############################################################################

if __name__ == '__main__':
    main()    
    
##############################################################################
