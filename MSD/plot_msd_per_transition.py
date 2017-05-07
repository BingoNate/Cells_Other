
""" Plot mean square displacement per transition 
    in the phase diagram"""

### example command line arguments: 
###    -e=1.0 -f=0.5
###    -fl=/local/duman/SIMULATIONS/Cells_in_LAMMPS/density_0.8/ 
###         -sb=/usr/users/iff_th2/duman/Cells_in_LAMMPS/PLOTS/
###             -sf=MSD

##############################################################################

import sys
sys.path.append('../Utility')

import itertools
import argparse
import numpy as np
import os
import matplotlib as mpl
mpl.use('Agg', warn=False)
import matplotlib.pyplot as plt
import read_write
import misc_tools 
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

def th_msd(t, sim):
    """ theoretical expression of MSD for a single active microswimmer"""
        
    return 4*sim.Dt*t + (2*sim.v0**2/sim.Dr**2)*(t*sim.Dr + np.exp(-t*sim.Dr) - 1)
    
##############################################################################
    
def list_from_cycle(cycle):
    """ get the entire color cycle as a list without changing the state of the cycle"""
    
    ### populate a list of results
    
    first = next(cycle)
    result = [first]
    for current in cycle:
        if current == first:
            break
        result.append(current)

    ### reset iterator state
        
    for current in cycle:
        if current == result[-1]:
            break
        
    return result
    
##############################################################################
     
def plot_data(x, data, param_choice, sims, savebase, savefolder):
    """ plot the data as a function of the chosen parameter"""

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
    
    cnt = 0
    for p in data.keys():
        
        sim = sims[p]           # simulation information
        y = data[p]             # analysis data
        yth = th_msd(x, sim)    # theoretical result for MSD
    
        label = r'$\epsilon=$' + str(sim.eps) + '$,f_{m}=$' + str(sim.fp) + \
            '$,\kappa_{A}=$' + str(sim.areak)
        color = list_from_cycle(ax0._get_lines.prop_cycler)[cnt]['color']
        line0 = ax0.loglog(x/sim.tau_D, y/sim.r_avg**2, \
                         linewidth=2.0, label=label, color=color)
#        line1 = ax0.loglog(x/sim.tau_D, yth/sim.r_avg**2, '--', \
#                         linewidth=1.0, label='_nolegend_', color=color)
        cnt += 1
        
    ax0.loglog(x/sim.tau_D, x/sim.tau_D, '--', linewidth=1.0, color='grey')
    ax0.loglog(x/sim.tau_D, (x/sim.tau_D)**2, '--', linewidth=1.0, color='grey')
    
#    ax0.axvline(1./sim.Dr/sim.tau_D, color='black', alpha=0.5)
#    ax0.axhline(1.0, color='black', alpha=0.5)
    
    ### title
    
#    ax0.set_title("$t/\\tau_{D}$ = " + "{0:.2f}".format(time/sim.tau_D) + \
#        ", $t/\\tau_{A}$ = " + "{0:.2f}".format(time/sim.tau_A), fontsize=30)
    
    ### labels
        
    ax0.set_xlabel(r"$t/\tau_{D}$", fontsize=30)
    ax0.set_ylabel(r"$\Delta r^{2}/R^{2}$", fontsize=30)

    ### limits

#    ax0.set_xlim((1e-1, 5e2))
    ax0.set_ylim((1e-1, 1e6))
    ax0.set_ylim((1e-1, 1e10))
    
    ### ticks
    
#    ax0.xaxis.set_ticks([1e-1, 1e0, 1e1, 1e2, 5e2])
    ax0.yaxis.set_ticks([1e-2, 1e0, 1e2, 1e4, 1e6])
    ax0.yaxis.set_ticks([1e-2, 1e0, 1e2, 1e4, 1e6, 1e8, 1e10])
    ax0.tick_params(axis='both', which='major', labelsize=30)
    
    ### legend

    ax0.legend(bbox_to_anchor=(0.005, 0.,0.65, 1.), loc=2, borderaxespad=0., \
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

    analysistype = "MSD_subt"
    analysisdatabase = '/usr/users/iff_th2/duman/Cells_in_LAMMPS/DATA/' 
    analysisdatabase += analysistype + '/' 
    
    datafolderbase = '/local/duman/SIMULATIONS/Cells_in_LAMMPS/density_0.8/'
    
    savebase = "/usr/users/iff_th2/duman/Cells_in_LAMMPS/PLOTS/"
    savefolder = "MSD_subt"
    
    ### make the parameter choice
    
    # motility
    fp = [1.0, 5.0, 10.0]
    eps = 1.0
    areak = 10.0
    param = fp
    param_choice = 'fp'
    
    # compressibility
#    fp = 5.0
#    eps = 5.0
#    areak = [1.0, 10.0, 100.0]
#    param = areak
#    param_choice = 'areak'
    
    # adhesion
#    fp = 3.0
#    areak = 10.0
#    eps = [0.05, 1.0, 10.0]
#    param = eps
#    param_choice = 'eps'

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
