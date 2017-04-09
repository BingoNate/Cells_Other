
""" Plot mean square displacement per area constraint potential parameter"""

### example command line arguments: 
###    -e=1.0 -f=0.5
###    -fl=/local/duman/SIMULATIONS/Cells_in_LAMMPS/density_0.8/ 
###         -sb=/usr/users/iff_th2/duman/Cells_in_LAMMPS/PLOTS/
###             -sf=MSD

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
   
##############################################################################
     
def plot_data(x, data, sims, savebase, savefolder, save_eps):
    """ plot the data as a function of the chosen parameter"""

    ### set general plot properties

    base = savebase + savefolder + '/'
    os.system("mkdir -p " + base)
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
    
    for a in data.keys():
        
        sim = sims[a]       # simulation information
        y = data[a]         # analysis data
    
        label = '$\\epsilon=$' + str(sim.eps) + '$,f=$' + str(sim.fp) + '$,\\kappa_{A}=$' + str(sim.areak)
        line0 = ax0.loglog(x/sim.tau_D, y/sim.r_avg**2, \
                         linewidth=2.0, label=label)
    
    ### title
    
#    ax0.set_title("$t/\\tau_{D}$ = " + "{0:.2f}".format(time/sim.tau_D) + \
#        ", $t/\\tau_{A}$ = " + "{0:.2f}".format(time/sim.tau_A), fontsize=30)
    
    ### labels
        
    ax0.set_xlabel("$t/\\tau_{D}$", fontsize=40)
    ax0.set_ylabel("$\\Delta r^{2}/R^{2}$", fontsize=40)

    ### limits

    #ax0.set_xlim((-1, 15))
    #ax0.set_ylim((downlim, uplim))
    
    ### ticks
    
    #ax0.xaxis.set_ticks(np.linspace(0, 15, num_ticks, endpoint=True))
    #ax0.yaxis.set_ticks(np.linspace(0, uplim, num_ticks, endpoint=True))
    ax0.tick_params(axis='both', which='major', labelsize=30)
    
    ### legend

    ax0.legend(bbox_to_anchor=(1.005, 0.,0.65, 1.), loc=2, borderaxespad=0., \
        prop={'size': 20}, mode="expand", frameon=False)
    
    ### save

    savepath1 = base + savefolder + "_per_AREAK_eps_" + str(sim.eps) + "_fp_" + str(sim.fp) + ".png"
    if save_eps:
        savepath2 = base + savefolder + "_per_AREAK_eps_" + str(sim.eps) + "_fp_" + str(sim.fp) + ".eps"            
    plt.savefig(savepath1, dpi=200, bbox_inches='tight', pad_inches=0.08)
    if save_eps:
        plt.savefig(savepath2, dpi=200, bbox_inches='tight', pad_inches=0.08)        
    fig.clf()                       
        
    return
        
##############################################################################

def main():

    ### get the data folder
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--eps", type=float, help="Strength of LJ potential")
    parser.add_argument("-f", "--fp", type=float, help="Propulsion force")
    parser.add_argument("-fl", "--folder", nargs="?", \
                        const='/local/duman/SIMULATIONS/Cells_in_LAMMPS/density_0.8/', \
                        help="Folder containing data, as in /local/duman/SIMULATIONS/Cells_in_LAMMPS/density_0.8/")    
    parser.add_argument("-sb", "--savebase", nargs="?", \
                        const = "/usr/users/iff_th2/duman/Cells_in_LAMMPS/PLOTS/", \
                        help="Folder to save the data, as in /usr/users/iff_th2/duman/Cells_in_LAMMPS/PLOTS/") 
    parser.add_argument("-sf", "--savefolder", nargs="?", const="MSD", \
                        help="Specific folder for saving, as in MSD")   
    parser.add_argument("-s","--save_eps", action="store_true", help="Decide whether to save in eps or not")            
    args = parser.parse_args()
    
    ### accummulate data as a function of areak -with a dictionary-

    analysisdatabase = '/usr/users/iff_th2/duman/Cells_in_LAMMPS/DATA/'  
    analysisdatabase += args.savefolder + '/'    
    areak = [1.0, 10.0, 100.0]
    data = {}       # carries the data per parameter set
    sims = {}       # carries the simulation information per parameter set

    for a in areak:
        
        datafolder, analysisfile = read_write.gen_folders(args.eps, args.fp, a, args.savefolder, 
                                               args.folder, analysisdatabase)
        sims[a] = read_write.read_sim_info(datafolder)
        x, y = read_write.read_2d_analysis_data(analysisfile)
        data[a] = y
        
    ### plot the data as a function of the parameter
    
    plot_data(x, data, sims, args.savebase, args.savefolder, args.save_eps)

    return
    
##############################################################################

if __name__ == '__main__':
    main()    
    
##############################################################################
