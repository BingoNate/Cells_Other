
""" Plot the centre of mass trajectories of cells per simulation"""

### example command line arguments: 
###    -fl=/local/duman/SIMULATIONS/Cells_in_LAMMPS/.../ 
###         -sb=/usr/users/iff_th2/duman/Cells_in_LAMMPS/PLOTS/
###             -sf=TRAJECTORY
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

##############################################################################

def plot_traj(x, sim, ti, tf, savebase, savefolder, save_eps):
    """ plot trajectory within the specified time window"""

    ### normalize variables for plotting purposes
    
    x /= sim.bl
    lx = sim.lx/sim.bl
    ly = sim.ly/sim.bl
    
    ### set general plot properties
    
    base = savebase + savefolder + '/'
    os.system("mkdir -p " + base)
    downlim = -2
    uplim = lx+2
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
    
    label = '$\\epsilon=$' + str(sim.eps) + '$,\\f=$' + str(sim.fp) + '$,\\kappa$' + str(sim.areak)
    line0 = ax0.plot(x[ti:tf, 0, :], x[ti:tf, 1, :], \
                     linewidth=2.0, label=label)
    
    ax0.axis('scaled')
    

    ### title
    
#    ax0.set_title("$t/\\tau_{D}$ = " + "{0:.2f}".format(time/sim.tau_D) + \
#        ", $t/\\tau_{A}$ = " + "{0:.2f}".format(time/sim.tau_A), fontsize=30)
    
    ### labels
        
    ax0.set_xlabel("$x/r_{0}$", fontsize=40)
    ax0.set_ylabel("$y/r_{0}$", fontsize=40)

    ### limits

    #ax0.set_xlim((downlim, uplim))
    #ax0.set_ylim((downlim, uplim))
    
    ### ticks
    
    #ax0.xaxis.set_ticks(np.linspace(0, uplim, num_ticks, endpoint=True))
    #ax0.yaxis.set_ticks(np.linspace(0, uplim, num_ticks, endpoint=True))
    ax0.tick_params(axis='both', which='major', labelsize=30)
    
    ### save

    savepath1 = base + savefolder + "_eps_" + str(sim.eps) + "_fp_" + str(sim.fp) + \
        "_areak_" + str(sim.areak) + ".png"
    if save_eps:
        savepath2 = base + savefolder + "_eps_" + str(sim.eps) + "_fp_" + str(sim.fp) + \
            "_areak_" + str(sim.areak) + ".eps"            
    plt.savefig(savepath1, dpi=200, bbox_inches='tight', pad_inches=0.08)
    if save_eps:
        plt.savefig(savepath2, dpi=200, bbox_inches='tight', pad_inches=0.08)        
    fig.clf()                

        
    return
        
##############################################################################

def main():
    

    ### get the data folder
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-fl", "--folder", \
                        help="Folder containing data, as in /local/duman/SIMULATIONS/Cells_in_LAMMPS/.../")
    parser.add_argument("-sb", "--savebase", nargs="?", \
                        const = "/usr/users/iff_th2/duman/Cells_in_LAMMPS/PLOTS/", \
                        help="Folder to save the data, as in /usr/users/iff_th2/duman/Cells_in_LAMMPS/PLOTS/") 
    parser.add_argument("-sf", "--savefolder", nargs="?", const="TRAJECTORY", \
                        help="Specific folder for saving, as in TRAJECTORY")   
    parser.add_argument("-ti","--init_time", nargs="?", const=10, type=int, \
                        help="First frame of the trajectory, you can also leave it empty")
    parser.add_argument("-tf","--fin_time", nargs="?", const=2000, type=int, \
                        help="Last frame of the trajectory, you can also leave it empty")    
    parser.add_argument("-s","--save_eps", action="store_true", help="Decide whether to save in eps or not")            
    args = parser.parse_args()
    
    ### read the data and general information from the folder
    
    sim, cells, beads = read_write.read_h5_file(args.folder)
    print "folder = ", args.folder, "cells.xu = \n", cells.xu
        
    ### plot the data in the given time window
    
    plot_traj(cells.xu, sim, args.init_time, args.fin_time, \
              args.savebase, args.savefolder, args.save_eps)
    
    return
    
##############################################################################

if __name__ == '__main__':
    main()    
    
##############################################################################
