
""" Plot energy spectrum per simulation"""

### example command line arguments: 
###    -e=1.0 -f=0.5 -k=10.0
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
from scipy.optimize import curve_fit

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
     
def plot_data(x, y, sim, savebase, savefolder, save_eps):
    """ plot the data as a function of the chosen parameter"""

    ### set normalization parameter 
    
    knorm = 2.*np.pi/sim.r_avg

    ### curve fitting
    
    popt, pcov = curve_fit(up_power_law, x, y)
    yfit = up_power_law(x, popt[0], popt[1])
    print popt[1]
                          
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
    
    label = '$\\epsilon=$' + str(sim.eps) + '$,f=$' + str(sim.fp) + '$,\\kappa_{A}=$' + str(sim.areak)
    line0 = ax0.plot(x/knorm, y, \
                     linewidth=2.0, label=label)
    line0 = ax0.loglog(x/knorm, yfit, '--', \
                     linewidth=2.0, label=label)    
    ### title
    
#    ax0.set_title("$t/\\tau_{D}$ = " + "{0:.2f}".format(time/sim.tau_D) + \
#        ", $t/\\tau_{A}$ = " + "{0:.2f}".format(time/sim.tau_A), fontsize=30)
    
    ### labels
        
    ax0.set_xlabel("$q/q_{R}$", fontsize=40)
    ax0.set_ylabel("$E_{q}$", fontsize=40)

    ### limits

    #ax0.set_xlim((-1, 15))
    #ax0.set_ylim((downlim, uplim))
    
    ### ticks
    
    #ax0.xaxis.set_ticks(np.linspace(0, 15, num_ticks, endpoint=True))
    #ax0.yaxis.set_ticks(np.linspace(0, uplim, num_ticks, endpoint=True))
    ax0.tick_params(axis='both', which='major', labelsize=30)
    
    ### legend

#    ax0.legend(bbox_to_anchor=(1.005, 0.,0.65, 1.), loc=2, borderaxespad=0., \
#        prop={'size': 20}, mode="expand", frameon=False)
    
    ### save

    savepath1 = base + savefolder + "_eps_" + str(sim.eps) + "_fp_" + str(sim.fp) + \
        '_areak_' + str(sim.areak) + ".png"
    if save_eps:
        savepath2 = base + savefolder + "_eps_" + str(sim.eps) + "_fp_" + str(sim.fp) + \
            '_areak_' + str(sim.areak) + ".eps"            
    plt.savefig(savepath1, dpi=300, bbox_inches='tight', pad_inches=0.08)
    if save_eps:
        plt.savefig(savepath2, dpi=300, bbox_inches='tight', pad_inches=0.08)        
    fig.clf()                       
        
    return
        
##############################################################################

def main():
    
    ### get the data folder
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--eps", type=float, help="Strength of LJ potential")
    parser.add_argument("-f", "--fp", type=float, help="Propulsion force")
    parser.add_argument("-k", "--areak", type=float, help="Strength of area constraint potential")
    parser.add_argument("-fl", "--folder", nargs="?", \
                        const='/local/duman/SIMULATIONS/Cells_in_LAMMPS/density_0.8/', \
                        help="Folder containing data, as in /local/duman/SIMULATIONS/Cells_in_LAMMPS/density_0.8/")    
    parser.add_argument("-sb", "--savebase", nargs="?", \
                        const = "/usr/users/iff_th2/duman/Cells_in_LAMMPS/PLOTS/", \
                        help="Folder to save the data, as in /usr/users/iff_th2/duman/Cells_in_LAMMPS/PLOTS/") 
    parser.add_argument("-sf", "--savefolder", nargs="?", const="Energy_spectrum", \
                        help="Specific folder for saving, as in Energy_spectrum") 
    parser.add_argument("-s","--save_eps", action="store_true", help="Decide whether to save in eps or not")            
    args = parser.parse_args()
    
    ### read general information from the folder
    
    path1 = 'eps_' + str(args.eps) + '/fp_' + str(args.fp) + '/areak_' + str(args.areak) + '/' 
    datafolder = args.folder + path1
    sim = read_write.read_sim_info(datafolder)
    
    ### read the data
    
    path2 = 'eps_' + str(args.eps) + '_fp_' + str(args.fp) + '_areak_' + str(args.areak) 
    analysisdatabase = '/usr/users/iff_th2/duman/Cells_in_LAMMPS/DATA/'
    analysisdatabase += args.savefolder + '/'
    analysisdata = analysisdatabase + args.savefolder + '_' + path2 + '.txt'
    x, y = read_write.read_2d_analysis_data(analysisdata)
        
    ### plot the data in the given time window
    
    plot_data(x, y, sim, args.savebase, args.savefolder, args.save_eps)
    
    return
    
##############################################################################

if __name__ == '__main__':
    main()    
    
##############################################################################
