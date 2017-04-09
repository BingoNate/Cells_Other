
""" Plot global lattice order per area constraint potential parameter"""

### example command line arguments: 
###    -a=100.0  
###    -fl=/local/duman/SIMULATIONS/Cells_in_LAMMPS/density_0.8/ 
###         -sb=/usr/users/iff_th2/duman/Cells_in_LAMMPS/PLOTS/
###             -sf=Glob_lattice_order

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
    
class Data:
    """ hold the simulation data as a function of two parameters"""
    
    def __init__(self, p1, p2, data, sim):
        self.p1 = p1
        self.p2 = p2
        self.data = data
        self.sim = sim
        
        return
        
##############################################################################

def plot_data(data, a, savebase, savefolder, save_eps):
    """ plot the data as a function of the chosen parameter
    z encapsulates the data and sim info as a function of two parameters"""

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
    
    ### collect the values to be plotted in a numpy array
    
    x = np.zeros((len(data)), dtype=np.float32)
    y = np.zeros((len(data)), dtype=np.float32)
    z = np.zeros((len(data)), dtype=np.float32)
    
    for j, key in enumerate(data.keys()):
            
        #e, f = key                    # parameters
        #sim = data[key].sim           # simulation info for this parameter set
        z[j] = data[key].data          # analysis data
        x[j] = data[key].p1            # x value in the plot
        y[j] = data[key].p2            # y value in the plot 
        

    ### get the size and color code
    
    print z
    size = 1e+8 * np.pi * z**2
    quant_steps = 2056
    norm = mpl.colors.Normalize(vmin=0.0, vmax=max(z))       
    line0 = ax0.scatter(x, y, c=z, \
                        cmap=plt.cm.get_cmap('jet',quant_steps), \
                        edgecolors='None', alpha=1.0, vmin=0.0, vmax=max(z), \
                        norm=norm, rasterized=True) 
    fig.colorbar(line0, ax=ax0)
    

    ### title
    
#    ax0.set_title("$t/\\tau_{D}$ = " + "{0:.2f}".format(time/sim.tau_D) + \
#        ", $t/\\tau_{A}$ = " + "{0:.2f}".format(time/sim.tau_A), fontsize=30)
    
    ### labels
        
    ax0.set_xlabel("$\\epsilon$", fontsize=40)
    ax0.set_ylabel("$f_{m}$", fontsize=40)

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

    savepath1 = base + savefolder + "_per_AREAK_" + str(a) + ".png"
    if save_eps:
        savepath2 = base + savefolder + "_per_AREAK_" + str(a) + ".eps"            
    plt.savefig(savepath1, dpi=200, bbox_inches='tight', pad_inches=0.08)
    if save_eps:
        plt.savefig(savepath2, dpi=200, bbox_inches='tight', pad_inches=0.08)        
    fig.clf()                

        
    return
        
##############################################################################

def main():

    ### get the data folder
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--areak", type=float, help="Area constraint potential strength")    
    parser.add_argument("-fl", "--folder", nargs="?", \
                        const='/local/duman/SIMULATIONS/Cells_in_LAMMPS/density_0.8/', \
                        help="Folder containing data, as in /local/duman/SIMULATIONS/Cells_in_LAMMPS/density_0.8/")    
    parser.add_argument("-sb", "--savebase", nargs="?", \
                        const = "/usr/users/iff_th2/duman/Cells_in_LAMMPS/PLOTS/", \
                        help="Folder to save the data, as in /usr/users/iff_th2/duman/Cells_in_LAMMPS/PLOTS/") 
    parser.add_argument("-sf", "--savefolder", nargs="?", const="Glob_lattice_order", \
                        help="Specific folder for saving, as in Glob_lattice_order")   
    parser.add_argument("-s","--save_eps", action="store_true", help="Decide whether to save in eps or not")            
    args = parser.parse_args()
    
    ### accummulate data as a function of areak -with a dictionary-

    analysisdatabase = '/usr/users/iff_th2/duman/Cells_in_LAMMPS/DATA/'  
    analysisdatabase += args.savefolder + '/'    
    eps = [0.05, 0.5, 1.0, 5.0, 10.0, 20.0]
    fp = [0.0, 0.5, 1.0, 3.0, 5.0, 10.0]
    #areak = [1.0, 10.0, 100.0]
    alldata = {}       # carries the data per parameter set

    for e in eps:
        for f in fp:
                      
            key = (e, f)
            datafolder, analysisfile = read_write.gen_folders( \
                e, f, args.areak, args.savefolder, 
                args.folder, analysisdatabase)
            sim_info = read_write.read_sim_info(datafolder)
            data_info = read_write.read_single_analysis_data(analysisfile)
            alldata[key] = Data(e, f, data_info, sim_info)
                
        
    ### plot the data as a function of the parameter
    
    plot_data(alldata, args.areak, args.savebase, args.savefolder, args.save_eps)

    return
    
##############################################################################

if __name__ == '__main__':
    main()    
    
##############################################################################
