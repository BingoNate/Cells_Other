
""" Generate a movie of the beads of cells"""

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
        
##############################################################################

def plot_frames(beads, sim, ti, tf, savebase, save_eps):
    """ plot frames within the specified time window"""
    
    ### normalize variables for plotting purposes
    
    lx = sim.lx/sim.bl
    ly = sim.ly/sim.bl
        
    ### set general plot properties

    savebase += 'eps_' + str(sim.eps) + '_fp_' + str(sim.fp) + '_areak_' + str(sim.areak) + '/'
    os.system("mkdir -p " + savebase)
    quant_steps = 2056
    norm = mpl.colors.Normalize(vmin=0, vmax=sim.ncells)  
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
    
    ### plot the frames
    
    for step in range(ti, tf):
        
        ### normalize variables for plotting purposes
        
        beads.xi[step, 0, :] /= sim.bl
        beads.xi[step, 1, :] /= sim.bl
        
        time = step*sim.dt
        print 'Step / Total : ', step, tf
        
        ### plot 

        subp = misc_tools.Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x) 
        ax0 = subp.addSubplot()
        
        line0 = ax0.scatter(beads.xi[step, 0, :], beads.xi[step, 1, :], s=1, c=beads.cid, \
                            cmap=plt.cm.get_cmap('jet',quant_steps), \
                            edgecolors='None', alpha=0.7, vmin=0, vmax=sim.ncells, \
                            norm=norm, rasterized=True)
        
        ax0.axis('scaled')
    
        ### title
        
        ax0.set_title("$t/\\tau_{D}$ = " + "{0:.2f}".format(time/sim.tau_D) + \
            ", $t/\\tau_{A}$ = " + "{0:.2f}".format(time/sim.tau_A), fontsize=30)
        
        ### labels
            
        ax0.set_xlabel("$x/r_{0}$", fontsize=40)
        ax0.set_ylabel("$y/r_{0}$", fontsize=40)

        ### limits

        ax0.set_xlim((downlim, uplim))
        ax0.set_ylim((downlim, uplim))
        
        ### ticks
        
        ax0.xaxis.set_ticks(np.linspace(0, uplim, num_ticks, endpoint=True))
        ax0.yaxis.set_ticks(np.linspace(0, uplim, num_ticks, endpoint=True))
        ax0.tick_params(axis='both', which='major', labelsize=30)
        
        ### save

        savepath1 = savebase + "frame-" + "{0:05d}".format(int(step)) + ".png"
        if save_eps:
            savepath2 = savebase + "frame-" + "{0:05d}".format(int(step)) + ".eps"
            
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
                        const = "/usr/users/iff_th2/duman/Cells_in_LAMMPS/MOVIES/", \
                        help="Folder to save the data, as in /usr/users/iff_th2/duman/Cells_in_LAMMPS/MOVIES/")     
    parser.add_argument("-ti","--init_time", nargs="?", const=10, type=int, \
                        help="First frame of the video (in terms of frame number), you can also leave it empty")
    parser.add_argument("-tf","--fin_time", nargs="?", const=2000, type=int, \
                        help="Last frame of the video (in terms of frame number), you can also leave it empty")
    parser.add_argument("-s","--save_eps", action="store_true", help="Decide whether to save in eps or not")            
    args = parser.parse_args()
    
    ### read the data and general information from the folder
    
    sim, cells, beads = read_write.read_h5_file(args.folder)
    beads.get_img_pos(sim.lx)
    print "folder = ", args.folder
        
    ### plot the data in the given time window
    
    plot_frames(beads, sim, args.init_time, args.fin_time, args.savebase, args.save_eps)
    
    return
    
##############################################################################

if __name__ == '__main__':
    main()    
    
##############################################################################
