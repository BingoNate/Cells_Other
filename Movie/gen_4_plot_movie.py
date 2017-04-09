
""" Generate a detailed movie for cells"""

### example command line arguments: 
###    -fl=/local/duman/SIMULATIONS/Cells_in_LAMMPS/.../ 
###         -sb=/usr/users/iff_th2/duman/Cells_in_LAMMPS/MOVIES/detailed/
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
from collections import deque
import seaborn as sns

sns.set(style="white",context='paper',
        font_scale=1.2,font="Open Sans",
        rc={'mathtext.default': 'regular','font.size': 30, 
            'font.family': 'sans',"figure.dpi":300,
            "xtick.major.size": 8, "ytick.major.size": 8,
            'grid.linestyle': '--'})

#mpl.rc({'mathtext.default':'regular','font.size':30, 
#            'figure.figsize':(5,5),
#            'font.family':'sans', "figure.dpi":300,
#            "xtick.major.size": 8, "ytick.major.size": 8})

#sns.set_palette(sns.color_palette("cubehelix"))

##############################################################################
  
def cell2bead_vertical_pos(beads, cells, sim, nslices, ti):
    """ assign colors from cell indices to bead indices"""
    
    ### assignment based on the vertical positions of cells 
    
    slice_width = sim.lx/sim.r_avg/nslices
    vertical_pos_colors_per_bead = np.zeros((sim.nbeads), dtype=np.int32) - 1
    k = 0
    for i in xrange(sim.ncells):
        for j in xrange(sim.nbpc[i]):
            vertical_pos_colors_per_bead[k] = \
                np.round((cells.xi[ti, 1, i]/sim.r_avg) / slice_width)
            k += 1

    return vertical_pos_colors_per_bead
    
##############################################################################

def plot_frames(beads, cells, sim, ti, tf, savebase, save_eps):
    """ plot frames within the specified time window"""
    
    ### normalize variables for plotting purposes
    
    lx = sim.lx/sim.r_avg
    ly = sim.ly/sim.r_avg
        
    ### set general plot properties

    savebase += 'eps_' + str(sim.eps) + '_fp_' + str(sim.fp) + \
        '_areak_' + str(sim.areak) + '/'
    os.system("mkdir -p " + savebase)
    quant_steps = 2056
    
    # limits
    full_box_downlim = -2
    full_box_uplim = lx+2
    full_box_ticks = [0, 35, 70, 105, 135]

    half_box_downlim = 43
    half_box_uplim = 92
    half_box_ticks = [45, 90]

    num_ticks = 5
    
    ax_len = 2.0                         # Length of one subplot square box
    ax_b = 0.0                           # Beginning/offset of the subplot in the box
    ax_sep = 0.2                         # Separation length between two subplots
    total_subplots_in_x = 2              # Total number of subplots    
    fig = plt.figure()
    subp = misc_tools.Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x) 
    ax2 = subp.addSubplot()
    ax3 = subp.addSubplot()
    ax0 = subp.addSubplot()
    ax1 = subp.addSubplot()
    
    ### set subplot properties 
    
    nslices = sim.ncells
#    nslices = 100                          # Number of slices to take in the y dir
    norm_ax0 = mpl.colors.Normalize(vmin=0, vmax=nslices)  
#    vertical_pos_colors_per_bead = cell2bead_vertical_pos(beads, cells, sim, nslices, ti)
#    cmap_ax0 = sns.dark_palette((260, 75, 60), n_colors=nslices, input="husl", as_cmap=True)
#    cmap_ax0 = sns.diverging_palette(255, 133, l=60, n=nslices, center="dark", as_cmap=True)
    cmap_ax0 = plt.cm.get_cmap('jet', quant_steps)
    
    ### plot the frames
    
    comx = deque()
    comy = deque()
    nfreq = 5
    
    for step in range(ti, tf):

        time = step*sim.dt
        print 'Step / Total : ', step, tf        
        
        ### normalize variables for plotting purposes
        
        beads.xi[step, 0, :] /= sim.r_avg
        beads.xi[step, 1, :] /= sim.r_avg
        cells.xi[step, 0, :] /= sim.r_avg
        cells.xi[step, 1, :] /= sim.r_avg

        ### keep the center of mass trajectory
        
        comx.append(cells.xi[step, 0, :])
        comy.append(cells.xi[step, 1, :])
        
        if step > ti+nfreq:
            comx.popleft()
            comy.popleft()
                        
        ### plot 

        subp = misc_tools.Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x) 
        
        ax2 = subp.addSubplot()
        ax3 = subp.addSubplot()
        ax0 = subp.addSubplot()
        ax1 = subp.addSubplot()
        
        ### AX0

        #ax0 = subp.addSubplot()
        
        line0 = ax0.scatter(beads.xi[step, 0, :], beads.xi[step, 1, :], s=4.0, \
                            c=beads.cid, \
                            cmap=cmap_ax0, \
                            edgecolors='None', alpha=0.8, vmin=0, vmax=nslices, \
                            norm=norm_ax0, rasterized=True)
        
#        ax0.quiver(cells.xi[step, 0, :], cells.xi[step, 1, :], \
#                   10*np.cos(cells.pol[step]), 10*np.sin(cells.pol[step]), \
#                   headwidth=5, headlength=6, headaxislength=3.5, alpha=0.7)
        
        ax0.axis('scaled')
    
        ### title
        
        ax0.set_title(r"$t/\tau_{D}$ = " + "{0:.2f}".format( time/sim.tau_D) + \
            r", $t/\tau_{A}$ = " + "{0:.2f}".format(time/sim.tau_A), fontsize=40)
        
        ### labels
            
        ax0.set_ylabel(r"$y/R$", fontsize=40)

        ### limits

        ax0.set_xlim((full_box_downlim, full_box_uplim))
        ax0.set_ylim((full_box_downlim, full_box_uplim))
        
        ### ticks
        
        plt.setp(ax0.get_xticklabels(),visible=False)  
        ax0.yaxis.set_ticks(full_box_ticks)
        ax0.tick_params(axis='both', which='major', labelsize=40)
        
        
        
        
        ### AX1
        
        #ax1 = subp.addSubplot()
        
        line1 = ax1.scatter(cells.xi[step, 0, :], cells.xi[step, 1, :], s=5.0, \
                            c=np.arange(nslices), \
                            cmap=cmap_ax0, \
                            edgecolors='None', alpha=1.0, vmin=0, vmax=nslices, \
                            norm=norm_ax0, rasterized=True)
        
        line2 = ax1.scatter(list(comx), list(comy), s=4.0, \
                            c=np.ones(np.shape(list(comx)))*np.arange(nslices), \
                            cmap=cmap_ax0, \
                            edgecolors='None', alpha=0.5, vmin=0, vmax=nslices, \
                            norm=norm_ax0, rasterized=True)

        ax1.axis('scaled')
    
        ### title
        
        ax1.set_title(r"Center of Mass Trajectories", fontsize=40)

        ### limits

        ax1.set_xlim((full_box_downlim, full_box_uplim))
        ax1.set_ylim((full_box_downlim, full_box_uplim))
        
        ### ticks

        plt.setp(ax1.get_xticklabels(),visible=False)       
        plt.setp(ax1.get_yticklabels(),visible=False)                
#        ax1.xaxis.set_ticks([])
#        ax1.yaxis.set_ticks([])
        ax1.tick_params(axis='both', which='major', labelsize=40)

        
        
        ### AX2
        
        #ax2 = subp.addSubplot()
        
        line3 = ax2.quiver(cells.xi[step, 0, :], cells.xi[step, 1, :], \
                   100*np.cos(cells.pol[step]), 100*np.sin(cells.pol[step]), \
                   headwidth=5, headlength=6, headaxislength=3.5, alpha=0.7)
        

        ax2.axis('scaled')
    
        ### title
        
        ax2.set_title(r"Polarities", fontsize=40)

        ### labels
            
        ax2.set_xlabel(r"$x/R$", fontsize=40)
        ax2.set_ylabel(r"$y/R$", fontsize=40)
        
        ### limits

        ax2.set_xlim((full_box_downlim, full_box_uplim))
        ax2.set_ylim((full_box_downlim, full_box_uplim))
        
        ### ticks
        
        ax2.xaxis.set_ticks(full_box_ticks)
        ax2.yaxis.set_ticks(full_box_ticks)
        ax2.tick_params(axis='both', which='major', labelsize=40)

        
        
        ### AX3
        
        #ax3 = subp.addSubplot()
        
        line3 = ax3.quiver(cells.xi[step, 0, :], cells.xi[step, 1, :], \
                   100*np.cos(cells.pol[step]), 100*np.sin(cells.pol[step]), \
                   headwidth=5, headlength=6, headaxislength=3.5, alpha=0.7)
        

        ax3.axis('scaled')
    
        ### title
        
        ax3.set_title(r"Polarities", fontsize=40)

        ### limits

        ### labels
            
        ax3.set_xlabel(r"$x/R$", fontsize=40)
        #ax3.set_ylabel(r"$y/R$", fontsize=40)

        ### limits

        ax3.set_xlim((full_box_downlim, full_box_uplim))
        ax3.set_ylim((full_box_downlim, full_box_uplim))
        
        ### ticks
        
        ax3.xaxis.set_ticks(full_box_ticks)
        plt.setp(ax3.get_yticklabels(),visible=False)  
        ax3.tick_params(axis='both', which='major', labelsize=40)
         
        
        
        ### save

        savepath1 = savebase + "frame-" + "{0:05d}".format(int(step)) + ".png"
        if save_eps:
            savepath2 = savebase + "frame-" + "{0:05d}".format(int(step)) + ".eps"
            
        plt.savefig(savepath1, dpi=300, bbox_inches='tight', pad_inches=0.08)
        if save_eps:
            plt.savefig(savepath2, dpi=300, bbox_inches='tight', pad_inches=0.08)        
        fig.clf()                
        
    return
        
##############################################################################

def main():

    ### get the data folder
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-fl", "--folder", \
                        help="Folder containing data, as in /local/duman/SIMULATIONS/Cells_in_LAMMPS/.../")
    parser.add_argument("-sb", "--savebase", nargs="?", \
                        const = "/usr/users/iff_th2/duman/Cells_in_LAMMPS/MOVIES/detailed/", \
                        help="Folder to save the data, as in /usr/users/iff_th2/duman/Cells_in_LAMMPS/MOVIES/detailed/")     
    parser.add_argument("-ti","--init_time", nargs="?", const=10, type=int, \
                        help="First frame of the video (in terms of frame number), you can also leave it empty")
    parser.add_argument("-tf","--fin_time", nargs="?", const=12, type=int, \
                        help="Last frame of the video (in terms of frame number), you can also leave it empty")
    parser.add_argument("-s","--save_eps", action="store_true", help="Decide whether to save in eps or not")            
    args = parser.parse_args()
    
    ### read the data and general information from the folder
    
    sim, cells, beads = read_write.read_h5_file(args.folder)
    beads.get_img_pos(sim.lx)
    cells.get_img_pos(sim.lx)
    print "folder = ", args.folder
        
    ### plot the data in the given time window
    
    plot_frames(beads, cells, sim, args.init_time, args.fin_time, args.savebase, args.save_eps)
    
    return
    
##############################################################################

if __name__ == '__main__':
    main()    
    
##############################################################################
