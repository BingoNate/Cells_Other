
""" Plot the vorticity field per simulation"""

### example command line arguments: 
###    -fl=/local/duman/SIMULATIONS/Cells_in_LAMMPS/.../ 
###         -sb=/usr/users/iff_th2/duman/Cells_in_LAMMPS/PLOTS/Vorticity_field/
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
from numpy import ma
from matplotlib import cbook
from matplotlib.colors import Normalize
from scipy import interpolate

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

class MidPointNorm(Normalize):    
    def __init__(self, midpoint=0, vmin=None, vmax=None, clip=False):
        Normalize.__init__(self,vmin, vmax, clip)
        self.midpoint = midpoint

    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)

        self.autoscale_None(result)
        vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

        if not (vmin < midpoint < vmax):
            raise ValueError("midpoint must be between maxvalue and minvalue.")       
        elif vmin == vmax:
            result.fill(0) # Or should it be all masked? Or 0.5?
        elif vmin > vmax:
            raise ValueError("maxvalue must be bigger than minvalue")
        else:
            vmin = float(vmin)
            vmax = float(vmax)
            if clip:
                mask = ma.getmask(result)
                result = ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                  mask=mask)

            # ma division is very slow; we can take a shortcut
            resdat = result.data

            #First scale to -1 to 1 range, than to from 0 to 1.
            resdat -= midpoint            
            resdat[resdat>0] /= abs(vmax - midpoint)            
            resdat[resdat<0] /= abs(vmin - midpoint)

            resdat /= 2.
            resdat += 0.5
            result = ma.array(resdat, mask=result.mask, copy=False)                

        if is_scalar:
            result = result[0]            
        return result

    def inverse(self, value):
        if not self.scaled():
            raise ValueError("Not invertible until scaled")
        vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

        if cbook.iterable(value):
            val = ma.asarray(value)
            val = 2 * (val-0.5)  
            val[val>0]  *= abs(vmax - midpoint)
            val[val<0] *= abs(vmin - midpoint)
            val += midpoint
            return val
        else:
            val = 2 * (val - 0.5)
            if val < 0: 
                return  val*abs(vmin-midpoint) + midpoint
            else:
                return  val*abs(vmax-midpoint) + midpoint

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
  
def calc_displacement_magnitudes(cells, step, delta, sim):
    """ calculate the displacements and displacement magnitudes"""
    
    ### NOTE THAT displacements are calculated from unwrapped coordinates
    
    dx = np.zeros((sim.ncells), dtype=np.float32)
    dy = np.zeros((sim.ncells), dtype=np.float32)
    dr = np.zeros((sim.ncells), dtype=np.float32)
    
    dx = (cells.xu[step+delta, 0, :] - cells.xu[step, 0, :])/sim.dt
    dy = (cells.xu[step+delta, 1, :] - cells.xu[step, 1, :])/sim.dt
    dr = np.sqrt(dx**2 + dy**2)
    
    return dx, dy, dr
    
##############################################################################

def plot_frames(data, beads, cells, sim, ti, tf, savebase, save_eps):
    """ plot frames within the specified time window"""
    
    ### normalize variables for plotting purposes
    
    lx = sim.lx/sim.r_avg
    ly = sim.ly/sim.r_avg
    
    ### vorticity information
    
    steps, xbins, ybins, w, vx, vy, v = data
    nwbins = int(max(xbins)+1)
    xlin = np.linspace(0., lx, nwbins)
    ylin = np.linspace(0., ly, nwbins)
    xgrid, ygrid = np.meshgrid(xlin, ylin)    
      
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
    
    ax_len = 2.2                        # Length of one subplot square box
    ax_b = 0.01                         # Beginning/offset of the subplot in the box
    ax_sep = 0.15                        # Separation length between two subplots
    total_subplots_in_x = 2             # Total number of subplots    
    fig = plt.figure()
    subp = misc_tools.Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x) 
    multi = True
    ax1 = subp.addSubplot(multi)    
    ax3 = subp.addSubplot(multi) 
    ax0 = subp.addSubplot(multi)           
    ax2 = subp.addSubplot(multi)    
    
    ### set subplot properties 
    
    nslices = sim.ncells
    norm_ax0 = mpl.colors.Normalize(vmin=0, vmax=nslices)  
    cmap_ax0 = plt.cm.get_cmap('jet', quant_steps)
    
    ### plot the frames
    
    comx = deque()
    comy = deque()
    ntrace = 6
    delta = 3

    ### normalize central box positions for plotting purposes
    
    beads.xi /= sim.r_avg
    cells.xi /= sim.r_avg
    
    if tf+delta > sim.nsteps:
        tf -= delta
        
    for step in range(ti, tf):

        time = step*sim.dt
        print 'Step / Total : ', step, tf        

        ### calculate the displacement magnitudes
        
        dx, dy, dr = calc_displacement_magnitudes(cells, step, delta, sim)
        dx /= sim.r_avg
        dy /= sim.r_avg
        dr /= sim.r_avg
        #norm_ax1 = mpl.colors.Normalize(vmin=min(dr), vmax=max(dr))
               
        ### keep the center of mass trajectory
        
        comx.append(cells.xi[step, 0, :])
        comy.append(cells.xi[step, 1, :])
        
        if step > ti+ntrace:
            comx.popleft()
            comy.popleft()
            
        ### get the vorticity and velocity information
        
        vs = v[steps==step]
        vmean = np.mean(vs)
        vn = np.zeros((nwbins, nwbins), dtype=np.float32)   
        
        ws = w[steps==step]
        wmean = np.mean(np.abs(ws))
        wn = np.zeros((nwbins, nwbins), dtype=np.float32)
        
        for xi, yi in zip(xbins, ybins):
            xi = int(xi)
            yi = int(yi)
            wn[xi, yi] = ws[xi*nwbins+yi]
            vn[xi, yi] = vs[xi*nwbins+yi]     

        vn /= vmean       
        vmin = np.min(vn)
        vmax = np.max(vn)
        norm_ax2 = mpl.colors.Normalize(vmin=vmin, vmax=vmax)       
        
        wn /= wmean
        wmin = np.min(wn)
        wmax = np.max(wn)
        norm_ax3 = mpl.colors.Normalize(vmin=wmin, vmax=wmax)   
            
        
        ### plot 

        subp = misc_tools.Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x) 
        ax1 = subp.addSubplot(multi)    
        ax3 = subp.addSubplot(multi) 
        ax0 = subp.addSubplot(multi)           
        ax2 = subp.addSubplot(multi)    
        
        text = r"$t/\tau_{D}$ = " + "{0:.2f}".format( time/sim.tau_D) + \
            r", $t/\tau_{A}$ = " + "{0:.2f}".format(time/sim.tau_A)
            
            
        
        ### AX0
        
        line0 = ax0.scatter(beads.xi[step, 0, :], beads.xi[step, 1, :], s=4.0, \
                            c=beads.cid, \
                            cmap=cmap_ax0, \
                            edgecolors='None', alpha=1.0, vmin=0, vmax=nslices, \
                            norm=norm_ax0, rasterized=True)
        
        
        ax0.axis('scaled')
    
        ### labels

        ax0.set_ylabel(r"$y/R$", fontsize=40)

        ### limits

        ax0.set_xlim((full_box_downlim, full_box_uplim))
        ax0.set_ylim((full_box_downlim, full_box_uplim))
        
        ### ticks
        
        ax0.xaxis.set_ticks(full_box_ticks)
        ax0.yaxis.set_ticks(full_box_ticks)   
        plt.setp(ax0.get_xticklabels(),visible=False)                        
        ax0.tick_params(axis='both', which='major', labelsize=40)
        
        
        
        ### AX1
                
        line1 = ax1.scatter(cells.xi[step, 0, :], cells.xi[step, 1, :], s=6.0, \
                            c=np.arange(nslices), \
                            #c=dr,
                            cmap=cmap_ax0, \
                            edgecolors='None', alpha=1.0, vmin=0, vmax=nslices, \
                            norm=norm_ax0, rasterized=True)
        
        line2 = ax1.scatter(list(comx), list(comy), s=5.0, \
                            c=np.ones(np.shape(list(comx)))*np.arange(nslices), \
                            #c=np.ones(np.shape(list(comx)))*dr,
                            cmap=cmap_ax0, \
                            edgecolors='None', alpha=0.5, vmin=0, vmax=nslices, \
                            norm=norm_ax0, rasterized=True)

        ax1.axis('scaled')

        ### labels

        ax1.set_xlabel(r"$x/R$", fontsize=40)
        ax1.set_ylabel(r"$y/R$", fontsize=40)
        
        ### limits

        ax1.set_xlim((full_box_downlim, full_box_uplim))
        ax1.set_ylim((full_box_downlim, full_box_uplim))
        
        ### ticks

        ax1.xaxis.set_ticks(full_box_ticks)
        ax1.yaxis.set_ticks(full_box_ticks)
        #plt.setp(ax1.get_yticklabels(),visible=False)                
        ax1.tick_params(axis='both', which='major', labelsize=40)


        
        ### AX2
                
        line3 = ax2.contourf(xgrid, ygrid, vn.transpose(), cmap=cmap_ax0,
                             #norm=norm_ax2, vmin=vmin, vmax=vmax)  
                             )
        
        line4 = ax2.quiver(cells.xi[step, 0, :], cells.xi[step, 1, :], \
                   dx, dy, \
                   headwidth=5, headlength=6, headaxislength=3.5, alpha=0.7)    
        
        ax2.axis('scaled')
        
        cax2 = plt.colorbar(line3, ax=ax2)
        #plt.colorbar(line5, cax=cax3, ticks=[])
        #cax2.ax.set_yticks([0.0, 0.3, 0.6, 0.9, 1.2 ,1.5])
        #cax3.set_yticklabels(['0', '0.7'])  
        cax2.ax.tick_params(labelsize=40) 
        cax2.ax.set_title(r"$|v|/<|v|>$",fontsize=40)
        
        ### labels

        #ax2.set_xlabel(r"$x/R$", fontsize=40)
        
        ### limits

        ax2.set_xlim((full_box_downlim, full_box_uplim))
        ax2.set_ylim((full_box_downlim, full_box_uplim))
        
        ### ticks

        #ax2.xaxis.set_ticks(full_box_ticks)
        ax2.xaxis.set_ticks(full_box_ticks)
        ax2.yaxis.set_ticks(full_box_ticks)           
        plt.setp(ax2.get_xticklabels(),visible=False)                        
        plt.setp(ax2.get_yticklabels(),visible=False)                
        ax2.tick_params(axis='both', which='major', labelsize=40)
        
        
        
        ### AX3 

        line6 = ax3.contourf(xgrid, ygrid, wn.transpose(), cmap=cmap_ax0,
                             #norm=norm_ax3, vmin=wmin, vmax=wmax)
                            )
        
        line7 = ax3.quiver(cells.xi[step, 0, :], cells.xi[step, 1, :], \
                   dx, dy, \
                   headwidth=5, headlength=6, headaxislength=3.5, alpha=1.0)   
        
        ax3.axis('scaled')
        
        cax3 = plt.colorbar(line6, ax=ax3)
        #plt.colorbar(line5, cax=cax3, ticks=[])
        #cax3.ax.set_yticks([-6.0, -3.0, 0.0 , 3.0, 6.0])
        #cax3.set_yticklabels(['0', '0.7'])  
        cax3.ax.tick_params(labelsize=40)    
        cax3.ax.set_title(r"$\omega/<|\omega|>$",fontsize=40)
        
        ### labels

        ax3.set_xlabel(r"$x/R$", fontsize=40)
        
        ### limits

        ax3.set_xlim((full_box_downlim, full_box_uplim))
        ax3.set_ylim((full_box_downlim, full_box_uplim))
        
        ### ticks

        ax3.xaxis.set_ticks(full_box_ticks)
        ax3.yaxis.set_ticks(full_box_ticks)        
        plt.setp(ax3.get_yticklabels(),visible=False)                
        ax3.tick_params(axis='both', which='major', labelsize=40)
        
        
        ### text
        
        plt.figtext(subp.xbeg-1.1*ax_sep, subp.ybeg+ax_len+0.1*ax_sep, text, fontsize=40)
        
        
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
    
def plot_vorticity(data, beads, cells, sim, ti, tf, savebase, save_eps):
    """ plot vorticity field within the specified time window"""
    
    ### normalize variables for plotting purposes
    
    lx = sim.lx/sim.r_avg
    ly = sim.ly/sim.r_avg
    
    ### vorticity information
    
    steps, xbins, ybins, w, vx, vy, v = data
    nwbins = int(max(xbins)+1)
    xlin = np.linspace(0., lx, nwbins)
    ylin = np.linspace(0., ly, nwbins)
    xgrid, ygrid = np.meshgrid(xlin, ylin)    
      
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
    
    ax_len = 2.2                        # Length of one subplot square box
    ax_b = 0.01                         # Beginning/offset of the subplot in the box
    ax_sep = 0.15                       # Separation length between two subplots
    total_subplots_in_x = 1             # Total number of subplots    
    fig = plt.figure()
    subp = misc_tools.Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x) 
    multi = True
    ax0 = subp.addSubplot(multi)    
    ax1 = subp.addSubplot(multi) 
    
    ### set subplot properties 
    
    nslices = sim.ncells
    norm_ax0 = mpl.colors.Normalize(vmin=0, vmax=nslices)  
    cmap_ax0 = plt.cm.get_cmap('jet', quant_steps)
    
    ### plot the frames
    
    comx = deque()
    comy = deque()
    ntrace = 5
    delta = 5

    ### normalize central box positions for plotting purposes
    
    beads.xi /= sim.r_avg
    cells.xi /= sim.r_avg
    
    if tf+delta > sim.nsteps:
        tf -= delta
        
    for step in range(ti, tf):

        time = step*sim.dt
        print 'Step / Total : ', step, tf        

        ### calculate the displacement magnitudes
        
        dx, dy, dr = calc_displacement_magnitudes(cells, step, delta, sim)
        dx /= sim.r_avg
        dy /= sim.r_avg
        dr /= sim.r_avg
        #norm_ax1 = mpl.colors.Normalize(vmin=min(dr), vmax=max(dr))
            
        ### get the vorticity and velocity information
        
        vs = v[steps==step]
        vmean = np.mean(vs)
        vn = np.zeros((nwbins, nwbins), dtype=np.float32)   
        
        ws = w[steps==step]
        wmean = np.mean(np.abs(ws))
        wn = np.zeros((nwbins, nwbins), dtype=np.float32)
        
        for xi, yi in zip(xbins, ybins):
            xi = int(xi)
            yi = int(yi)
            wn[xi, yi] = ws[xi*nwbins+yi]
            vn[xi, yi] = vs[xi*nwbins+yi]     

        vn /= vmean       
        vmin = np.min(vn)
        vmax = np.max(vn)
        norm_ax2 = mpl.colors.Normalize(vmin=vmin, vmax=vmax)       
        
        wn /= wmean
        wmin = np.min(wn)
        wmax = np.max(wn)
        norm_ax3 = mpl.colors.Normalize(vmin=wmin, vmax=wmax)   
            
        
        ### plot 

        subp = misc_tools.Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x) 
        ax0 = subp.addSubplot(multi)    
        ax1 = subp.addSubplot(multi)  
        
        
        text = r"$t/\tau_{D}$ = " + "{0:.2f}".format( time/sim.tau_D) + \
            r", $t/\tau_{A}$ = " + "{0:.2f}".format(time/sim.tau_A)

        
        ### AX1
                
        line2 = ax1.contourf(xgrid, ygrid, vn.transpose(), cmap=cmap_ax0,
                             #norm=norm_ax2, vmin=vmin, vmax=vmax)  
                             )
        
        line3 = ax1.quiver(cells.xi[step, 0, :], cells.xi[step, 1, :], \
                   dx, dy, \
                   headwidth=5, headlength=6, headaxislength=3.5, alpha=0.7)    
        
        ax1.axis('scaled')
        
        cax1 = plt.colorbar(line2, ax=ax1)
        #plt.colorbar(line5, cax=cax3, ticks=[])
        #cax2.ax.set_yticks([0.0, 0.3, 0.6, 0.9, 1.2 ,1.5])
        #cax3.set_yticklabels(['0', '0.7'])  
        cax1.ax.tick_params(labelsize=40) 
        cax1.ax.set_title(r"$|v|/<|v|>$",fontsize=40)
        
        ### labels

        ax1.set_ylabel(r"$y/R$", fontsize=40)
        
        ### limits

        ax1.set_xlim((full_box_downlim, full_box_uplim))
        ax1.set_ylim((full_box_downlim, full_box_uplim))
        
        ### ticks

        ax1.xaxis.set_ticks(full_box_ticks)
        ax1.yaxis.set_ticks(full_box_ticks)           
        plt.setp(ax1.get_xticklabels(),visible=False)                        
        #plt.setp(ax1.get_yticklabels(),visible=False)                
        ax1.tick_params(axis='both', which='major', labelsize=40)
        
        
        
        ### AX0

        line0 = ax0.contourf(xgrid, ygrid, wn.transpose(), cmap=cmap_ax0,
                             #norm=norm_ax3, vmin=wmin, vmax=wmax)
                            )
        
        line1 = ax0.quiver(cells.xi[step, 0, :], cells.xi[step, 1, :], \
                   dx, dy, \
                   headwidth=5, headlength=6, headaxislength=3.5, alpha=1.0)   
        
        ax0.axis('scaled')
        
        cax0 = plt.colorbar(line0, ax=ax0)
        #plt.colorbar(line5, cax=cax3, ticks=[])
        #cax3.ax.set_yticks([-6.0, -3.0, 0.0 , 3.0, 6.0])
        #cax3.set_yticklabels(['0', '0.7'])  
        cax0.ax.tick_params(labelsize=40)    
        cax0.ax.set_title(r"$\omega/<|\omega|>$",fontsize=40)
        
        ### labels
        
        ax0.set_xlabel(r"$x/R$", fontsize=40)
        ax0.set_ylabel(r"$y/R$", fontsize=40)
        
        ### limits

        ax0.set_xlim((full_box_downlim, full_box_uplim))
        ax0.set_ylim((full_box_downlim, full_box_uplim))
        
        ### ticks

        ax0.xaxis.set_ticks(full_box_ticks)
        ax0.yaxis.set_ticks(full_box_ticks)        
        #plt.setp(ax0.get_yticklabels(),visible=False)                
        ax0.tick_params(axis='both', which='major', labelsize=40)
        
        
        ### text
        
        #plt.figtext(subp.xbeg-1.1*ax_sep, subp.ybeg+ax_len+0.1*ax_sep, text, fontsize=40)
        
        
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
    parser.add_argument("-e", "--eps", type=float, nargs="?", const=-1, \
                        help="Strength of LJ potential")
    parser.add_argument("-f", "--fp", type=float, nargs="?", const=-1, \
                        help="Propulsion force")
    parser.add_argument("-a", "--areak", type=float, nargs="?", const=-1, \
                        help="Area constraint potential strength")
    parser.add_argument("-fl", "--folder", nargs="?", \
                        const='/local/duman/SIMULATIONS/Cells_in_LAMMPS/density_0.8/', \
                        help="Folder containing data, as in /local/duman/SIMULATIONS/Cells_in_LAMMPS/density_0.8/")    
    parser.add_argument("-sb", "--savebase", nargs="?", \
                        const = "/usr/users/iff_th2/duman/Cells_in_LAMMPS/PLOTS/", \
                        help="Folder to save the data, as in /usr/users/iff_th2/duman/Cells_in_LAMMPS/PLOTS/") 
    parser.add_argument("-sf", "--savefolder", nargs="?", const="Vorticity_field", \
                        help="Specific folder for saving, as in Vorticity_field")      
    parser.add_argument("-ti","--init_time", nargs="?", const=10, type=int, \
                        help="First frame of the video (in terms of frame number), you can also leave it empty")
    parser.add_argument("-tf","--fin_time", nargs="?", const=1000, type=int, \
                        help="Last frame of the video (in terms of frame number), you can also leave it empty")
    parser.add_argument("-s","--save_eps", action="store_true", help="Decide whether to save in eps or not")            
    args = parser.parse_args()
    
    ### read the data and general information from the folder
    
    folder = args.folder + "eps_" + str(args.eps) + "/fp_" + str(args.fp) + \
        "/areak_" + str(args.areak) + "/"
    sim, cells, beads = read_write.read_h5_file(folder)
    beads.get_img_pos(sim.lx)
    cells.get_img_pos(sim.lx)
    print "folder = ", args.folder
    
    ### read the vorticity data

    analysisdatabase = '/usr/users/iff_th2/duman/Cells_in_LAMMPS/DATA/'  
    analysisdatabase += args.savefolder + '/'     
    datafolder, analysisfile = read_write.gen_folders(args.eps, args.fp, args.areak, \
                                                      args.savefolder, args.folder, \
                                                      analysisdatabase) 
    
    sim = read_write.read_sim_info(datafolder)
    data = read_write.read_multid_analysis_data(analysisfile)

    ### plot the data in the given time window
    
    savebase = args.savebase + args.savefolder + "/"
    os.system("mkdir -p " + savebase)
    plot_vorticity(data, beads, cells, sim, args.init_time, args.fin_time, savebase, args.save_eps)
    
    return
    
##############################################################################

if __name__ == '__main__':
    main()    
    
##############################################################################
