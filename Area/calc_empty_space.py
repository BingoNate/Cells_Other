
""" Calculate empty spaces per simulation parameter"""

### example command line arguments: 
###    -fl=/local/duman/SIMULATIONS/Cells_in_LAMMPS/.../ 
###         -sb=/usr/users/iff_th2/duman/Cells_in_LAMMPS/DATA/
###             -sf=Empty_space

##############################################################################

import sys
sys.path.append('../Utility')

import argparse
import numpy as np
import os
import read_write
import misc_tools
    
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.collections import PolyCollection
import matplotlib.cm as cm
import matplotlib.colors as mplcolors

import matplotlib as mpl
mpl.use('Agg', warn=False)
import matplotlib.pyplot as plt
  
import PIL
from PIL import Image

##############################################################################

def pbc_dist_around_com(x1, x2, lx):
    """ correct the periodic boundary corrected distance between the positions x1-x2"""
    
    dx = x1-x2 
    
    if dx > lx/2.:
        dx = -lx
    elif dx < -lx/2.:
        dx = lx
    else:
        dx = 0.
        
    return dx

##############################################################################

def correct_pbc_around_com(xi, yi, xcom, ycom, nbeads, ncells, nbpc, lx, ly):
    """ correct bead positions according to cell center of mass positions"""
    
    xpbc = np.zeros((nbeads), dtype=np.float32)
    ypbc = np.zeros((nbeads), dtype=np.float32)

    
    k = 0
    for n in xrange(ncells):
        for j in xrange(nbpc[n]):
            dx = pbc_dist_around_com(xi[k], xcom[n], lx)
            dy = pbc_dist_around_com(yi[k], ycom[n], ly)
            xpbc[k] = xi[k] + dx
            ypbc[k] = yi[k] + dy
            k += 1
        
    return xpbc, ypbc   
 
##############################################################################
    
def plot_png(xpbc, ypbc, nbpc, ncells, lx, ly, frame):
    """ plot the frame as colored png image"""
 
    ax_len = 1.0                          # Length of one subplot square box
    ax_b = 0.0                            # Beginning/offset of the subplot in the box
    ax_sep = 0.0                          # Separation length between two subplots
    total_subplots_in_x = 1               # Total number of subplots    
    fig = plt.figure()
    subp = misc_tools.Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x) 
    ax0 = subp.addSubplot()
    
    patches = []
    k = 0
    for n in range(ncells):
        r = np.vstack((xpbc[k:k+nbpc[n]]/0.5, ypbc[k:k+nbpc[n]]/0.5)).T
        k += nbpc[n]
        polygon = Polygon(r, True)
        patches.append(polygon)
        
    p = PatchCollection(patches, facecolor='black', alpha=1.0)
    
    ax0.add_collection(p)
    ax0.autoscale_view()
    
    ax0.scatter(xpbc/0.5, ypbc/0.5, s=1.0, c='black', alpha=1.0)
    
    ax0.axis('equal')
    ax0.set_xlim((-10, lx/0.5+10))
    ax0.set_ylim((-10, ly/0.5+10))
    ax0.axis('off')
    
    savepath = "area_polygons_" + str(frame) + ".png"          
    plt.savefig(savepath, dpi=300, bbox_inches='tight', pad_inches=0.0)
    plt.clf()
    
    return savepath

##############################################################################

def convert_binary(pngpath):
    """ convert rgb png image to binary image and save it as numpy array"""

    im = Image.open(pngpath)
    bim = im.convert('1') 
    os.system('rm ' + pngpath)

    return np.ravel(np.asarray(bim))
    
##############################################################################
    
def calc_empty_space_ratio(bpos, cpos, nbpc, sim):
    """ calculate the ratio of empty space to the total space via image analysis"""
    
    empty_space_ratio = 0.
    for step in xrange(sim.nsteps):
        xpbc, ypbc = correct_pbc_around_com(bpos.xi[step, 0, :], bpos.xi[step, 1, :], 
                                            cpos.xi[step, 0, :], cpos[step, 1, :], 
                                            sim.nbeads, sim.ncells, nbpc, 
                                            sim.lx, sim.ly)
        pngpath = plot_png(xpbc, ypbc,
                           nbpc, sim.ncells, sim.lx, sim.ly, step)
        binary_data = convert_binary(pngpath)
        empty_space_ratio += float(len(binary_data[binary_data==False]))/float(len(binary_data))
        print "step / nsteps : ", str(step), " / ", str(sim.nsteps)
    empty_space_ratio /= sim.nsteps
        
    return empty_space_ratio

##############################################################################

def main():

    ### get the data folder
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-fl", "--folder", \
                        help="Folder containing data, as in /local/duman/SIMULATIONS/Cells_in_LAMMPS/.../")
    parser.add_argument("-sb", "--savebase", nargs="?", \
                        const = "/usr/users/iff_th2/duman/Cells_in_LAMMPS/DATA/", \
                        help="Folder to save the data, as in /usr/users/iff_th2/duman/Cells_in_LAMMPS/DATA/") 
    parser.add_argument("-sf", "--savefolder", nargs="?", \
                        const="Empty_space", \
                        help="Specific folder for saving, as in Empty_space")    
    parser.add_argument("-b","--bending", action="store_true", help="Decide whether to include bending or not")                
    args = parser.parse_args()
    
    ### read the data and general information from the folder
    
    sim, cells, beads = read_write.read_h5_file(args.folder, args.bending)
    cells.get_img_pos(sim.lx)
    beads.get_img_pos(sim.lx)
    print "folder = ", args.folder
        
    ### calculate the empty spaces
    
    empty_space_ratio = calc_empty_space_ratio(beads.xi, cells.xi, cells.nbpc, sim)
    print "folder = ", args.folder
    
    ### write the data to the corresponding file
   
    read_write.write_single_analysis_data(empty_space_ratio, args.savebase,
                                          args.savefolder, sim, args.bending)
    
    return
    
##############################################################################

if __name__ == '__main__':
    main()    
    
##############################################################################
