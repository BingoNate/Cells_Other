
""" Calculate the overlap function of displacements of the centre of mass of cells per parameter"""

### example command line arguments: 
###    -fl=/local/duman/SIMULATIONS/Cells_in_LAMMPS/.../ 
###         -sb=/usr/users/iff_th2/duman/Cells_in_LAMMPS/DATA/
###             -sf=Overlap_fnc

##############################################################################

import sys
sys.path.append('../Utility')

import argparse
import numpy as np
import numpy.ma as ma
import read_write

##############################################################################
     
def calculate_overlap_fnc(xu, sim, threshold_amplitude):
    """ calculate and average the overlap fnc (time and ensemble avgs)"""
        
    ndelay = int(sim.nsteps/2)
    delay = np.zeros((ndelay), dtype=np.int64)
    overlap = np.zeros((ndelay), dtype=np.float64)
    overlap[0] = 1.0
    
    for d in range(1, ndelay):
        delay[d] = d*sim.dt

        ### subtract the mean displacement at the given delay time
        
        xd = xu[d:,0,:]-xu[:-d,0,:]
        xd_mean = np.mean(xd, axis=1)
        yd = xu[d:,1,:]-xu[:-d,1,:]
        yd_mean = np.mean(yd, axis=1)
        xd -= xd_mean[:,None]
        yd -= yd_mean[:, None]

        ### calculate the displacements
        
        displacements = np.sqrt(xd**2 + yd**2)
        
        ### calculate the overlap function
        
        masked = ma.masked_greater(displacements, threshold_amplitude)
        masked = masked[~masked.mask]
        low_displacements = float(len(masked.data))
        total_displacements = float(np.size(displacements))
        overlap[d] = low_displacements / total_displacements

    return delay, overlap        
        
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
                        const='Overlap_fnc', \
                        help="Specific folder for saving, as in Overlap_fnc")     
    parser.add_argument("-s","--save_eps", action="store_true", help="Decide whether to save in eps or not") 
    args = parser.parse_args()
    
    ### read the data and general information from the folder
    
    sim, cells, beads = read_write.read_h5_file(args.folder)
    print "folder = ", args.folder, "cells.xu = \n", cells.xu
        
    ### calculate the overlap fnc of displacements of the centre of mass of cells

    threshold_amplitude = 4.0
    delay, overlap = calculate_overlap_fnc(cells.xu, sim, threshold_amplitude)
    
    ### write the overlap fnc data to the corresponding file
    
    read_write.write_2d_analysis_data(delay, overlap, args.savebase, args.savefolder, sim)
    
    return
    
##############################################################################

if __name__ == '__main__':
    main()    
    
##############################################################################
