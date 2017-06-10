
""" Calculate the 4 point susceptibility of the centre of mass of cells per parameter"""

### example command line arguments: 
###    -fl=/local/duman/SIMULATIONS/Cells_in_LAMMPS/.../ 
###         -sb=/usr/users/iff_th2/duman/Cells_in_LAMMPS/DATA/
###             -sf=4_pt_susceptibility

##############################################################################

import sys
sys.path.append('../Utility')

import argparse
import numpy as np
import numpy.ma as ma
import read_write

##############################################################################
     
def calc_4_pt_suscp(xu, sim, threshold_amplitude):
    """ calculate and average the susceptibility"""
        
    ndelay = int(sim.nsteps/2)
    delay = np.zeros((ndelay), dtype=np.int64)
    overlap = np.zeros((ndelay), dtype=np.float64)
    overlap_sq = np.zeros((ndelay), dtype=np.float64)  
    suscp = np.zeros((ndelay), dtype=np.float64)  
    overlap[0] = 1.0
    
    for d in range(1, ndelay):
        delay[d] = d*sim.dt

        ### subtract the mean displacement from cell displacements at the given delay time
        
        xd = xu[d:,0,:]-xu[:-d,0,:]
        xd_mean = np.mean(xd, axis=1)
        yd = xu[d:,1,:]-xu[:-d,1,:]
        yd_mean = np.mean(yd, axis=1)
        xd -= xd_mean[:,None]
        yd -= yd_mean[:,None]

        ### calculate the displacement magnitudes per cell per time
        
        displacements = np.sqrt(xd**2 + yd**2)
        
        for t0 in range(sim.nsteps-ndelay):
        
            ### calculate the ensemble averaged overlap function
            
            masked = ma.masked_greater(displacements[t0], threshold_amplitude)  # filter large displacements
            masked = masked[~masked.mask]                                       # get the low displacements
            low_displacements = float(len(masked.data))                         # get the number of low displacements
            overlap_ens_avg = low_displacements / sim.ncells                    # ens. averaged overlap fnc.
            overlap[d] += overlap_ens_avg
            overlap_sq[d] += overlap_ens_avg**2

        overlap[d] /= (sim.nsteps-ndelay)
        overlap_sq[d] /= (sim.nsteps-ndelay)
        suscp[d] = sim.ncells*(overlap_sq[d] - overlap[d]**2)

    return delay, suscp   
        
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
    parser.add_argument("-b","--bending", action="store_true", help="Decide whether to include bending or not")                    
    parser.add_argument("-s","--save_eps", action="store_true", help="Decide whether to save in eps or not") 
    args = parser.parse_args()
    
    ### read the data and general information from the folder
    
    sim, cells, beads = read_write.read_h5_file(args.folder, args.bending)
    print "folder = ", args.folder, " bending rigidity = ", sim.kappa
        
    ### calculate the 4 point susceptibility of displacements of the centre of mass of cells

    threshold_amplitude = 4.0
    delay, suscp = calc_4_pt_suscp(cells.xu, sim, threshold_amplitude)
    
    ### write the 4 point susceptibility data to the corresponding file
    
    read_write.write_2d_analysis_data(delay, suscp, args.savebase, args.savefolder, sim, args.bending)
    
    return
    
##############################################################################

if __name__ == '__main__':
    main()    
    
##############################################################################
