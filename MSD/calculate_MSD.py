
""" Calculate the MSD of the centre of mass of cells per simulation parameter"""

### example command line arguments: 
###    -fl=/local/duman/SIMULATIONS/Cells_in_LAMMPS/.../ 
###         -sb=/usr/users/iff_th2/duman/Cells_in_LAMMPS/DATA/
###             -sf=MSD

##############################################################################

import sys
sys.path.append('../Utility')

import argparse
import numpy as np
import os
import read_write
        
##############################################################################

def calculate_MSD(xu, sim):
    """ calculate and average the MSD (time and ensemble avgs)"""
        
    ndelay = int(sim.nsteps/2)
    delay = np.zeros((ndelay), dtype=np.int64)
    msd = np.zeros((ndelay), dtype=np.float64)
    
    for d in range(1, ndelay):
        delay[d] = d*sim.dt
        msd[d] = np.mean(np.mean( \
            (xu[d:,0,:]-xu[:-d,0,:])**2 + \
                (xu[d:,1,:]-xu[:-d,1,:])**2, axis=1) )
        
    return delay, msd        
   

##############################################################################

def calculate_drift_subtracted_MSD(xu, sim):
    """ calculate and average the drift subtracted MSD (time and ensemble avgs)"""
        
    ndelay = int(sim.nsteps/2)
    delay = np.zeros((ndelay), dtype=np.int64)
    msd = np.zeros((ndelay), dtype=np.float64)
    
    for d in range(1, ndelay):
        delay[d] = d*sim.dt
        xd = xu[d:,0,:]-xu[:-d,0,:]
        xd_mean = np.mean(xd, axis=1)
        yd = xu[d:,1,:]-xu[:-d,1,:]
        yd_mean = np.mean(yd, axis=1)
        xd -= xd_mean[:,None]
        yd -= yd_mean[:, None]
        msd[d] = np.mean(np.mean( \
            xd**2 + yd**2, axis=1) )
        
    return delay, msd  
     
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
                        const="MSD", \
                        help="Specific folder for saving, as in MSD")    
    parser.add_argument("-b","--bending", action="store_true", help="Decide whether to include bending or not")                
    parser.add_argument("-s","--save_eps", action="store_true", help="Decide whether to save in eps or not")            
    args = parser.parse_args()
    
    ### read the data and general information from the folder
    
    sim, cells, beads = read_write.read_h5_file(args.folder, args.bending)
    print "folder = ", args.folder, " bending rigidity = ", sim.kappa
        
    ### calculate MSD of the centre of mass of cells
    
    delay, msd = calculate_drift_subtracted_MSD(cells.xu, sim)
    print "folder = ", args.folder
    
    ### write the MSD data to the corresponding file
    
    read_write.write_2d_analysis_data(delay, msd, args.savebase, args.savefolder, sim, args.bending)
    
    return
    
##############################################################################

if __name__ == '__main__':
    main()    
    
##############################################################################
