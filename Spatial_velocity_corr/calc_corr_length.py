
""" Calculate the correlation length from spatial velocity correlation 
    per simulation parameter"""

### example command line arguments: 
###    -e=1.0 -f=0.5 -a=1.0
###    -fl=/local/duman/SIMULATIONS/Cells_in_LAMMPS/density_0.8/ 
###         -sb=/usr/users/iff_th2/duman/Cells_in_LAMMPS/PLOTS/
###             -sf=Sp_velocity_corr

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

def calc_corr_len(data, sim):
    """ calculate the correlation length"""

    idx =  (np.abs(data-np.exp(-1))).argmin()
        
    return idx/sim.r_avg
        
##############################################################################

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--eps", type=float, \
                        help="Strength of LJ potential")
    parser.add_argument("-f", "--fp", type=float, \
                        help="Propulsion force")
    parser.add_argument("-a", "--areak", type=float, \
                        help="Area constraint potential strength")
    parser.add_argument("-k", "--kappa", type=float, \
                        help="Bending rigidity")    
    parser.add_argument("-fl", "--folder", nargs="?", \
                        const='/local/duman/SIMULATIONS/Cells_in_LAMMPS/density_0.8/', \
                        help="Folder containing data, as in /local/duman/SIMULATIONS/Cells_in_LAMMPS/density_0.8/")    
    parser.add_argument("-sb", "--savebase", nargs="?", \
                        const = "/usr/users/iff_th2/duman/Cells_in_LAMMPS/DATA/", \
                        help="Folder to save the data, as in /usr/users/iff_th2/duman/Cells_in_LAMMPS/DATA/") 
    parser.add_argument("-sf", "--savefolder", nargs="?", \
                        const="Corr_length", \
                        help="Specific folder for saving, as in Corr_length")    
    parser.add_argument("-s", "--save", action="store_true", 
                        help="Save for the phase diagram or not")
    args = parser.parse_args()
    
    ### read the data and general information from the folder
    
    print "Calculating correlation length for the following parameters : ", \
        args.eps, ", ", args.fp, ", ", args.areak, ", ", args.kappa
        
    ### get the simulation data
    
    folder = args.folder + "eps_" + str(args.eps) + \
        "/fp_" + str(args.fp) + "/areak_" + str(args.areak) + \
        "/kappa_" + str(args.kappa) + "/"
    if args.kappa == 100.0:
        sim = read_write.read_sim_info(folder, False)
    else:
        sim = read_write.read_sim_info(folder, True)
    
    ### calculate correlation length
    
    corr_len_folder = "/usr/users/iff_th2/duman/Cells_in_LAMMPS/DATA/"
    aname = "Sp_velocity_corr_subt"
    corr_len_folder += aname
    corr_len_file = corr_len_folder + "/" + aname + "_eps_" + str(args.eps) + \
        "_fp_" + str(args.fp) + "_areak_" + str(args.areak) + \
        "_kappa_" + str(args.kappa) + ".txt"
    x, y = read_write.read_2d_analysis_data(corr_len_file)
    corr_len = calc_corr_len(y, sim)
    
    ### write the data
    
    savebase = args.savebase + args.savefolder + "/"
    os.system("mkdir -p " + savebase)
    savefile = savebase + args.savefolder + "_eps_" + str(args.eps) + \
        "_fp_" + str(args.fp) + "_areak_" + str(args.areak) + \
        "_kappa_" + str(args.kappa) + ".txt"
    read_write.write_single_analysis_data(corr_len, savefile)
    
    ### collect the data for a phase diagram
    
    if args.save:
        aname = "PHASE_DIAGRAM"
        savebase = args.savebase + aname + "/"
        os.system("mkdir -p " + savebase)
        savefile = savebase + aname + ".txt"
        fl = open(savefile, "a")
        fl.write(str(sim.eps) + "\t" + str(sim.fp) + "\t" + str(sim.areak) \
            + "\t" + str(sim.kappa) + "\t" + str(corr_len) + "\n")
        fl.close()

    return
    
##############################################################################

if __name__ == '__main__':
    main()    
    
##############################################################################
