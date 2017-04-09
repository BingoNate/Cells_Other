
""" Calculate the enstrophy
    per simulation parameter"""

### example command line arguments: 
###    -e=1.0 -f=0.5 -a=1.0
###    -fl=/local/duman/SIMULATIONS/Cells_in_LAMMPS/density_0.8/ 
###         -sb=/usr/users/iff_th2/duman/Cells_in_LAMMPS/PLOTS/
###             -sf=Enstophy

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

def calc_ens(data, sim):
    """ calculate the enstrophy"""

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
    
    print "Calculating enstrophy for the following parameters : ", \
        args.eps, ", ", args.fp, ", ", args.areak
        
    ### get the simulation data
    
    folder = args.folder + "eps_" + str(args.eps) + \
        "/fp_" + str(args.fp) + "/areak_" + str(args.areak) + "/"
    sim = read_write.read_sim_info(folder)
    
    ### calculate correlation length
    
    ens_folder = "/usr/users/iff_th2/duman/Cells_in_LAMMPS/DATA/"
    aname = "Enstrophy"
    ens_folder += aname
    ens_file = ens_folder + "/" + aname + "_eps_" + str(args.eps) + \
        "_fp_" + str(args.fp) + "_areak_" + str(args.areak) + ".txt"
    ens = read_write.read_single_analysis_data(ens_file)
    
    ### write the data
    
    savebase = args.savebase + args.savefolder + "/"
    os.system("mkdir -p " + savebase)
    savefile = savebase + args.savefolder + "_eps_" + str(args.eps) + \
        "_fp_" + str(args.fp) + "_areak_" + str(args.areak) + ".txt"
    read_write.write_single_analysis_data(ens, savefile)
    
    ### collect the data for a phase diagram
    
    if args.save:
        aname = "PHASE_DIAGRAM"
        savebase = args.savebase + aname + "/"
        os.system("mkdir -p " + savebase)
        savefile = savebase + aname + ".txt"
        fl = open(savefile, "a")
        fl.write(str(sim.eps) + "\t" + str(sim.fp) + "\t" + str(sim.areak) \
            + "\t" + str(ens) + "\n")
        fl.close()

    return
    
##############################################################################

if __name__ == '__main__':
    main()    
    
##############################################################################
