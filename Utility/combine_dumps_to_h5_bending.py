
""" Combine multiple dump files from different restart instances into a single hdf5 file"""

### example command line arguments: 
###    -fl=/homea/ias2/duman/Cells_in_LAMMPS/ -t -d=0.8 -e=1.0 -f=0.5 -a=10.0 -k=1.0
###         -dt=0.001 -ns=50000 -b=0.5 -s=1.0 -nc=5000

##############################################################################

import argparse
import numpy as np
import os
import h5py

##############################################################################

def read_contextual_info():
    """ read the contextual information provided by the user"""
    
    ### get the data folder and the last timestep info
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-fl", "--folder", help="Folder containing data")
    parser.add_argument("-t", "--tstep", nargs="?", const="100000000", \
                            type=int, help="The last time step that is being searched for")
    parser.add_argument("-d", "--density", type=float, help="Packing fraction of the system")  
    parser.add_argument("-e", "--eps", type=float, help="Strength of LJ interaction")
    parser.add_argument("-f", "--fp", type=float, help="Propulsion force")
    parser.add_argument("-a", "--areak", type=float, help="Area constraint strength")
    parser.add_argument("-k", "--kappa", type=float, help="Bending rigidity")
    parser.add_argument("-dt", "--timestep", type=float, help="Timestep of the simulation")
    parser.add_argument("-ns", "--nsamp", type=int, help="Sampling rate of data")
    parser.add_argument("-b", "--bl", type=float, help="Bond length of the simulation")    
    parser.add_argument("-s", "--sigma", type=float, help="Lennard Jones length")     
    parser.add_argument("-nc", "--ncells", type=int, help="Number of cells")         
    args = parser.parse_args()
    
    ### generate folder path
    
    folder = args.folder + 'density_' + str(args.density) + '/eps_' + str(args.eps) + \
        '/fp_' + str(args.fp) + '/areak_' + str(args.areak) + '/kappa_' + str(args.kappa)
    fpath = folder + '/out1.dump'
    assert os.path.exists(fpath), "\nOUT1.DUMP DOES NOT EXIST FOR: " + folder 

    print fpath

    ### determine the total number of beads and the box size
    
    fl = open(fpath, 'r')
    fl.readline()
    fl.readline()
    fl.readline()
    line = fl.readline()
    line = line.split()
    nbeads = int(line[0])
    
    fl.readline()
    line = fl.readline()
    line = line.split()
    lx = float(line[1])
    line = fl.readline()
    line = line.split()
    ly = float(line[1])

    fl.close()
    
    ### total number of steps
    
    nsteps = args.tstep/args.nsamp
    nsteps += 1
            
    return folder, nbeads, nsteps, lx, ly, args

##############################################################################
    
def get_number_of_snaps(f, nbeads):
    """ determine the number of snapshots in the file"""
    
    os.system('wc -l ' + f + ' > tmp.txt')
    ifile = open('tmp.txt')
    line = ifile.readline()
    line = line.split()
    nlines = int(line[0])
    nsnaps = nlines/(nbeads+9)
    ifile.close()
    os.system('rm tmp.txt')
    
    return nsnaps

##############################################################################
    
def read_pos(fl, x, nbeads, nsnaps, lx, ly, checked, tstep_cnt, T, d, mid):
    """ read the position data from the file and return the last timestep at the end"""

    ### read the positions unique per each tstep in a single dump file

    already_checked = False
    for snap in range(nsnaps):
        
        ### read the headers to check the uniqueness of current tstep
        
        fl.readline()
        line = fl.readline()
        line = line.split()
        tstep = int(line[0])
        
        ### finish if the last tstep is exceeded already
        
        if tstep > T:
            return tstep_cnt, tstep
        
        ### make sure the current tstep is unique
        
        if tstep not in checked:
            checked.append(tstep)
            tstep_cnt += 1
            already_checked = False
        else:
            print "ALREADY CHECKED " + str(tstep)             
            already_checked = True
    
        ### read the remaining part of the header part
        
        for j in range(7):
            fl.readline()
         
        ### read the positions per bead if the tstep is unique
        
        for j in range(nbeads):
            line = fl.readline()
            if already_checked:
                continue
            line = line.split()
            if len(line) < 2:
                print line
                continue
            bid = int(line[0]) - 1
            mid[bid] = int(line[1]) - 1
            x[tstep_cnt, 0, bid] = float(line[2])
            x[tstep_cnt, 1, bid] = float(line[3])
            d[tstep_cnt, mid[bid]] = float(line[4])

        print tstep_cnt, tstep

    return tstep_cnt, tstep

##############################################################################
    
def read_pos_from_dump_files(folder, nbeads, ncells, nsteps, T, lx, ly):
    """ read the position data of each dump file until the last tstep is reached"""

    ### generate file path and the total number of snapshots in the file
    
    current_file_number = 0
    x = np.zeros((nsteps, 2, nbeads), dtype=np.float32)
    d = np.zeros((nsteps, ncells), dtype=np.int32)    
    mid = np.zeros((nbeads), dtype=np.int32)
    tstep_cnt = -1
    checked = []
    tstep = 0
    
    while tstep < T:
        
        current_file_number += 1
        fpath = folder + '/out' + str(current_file_number) + '.dump'
        assert os.path.exists(fpath), "out dump file does NOT exist for: " + fpath
        fl = open(fpath, 'r')
        nsnaps = get_number_of_snaps(fpath, nbeads)
        
        ### read the positions unique per each tstep in a single dump file
        
        tstep_cnt, tstep = read_pos(fl, x, nbeads, nsnaps, lx, ly, checked, tstep_cnt, T, d, mid)
        fl.close()
        
    return x, d, mid
    
##############################################################################
    
def write_h5_file(folder, x, d, mid, com, nbeads, nsteps, nbpc, lx, ly, args):
    """ write data to hdf5 file"""
    
    ### file path
    
    fpath = folder + '/out.h5'
    fl = h5py.File(fpath, 'w')
    
    ### positions of beads
    
    bead = fl.create_group('beads')
    bead.create_dataset('xu', (nsteps, 2, nbeads), data=x, dtype=np.float32, compression='gzip') 
    bead.create_dataset('cid', data=mid) 
    
    
    ### cell information
    
    cell = fl.create_group('cells')
    cell.create_dataset('comu', (nsteps, 2, args.ncells), data=com, dtype=np.float32, compression='gzip') 
    cell.create_dataset('pol', (nsteps, args.ncells), data=d, dtype=np.float32, compression='gzip')
    cell.create_dataset('nbpc', data=nbpc)
    
    ### simulation information
    
    info = fl.create_group('info')
    box = info.create_group('box')
    box.create_dataset('x', data=lx)
    box.create_dataset('y', data=ly)
    info.create_dataset('dt', data=args.timestep)
    info.create_dataset('nsteps', data=nsteps)
    info.create_dataset('ncells', data=args.ncells)
    info.create_dataset('nbeads', data=nbeads)
    info.create_dataset('nsamp', data=args.nsamp)
    
    ### simulation parameters
    
    param = fl.create_group('param')
    param.create_dataset('eps', data=args.eps)
    param.create_dataset('rho', data=args.density)
    param.create_dataset('fp', data=args.fp)
    param.create_dataset('areak', data=args.areak)
    param.create_dataset('kappa', data=args.kappa)
    param.create_dataset('bl', data=args.bl)
    param.create_dataset('sigma', data=args.sigma)
    
    fl.close()
    
    return

##############################################################################

def get_beads_per_cell_data(folder, nbeads, args):
    """ append the number of beads per cell data to the existing data file"""
    
    nbpc = np.zeros((args.ncells), dtype=int)
    
    ### read the input file in equib until the bonds section
    
    ifilepath = args.folder + 'density_' + str(args.density) + '/equib/input.data'
    ifile = open(ifilepath, 'r')
    
    for j in range(2*nbeads):
        line = ifile.readline()
        line = line.split()
        if len(line) > 0:
            if line[0] == 'Bonds':
                break
    
    ### check the bonds to determine number of beads per cell
    
    k = 0
    ifile.readline()
    for j in range(nbeads):
        line = ifile.readline()
        line = line.split()
        b1 = int(line[2])
        b2 = int(line[3])
        if b1 > b2:
            nbpc[k] = b1-b2+1
            k += 1
                
    ifile.close()

    return nbpc
 
##############################################################################
    
def calculate_com_of_cells(xu, nsteps, nbpc, args):   
    """ calculate the center of mass of cells"""
    
    com = np.zeros((nsteps, 2, args.ncells), dtype=np.float32)
    
    k = 0
    for j in range(args.ncells):
        com[:, :, j] = np.mean(xu[:, :, k:k+nbpc[j]], axis=2)
        k += nbpc[j]
    
    return com

##############################################################################
    
def calculate_com_of_cells_one_liner(xu, mid, nsteps, nbeads, nbpc):   
    """ calculate the center of mass of cells"""
    
    splitted = np.split(xu, np.cumsum(nbpc)[:-1], axis=2)
    r = np.array([np.mean(sfil, axis=2) for sfil in splitted])
    com = np.swapaxes(np.swapaxes(r, 0, 1), 1, 2)
    
    return com
    
##############################################################################
    
def main():

    folder, nbeads, nsteps, lx, ly, args = read_contextual_info()
    xu, d, mid = read_pos_from_dump_files(folder, nbeads, args.ncells, nsteps, args.tstep, lx, ly)
    nbpc = get_beads_per_cell_data(folder, nbeads, args)
    com = calculate_com_of_cells(xu, nsteps, nbpc, args)
    write_h5_file(folder, xu, d, mid, com, nbeads, nsteps, nbpc, lx, ly, args)
    
    return
    
##############################################################################

if __name__ == '__main__':
    main()    
    
##############################################################################
