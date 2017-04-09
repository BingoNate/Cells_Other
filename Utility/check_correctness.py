
""" Check the correctness of the dump files and report errors"""

### example command line arguments: 
###    -fl=/homea/ias2/duman/Cells_in_LAMMPS/ -t -d=0.8 -e=1.0 -f=0.5 -a=10.0 
###         -dt=0.001 -ns=50000 -b=0.5 -s=1.0 -nc=5000

##############################################################################

import argparse
import numpy as np
import os
import pandas as pd
from string import atof
import fnmatch

##############################################################################

def split_line(fl, read_line_cnt):
    """ read the line and split it into pieces"""
    
    line = fl.readline()
    read_line_cnt += 1
    line = line.split()
    
    return line, read_line_cnt

##############################################################################

def report_errors(flpath, nsamp, natoms):
    """ find the errors in the dump file if any 
        by checking sampling frequency of data and number of atoms"""
        
    print 'Checking ', flpath
    
    ### get the number of lines in the dump file
    
    os.system('wc -l ' + flpath + ' > tmp.txt')
    ifile = open('tmp.txt')
    line = ifile.readline()
    line = line.split()
    nlines = int(line[0])
    ifile.close()
    os.system('rm tmp.txt')
    
    ### read the file line by line
    
    fl = open(flpath, 'r')
    read_line_cnt = 0
    line, read_line_cnt = split_line(fl, read_line_cnt)
    line, read_line_cnt = split_line(fl, read_line_cnt)
    ti = int(line[0])
    atom_cnt = -4
    
    while read_line_cnt < nlines:
                
        line, read_line_cnt = split_line(fl, read_line_cnt)
                
        if len(line) > 0:
            
            if line[0] == 'ITEM:':
                                
                if line[1] == 'TIMESTEP':
                    if atom_cnt != natoms:
                        print "THERE IS A PROBLEM with atom count!", atom_cnt, natoms
                        
                    atom_cnt = 0
                    line, read_line_cnt = split_line(fl, read_line_cnt)
                    tj = int(line[0])
                    tdiff = tj-ti
                    if tdiff != nsamp:
                        print "THERE IS A PROBLEM with sampling!, ti = ", ti, " tj = ", tj
                    ti = tj
                    for i in range(7):
                        line, read_line_cnt = split_line(fl, read_line_cnt)
                        
            else:
                atom_cnt += 1
                    
        else:
            print "PROBLEM WITH length of line: ", line     
            
    
    fl.close()
        
    return
          
##############################################################################
    
def main():

    ### get the data folder
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--folder", nargs="?", \
                        const="/homea/ias2/duman/Cells_in_LAMMPS/density_0.8/", \
                        help="Folder containing data")
    parser.add_argument("-t", "--tstep", nargs="?", const="100000000", \
                            help="The last time step that is being searched for")
    parser.add_argument("-ns", "--nsamp", nargs="?", const="50000", \
                            help="The last time step that is being searched for")
    parser.add_argument("-n", "--natoms", nargs="?", const="197277", \
                            help="The last time step that is being searched for")    
    args = parser.parse_args()

    ### index the data
    
    eps = [0.05, 0.5, 1.0, 10.0, 20.0]
    fp = [0.0, 0.5, 1.0, 3.0, 5.0, 10.0]
    areak = [1.0, 10.0, 100.0]
    
    ### generate folder paths
    
    folders = []
    for e in eps:
        for f in fp:
            for a in areak:
                folders.append( args.folder + 'eps_' + str(e) + \
                               '/fp_' + str(f) + '/areak_' + str(a) )
                
                    
    ### run over all the folders containing set of dump files
    
    for folder in folders:  
        
        if os.path.exists(folder):
            
            print 'INSIDE ', folder    
            
            ### get the largest dump file for each folder
            
            dump_numbers = []
            for fl in os.listdir(folder):
                if fnmatch.fnmatch(fl, '*.dump'):
                    dump_numbers.append(int(fl.split('.')[0][3:]))
            if len(dump_numbers) != 0:
                max_dump_number = max(dump_numbers)   
            else:
                max_dump_number = 1
                
            ### run over each dump file

            for j in range(1, max_dump_number+1):
                
                flpath = folder + '/out' + str(j) + '.dump'
                report_errors(flpath, int(args.nsamp), int(args.natoms))
    
    return
    
##############################################################################

if __name__ == '__main__':
    main()    
    
##############################################################################