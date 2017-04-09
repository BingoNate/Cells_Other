
""" Generate restart jobs by checking the last timestep of LAMMPS simulations"""

##############################################################################

import argparse
import numpy as np
import os
import pandas as pd
from string import atof
import fnmatch

##############################################################################

def check_latest_dump(max_dump_number, last_dump_file, folder, last_dumps, finished_folders, tstep):
    
    os.system('grep ' + tstep + ' ' + last_dump_file + ' > tmp.txt')
    ifile = open('tmp.txt')
    line = ifile.readline()
    if len(line) > 0:
        finished_folders[folder] = 1
        ifile.close()
        os.system('rm tmp.txt')    
    else:
        ifile.close()
        os.system('rm tmp.txt') 
        max_dump_number = max_dump_number - 1
        if max_dump_number > 2:
            last_dump_file = folder + '/out' + str(max_dump_number) + '.dump'
            last_dumps[folder] = last_dump_file
            max_dump_number, last_dump_file = check_latest_dump(max_dump_number, last_dump_file, folder, last_dumps, finished_folders, tstep)
        else:
            return max_dump_number, last_dump_file
    
    return max_dump_number, last_dump_file
    
##############################################################################

def main(): 
    
    ### get the data folder
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--folder", help="Folder containing data")
    parser.add_argument("-t", "--tstep", nargs="?", const="100000000", \
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
                
                
    last_dumps = {}         # key->folder name, value->last dump folder  
    finished_folders = {}   # key->folder name, value->1(finished),0(unfinished)      
    cnt = 0
    
    for folder in folders:  
        
        if os.path.exists(folder):
            
            print 'Checking ', folder    
            
            ### index dump folders and get the largest dump file for each folder
            
            dump_numbers = []
            for fl in os.listdir(folder):
                if fnmatch.fnmatch(fl, '*.dump'):
                    dump_numbers.append(int(fl.split('.')[0][3:]))
            if len(dump_numbers) != 0:
                max_dump_number = max(dump_numbers)   
            else:
                continue
            last_dump_file = folder + '/out' + str(max_dump_number) + '.dump'
            last_dumps[folder] = last_dump_file
            
            ### check if the latest dump file containts the last time step,
            ### and mark the simulations as finished or not finished accordingly
            
            finished_folders[folder] = 0
            new_file_number = max_dump_number
            max_dump_number, last_dump_file = check_latest_dump(max_dump_number, last_dump_file, folder, last_dumps, finished_folders, args.tstep)
            if max_dump_number > 2:
                new_file_number = max_dump_number 
            print finished_folders[folder]

            ### if the simulation is not finished, create a restart job
            
            if finished_folders[folder] == 0:
                cnt += 1
                os.system('cp ' + folder + '/job' + str(new_file_number) + '.cmd ' + folder + '/job' + str(new_file_number+1) + '.cmd')
                os.system('cp ' + folder + '/input' + str(new_file_number) + '.lammps ' + folder + '/input' + str(new_file_number+1) + '.lammps')
                os.system('sed -i \"s/backup' + str(2*new_file_number-1) + '.rs/backup' + str(2*new_file_number+1) + '.rs/g\" ' + folder + '/input' + str(new_file_number+1) + '.lammps')
                os.system('sed -i \"s/backup' + str(2*new_file_number) + '.rs/backup' + str(2*new_file_number+2) + '.rs/g\" ' + folder + '/input' + str(new_file_number+1) + '.lammps')
                os.system('sed -i \"s/out' + str(new_file_number) + '.dump/out' + str(new_file_number+1) + '.dump/g\" ' + folder + '/input' + str(new_file_number+1) + '.lammps') 
                os.system('sed -i \"s/out' + str(new_file_number) + '.h5/out' + str(new_file_number+1) + '.h5/g\" ' + folder + '/input' + str(new_file_number+1) + '.lammps') 
                os.system('sed -i \"s/out' + str(new_file_number) + '.rs/out' + str(new_file_number+1) + '.rs/g\" ' + folder + '/input' + str(new_file_number+1) + '.lammps')       
                os.system('sed -i \"s/out' + str(new_file_number) + '.data/out' + str(new_file_number+1) + '.data/g\" ' + folder + '/input' + str(new_file_number+1) + '.lammps')             
                if new_file_number == 1:
                    os.system('sed -i \"s/read_data/read_restart/g\" ' + folder + '/input' + str(new_file_number+1) + '.lammps')
                    os.system('sed -i \"s/input.data/backup' + str(2*new_file_number) + '.rs/g\" ' + folder + '/input' + str(new_file_number+1) + '.lammps') 
                    os.system('sed -i \"s/100000000/50000000/g\" ' + folder + '/input' + str(new_file_number+1) + '.lammps')           
                else:
                    os.system('sed -i \"s/backup' + str(2*new_file_number-2) + '.rs' + '/backup' + str(2*new_file_number) + '.rs/g\" ' + folder + '/input' + str(new_file_number+1) + '.lammps')
                os.system('sed -i \"s/log' + str(new_file_number) + '.lammps/log' + str(new_file_number+1) + '.lammps/g\" ' + folder + '/job' + str(new_file_number+1) + '.cmd')             
                os.system('sed -i \"s/input' + str(new_file_number) + '.lammps/input' + str(new_file_number+1) + '.lammps/g\" ' + folder + '/job' + str(new_file_number+1) + '.cmd')   
                os.system('sed -i \"s/time=24:00:00/time=6:00:00/g\" ' + folder + '/job' + str(new_file_number+1) + '.cmd')   
                
                opath = 'joblist_restart/job' + str(cnt) + '.sh'
                ofile = open(opath, 'w')
                ofile.write('#!/bin/bash\n\n')
                ofile.write('cd ' + folder + '\n')
                ofile.write('chmod 755 job' + str(new_file_number+1) + '.cmd\n')
                ofile.write('sbatch job' + str(new_file_number+1) + '.cmd\n')
                ofile.write('cd -\n')
                ofile.close()
                
        else:
            print 'FOLDER DOES NOT EXIST, ', folder
            continue
            
##############################################################################        

if __name__ == '__main__':
    main()

##############################################################################

