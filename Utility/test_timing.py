
""" Test platform for profiling"""

import numpy as np
import timeit

##############################################################################

def wrapper(func, *args, **kwargs):
    def wrapped():
        return func(*args, **kwargs)
    return wrapped
    
##############################################################################
    
def calculate_com_of_cells(xu, nsteps, nbpc, ncells):   
    """ calculate the center of mass of cells"""
    
    com = np.zeros((nsteps, 2, ncells), dtype=np.float32)
    
    k = 0
    for j in range(ncells):
        com[:, :, j] = np.mean(xu[:, :, k:k+nbpc[j]], axis=2)
        k += nbpc[j]
    
    return com

##############################################################################
    
def calculate_com_of_cells_one_liner(xu, nsteps, nbpc):   
    """ calculate the center of mass of cells"""
    
    splitted = np.split(xu, np.cumsum(nbpc)[:-1], axis=2)
    r = np.array([np.mean(sfil, axis=2) for sfil in splitted])
    com = np.swapaxes(np.swapaxes(r, 0, 1), 1, 2)
    
    return com
    
##############################################################################
    
def main():

    base = np.linspace(1, 1000, num=1000)
    nsteps = 50
    nbeads = 10
    xu = np.reshape(base, (nsteps, 2, nbeads))
    nbpc = [2, 3, 4, 1]
    ncells = len(nbpc)

    com1 = calculate_com_of_cells(xu, nsteps, nbpc, ncells)
    com2 = calculate_com_of_cells_one_liner(xu, nsteps, nbpc)
    
#    print 'ARRAYS ARE: \n'
#    print com1
#    print com2

    print 'EQUALITY OF THE ARRAYS: ', np.array_equal(com1, com2), '\n\n'

    wrapped = wrapper(calculate_com_of_cells, xu, nsteps, nbpc, ncells)    
    time1 = timeit.timeit(wrapped, number=100)
    print 'com: ', time1

    wrapped = wrapper(calculate_com_of_cells_one_liner, xu, nsteps, nbpc)   
    time2 = timeit.timeit(wrapped, number=100)
    print 'com_one_liner: ', time2
    
    print 'speedup is: ', max(time1, time2)/min(time1, time2)
    
    return
    
##############################################################################

if __name__ == '__main__':
    main()    
    
##############################################################################