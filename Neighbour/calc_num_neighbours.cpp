
/* calculate average number of neighbours per cell per simulation parameters */

// COMPILATION AND RUN COMMANDS:
// g++ -O3 -Wl,-rpath=$HOME/hdf5/lib -L$HOME/hdf5/lib -I$HOME/hdf5/include ${spath}/calc_num_neighbours.cpp ${upath}/read_write.cpp -lhdf5 -lm -o calc_num_neighbours
// ./calc_num_neighbours out.h5 ${path}/Num_neighbour.txt

//////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <map>
#include <set>
#include <stack>
#include <algorithm>
#include "../Utility/read_write.hpp"
#include "../Utility/sim_info.hpp"

#define pi M_PI

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////

double calc_num_neighbours (const double * const *x, const double * const *y, SimInfo sim) {
  /* calculate densities per bin per frame */
  
  double avg_num_neigh = 0.;
  
  // build a neighbour list 
  
  for (int step = 0; step < sim.nsteps; step++) {
    
  }
  
  return avg_num_neigh;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

int main (int argc, char *argv[]) {

 // set data fetching properties and make allocations
 
 double **x;
 double **y;
 std::string cells_or_beads = "beads";	   // analysis will be done over cells
 bool get_image = false;		   // analysis will be done over images 
					   // instead of unwrapped coords
					    
 // get the file name by parsing and load simulation data as well as positions
  
  string filename = argv[1];
  SimInfo sim(filename, cells_or_beads, get_image, x, y);  
  
  // print info
   
  cout << "Calculating number of neighbours from " << cells_or_beads << 
    " for the following file: \n" << filename << endl;
  cout << "nsteps = " << sim.nsteps << endl;
  cout << "ncells = " << sim.ncells << endl;
  cout << "nbeads = " << sim.nbeads << endl;
    
  // perform the analysis
  
  
					    
  // write the computed data
  
  string outfilepath = argv[2];
  //cout << "Writing number of neighbours to the following file: \n" << outfilepath << endl;  
  //write_1d_vec_analysis_data(densities, densities.size(), outfilepath);
  
  // deallocate the arrays
  // 
  // TODO: I should update this to unique pointers to avoid this deallocation business
  
  deallocate_2d_array(x, sim.nsteps);
  deallocate_2d_array(y, sim.nsteps);
  
  return 0;
}  

//////////////////////////////////////////////////////////////////////////////////////////////////////////

