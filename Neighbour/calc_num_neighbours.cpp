
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


void build_linked_cell_list(const double * const *x, const double * const *y,
			    vector<int> & heads, vector<int> & llist, 
			    const double wbin, const int nboxes, const int nsize,
			    const int step, const int natoms, const double l) {
  /* build a linked cell/box list */
    
  for (int j = 0; j < natoms; j++) {
    int xbin = get_bin_number(get_single_img_pos(x[step][j], l), wbin, nboxes);
    int ybin = get_bin_number(get_single_img_pos(y[step][j], l), wbin, nboxes);
    int bin = xbin*nboxes + ybin;
    
    llist[j] = heads[bin];
    heads[bin] = j;
  }
      
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

double calc_num_neighbours (const double * const *x, const double * const *y, SimInfo sim) {
  /* calculate densities per bin per frame */
  
  double avg_num_neigh = 0.;
  
  // set linked cell/box list parameters
  
  double rcut = 2.*sim.sigma;
  int nboxes = int(floor(sim.lx/rcut));
  int nsize = nboxes*nboxes;
  rcut = sim.lx/nboxes;  
  vector<int> heads(nsize, -1);
  vector<int> llist(sim.nbeads*100, -1);

  vector<int> num_neighs_per_cell(sim.ncells, 0);

  for (int step = 0; step < sim.nsteps; step++) {
    cout << "step / nsteps : " << step << " / " << sim.nsteps << endl;
    
    // build the linked cell/box list
    
    build_linked_cell_list(x, y, heads, llist,
			  rcut, nboxes, nsize,
			  step, sim.nbeads, sim.lx);
			  
    // count the number of neighbours of each cell
    
    vector<int> num_neighs_per_cell_per_time(sim.ncells, 0);
    int k = 0;
    for (int n = 0; n < sim.ncells; n++) {
      for (int j = 0; j < sim.nbpc[j]; j++) {
        int xbin = get_bin_number(get_single_img_pos(x[step][k], sim.lx), rcut, nboxes);
        int ybin = get_bin_number(get_single_img_pos(y[step][k], sim.ly), rcut, nboxes);
        //int bin = xbin*nboxes + ybin;
        	
        // loop over the neighbouring bins
        
        for (int ix = -1; ix < 2; ix++) {
	  int xneighbin = wrap_to_range(ix+xbin, nboxes);
          
          for (int iy = -1; iy < 2; iy++) {
            int yneighbin = wrap_to_range(iy+ybin, nboxes);
            int neighbin = xneighbin*nboxes + yneighbin;
	    
            // restore the head bead
	    
	    int idx = heads[neighbin];
	    
	    // traverse the neighbouring bin
	    
	    while (idx != -1) {
	      
	      
	      // traverse the linked list 
	      
	      idx = llist[idx];
	      
	    }     // neighbouring bin
            
            
            
          }       // y neighbour bin loop
        }         // x neighbor bin loop
        k++;
      }           // nbpc loop
    }             // cell loop
    
  }               // timestep loop
  
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
  
  double avg_num_neigh = calc_num_neighbours(x, y, sim);
					    
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

