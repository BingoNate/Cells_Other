
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
#include <tuple>
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
    
  cout << "Building linked cell list" << endl;
  
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

set<int> get_neighs_of_a_cell (const int start, const int end,
			       const double * const *x, const double * const *y,
			       const int step, const double rcut, const int nboxes,
			       const double lx, const double ly,
			       const vector<int> & cid,
			       const vector<int> & heads, const vector<int> & llist) {
  /* get the neighbours of a single cell at a given timestep */

  // set variables related to the analysis
  
  set<int> neighs;
  const double rth2 = rcut*rcut;
  
  // run over each bead of the cell 
  
  for (int j = start; j < end; j++) {
    int host_cell_id = cid[j];
    
    int xbin = get_bin_number(get_single_img_pos(x[step][j], lx), rcut, nboxes);
    int ybin = get_bin_number(get_single_img_pos(y[step][j], ly), rcut, nboxes);
	    
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
	  
	  // calculate the distance between the beads if they do not belong to the same cell 

	  int neigh_cell_id = cid[idx];
	  if (host_cell_id != neigh_cell_id) {
	    
	    double dx = x[step][idx] - x[step][j];
	    dx = get_min_img_dist(dx, lx);
	    double dy = y[step][idx] - y[step][j];
	    dy = get_min_img_dist(dy, ly);
	    
	    if (dx*dx + dy*dy < rth2) 
	      neighs.insert(neigh_cell_id);
	    
	  }
	  
	  // traverse the linked list 
	  
	  idx = llist[idx];
	  
	}     // inside neighbouring bin linked list
	
      }       // y neighbour bins loop
    }         // x neighbour bins loop
  }           // bead loop
    
  return neighs;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

map<int, set<int> > build_neigh_list (const double * const *x, const double * const *y,
				      const vector<int> & cid, const vector<int> & nbpc,
				      const vector<int> & heads, const vector<int> & llist,
				      const int ncells, const double rcut, const int nsteps,
				      const int step, const int nboxes, 
				      const double lx, const double ly) {
  /* build the neighbour list of all cells at a given timestep */

  cout << "Building neighbour list" << endl;
  
  map<int, set<int> > neigh_list;
  
  int bead_idx = 0;
  for (int n = 0; n < ncells; n++) {
    
    // iterate over the beads of the cell and get the neighbours of the cell
    
    int last_bead_idx = bead_idx + nbpc[n];
    set<int> neighs_of_this_cell = get_neighs_of_a_cell(bead_idx, last_bead_idx, x, y, 
							step, rcut, nboxes, lx, ly, 
							cid, heads, llist);
    bead_idx = last_bead_idx;
    
    // add the neighbours of this cell to the neighbour list
    
    neigh_list.insert(pair <int, set<int> > (n, neighs_of_this_cell));
    
  }             // cell loop
  
  
  return neigh_list;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

tuple<double, double> calc_num_neighbours (const double * const *x, const double * const *y, SimInfo sim) {
  /* calculate densities per bin per frame */
  
  // set parameters related to the analysis 
  
  double avg_num_neigh = 0.;
  double avg_num_surface_cells = 0.;
  
  // set linked cell/box list parameters
  
  double rcut = 2.*sim.sigma;
  int nboxes = int(floor(sim.lx/rcut));
  int nsize = nboxes*nboxes;
  rcut = sim.lx/nboxes;  
  vector<int> heads(nsize, -1);
  vector<int> llist(sim.nbeads, -1);

  for (int step = 0; step < sim.nsteps; step++) {
    cout << "step / nsteps : " << step << " / " << sim.nsteps << endl;
    
    // build the linked cell/box list

    vector<int> heads(nsize, -1);
    vector<int> llist(sim.nbeads, -1);
    build_linked_cell_list(x, y, heads, llist,
			   rcut, nboxes, nsize,
			   step, sim.nbeads, sim.lx);
			  
    // build the neighbour list
    
    map<int, set<int> > neigh_list = build_neigh_list(x, y, sim.cid, sim.nbpc,
						      heads, llist, 
						      sim.ncells, rcut, sim.nsteps, step, 
						      nboxes, sim.lx, sim.ly);
				  
     // count the average number of neighbours of cells
     
     double avg_num_neigh_per_cell = 0.;
     for (const auto& kv : neigh_list) {
       int num_neighs = kv.second.size();
       avg_num_neigh_per_cell += num_neighs;  
       if (num_neighs < 4)
	 avg_num_surface_cells++;
     }
     avg_num_neigh_per_cell /= neigh_list.size();
     avg_num_neigh += avg_num_neigh_per_cell;
     avg_num_surface_cells /= neigh_list.size();
    
  }               // timestep loop
  
  avg_num_neigh /= sim.nsteps;
  avg_num_surface_cells /= sim.nsteps;
  
  
  return make_tuple(avg_num_neigh, avg_num_surface_cells);
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
  
  double avg_num_neigh, avg_num_surface_cells;
  tie(avg_num_neigh, avg_num_surface_cells) = calc_num_neighbours(x, y, sim);
					    
  // write the computed data
  
  string outfilepath = argv[2];
  cout << "Writing number of neighbours to the following file: \n" << outfilepath << endl;  
  write_single_analysis_data(avg_num_neigh, outfilepath);
  
  string outfilepath_2 = argv[3];
  cout << "Writing number of surface cells to the following file: \n" << outfilepath_2 << endl;  
  write_single_analysis_data(avg_num_surface_cells, outfilepath_2);
  
  // deallocate the arrays
  // 
  // TODO: I should update this to unique pointers to avoid this deallocation business
  
  deallocate_2d_array(x, sim.nsteps);
  deallocate_2d_array(y, sim.nsteps);
  
  return 0;
}  

//////////////////////////////////////////////////////////////////////////////////////////////////////////

