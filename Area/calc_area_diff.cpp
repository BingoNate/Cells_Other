
/* calculate the averaged area difference per cell and per simulation */

// COMPILATION AND RUN COMMANDS:
// g++ -O3 -Wl,-rpath=$HOME/hdf5/lib -L$HOME/hdf5/lib -I$HOME/hdf5/include ${spath}/calc_area_of_cells.cpp ${upath}/read_write.cpp -lhdf5 -lm -o calc_area_diff
// ./calc_area_diff out.h5 ${path}/Area_diff.txt

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

double calc_avg_area_diff (const double * const *x, const double * const *y,
                           vector<double> & init_areas,
                           const int nsteps, const int ncells, const int nbeads,
                           const vector<int> & nbpc) {
  /* calculate the area covered by cells per frame and return the averaged area over all the frames */
  
  double adiff = 0.;
  
  for (int step = 0; step < nsteps; step++) {
    
    cout << "steps / nsteps : " << step << " / " << nsteps << endl;
    
    double adiff_per_cell = 0.;
    int k = 0;
    
    for (int n = 0; n < ncells; n++) {
      double area = 0.;
      
      for (int j = 0; j < nbpc[n]-1; j++) {
        area += x[step][k]*y[step][k+1] - x[step][k+1]*y[step][k];
        k++;
      }		// beads of the cell
      
      area += x[step][k]*y[step][k-nbpc[n]+1] - x[step][k-nbpc[n]+1]*y[step][k];
      k++;
      area /= 2.;
      adiff_per_cell += area - init_areas[n];
      
    }	 	// cells
    
    adiff += adiff_per_cell/ncells;
    
  }     	// timesteps
  
  adiff /= nsteps;
  
  return adiff;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

vector<double> calc_initial_areas (const vector<int> & nbpc, const int ncells) {
  /* calculate the initial target areas of cells */
  
  vector<double> A0(ncells, 0.);
  for (int j = 0; j < ncells; j++) {
    double rad = nbpc[j]*0.5/(2.*pi);
    A0[j] = 0.9*pi*rad*rad;
  }
  
  return A0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

int main (int argc, char *argv[]) {

 // set data fetching properties and make allocations
 
 double **x;
 double **y;
 std::string cells_or_beads = "beads";	   // analysis will be done over cells or beads
 bool get_image = false;			   // analysis will be done over images
                                 // instead of unwrapped coords
					    
 // get the file name by parsing and load simulation data as well as positions
  
  string filename = argv[1];
  SimInfo sim(filename, cells_or_beads, get_image, x, y);  
  
  // print info
   
  cout << "Calculating area difference from " << cells_or_beads <<
    " for the following file: \n" << filename << endl;
  cout << "nsteps = " << sim.nsteps << endl;
  cout << "ncells = " << sim.ncells << endl;
  cout << "nbeads = " << sim.nbeads << endl;
  
  // perform the analysis

  vector<double> initial_areas = calc_initial_areas(sim.nbpc, sim.ncells);
  double avg_area_diff = calc_avg_area_diff(x, y, initial_areas,
					 sim.nsteps, sim.ncells, sim.nbeads, sim.nbpc);
					   
  // write the computed data
  
  string outfilepath = argv[2];
  cout << "Writing area difference of the cells to the following file: \n" << outfilepath << endl;
  write_single_analysis_data(avg_area_diff, outfilepath);
  
  // deallocate the arrays
  //
  // TODO: I should update this to unique pointers to avoid this deallocation business
  
  deallocate_2d_array(x, sim.nsteps);
  deallocate_2d_array(y, sim.nsteps);
  
  return 0;
}  
