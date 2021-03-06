
/* calculate the density histogram in terms of the center of mass of cells */

// COMPILATION AND RUN COMMANDS:
// g++ -O3 -Wl,-rpath=$HOME/hdf5/lib -L$HOME/hdf5/lib -I$HOME/hdf5/include ${spath}/calc_area_of_cells.cpp ${upath}/read_write.cpp -lhdf5 -lm -o calc_density_histo
// ./calc_density_histo out.h5 ${path}/Density_histo.txt

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

vector<double> calc_densities (const double * const *x, const double * const *y, SimInfo sim) {
  /* calculate densities per bin per frame */
  
  // set bin properties
  
  double avg_cell_diameter = 40.*0.5/(2.*pi);
  double bin_length = 4.*avg_cell_diameter;
  int nbins = int(floor(sim.lx/bin_length));
  cout << "number of bins : " << nbins << endl;
  int nsize = nbins*nbins;
  bin_length = sim.lx/nbins;
  vector<double> densities(nsize, 0.);
  
  // calculate the density per bin per frame
  
  for (int step = 0; step < sim.nsteps; step++) {
    cout << "step / nsteps : " << step << " / " << sim.nsteps << endl;
    
    for (int j = 0; j < sim.nbeads; j++) {
      int xbin = get_bin_number(x[step][j], bin_length, nbins);
      int ybin = get_bin_number(y[step][j], bin_length, nbins);
      int bin = xbin*nbins + ybin;
      densities[bin]++;
      
    }
  }
  
  // normalization
  
  for (int j = 0; j < nsize; j++) 
    densities[j] /= (sim.nsteps*bin_length*bin_length);
  
  return densities;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

int main (int argc, char *argv[]) {

 // set data fetching properties and make allocations
 
 double **x;
 double **y;
 std::string cells_or_beads = "beads";	   // analysis will be done over cells
 bool get_image = true;			   // analysis will be done over images 
					   // instead of unwrapped coords
					    
 // get the file name by parsing and load simulation data as well as positions
  
  string filename = argv[1];
  SimInfo sim(filename, cells_or_beads, get_image, x, y);  
  
  // print info
   
  cout << "Calculating density histogram for " << cells_or_beads << 
    " for the following file: \n" << filename << endl;
  cout << "nsteps = " << sim.nsteps << endl;
  cout << "ncells = " << sim.ncells << endl;
  cout << "nbeads = " << sim.nbeads << endl;
    
  // perform the analysis
  
  vector<double> densities = calc_densities(x, y, sim);
					    
  // write the computed data
  
  string outfilepath = argv[2];
  cout << "Writing density histogram to the following file: \n" << outfilepath << endl;  
  write_1d_vec_analysis_data(densities, densities.size(), outfilepath);
  
  // deallocate the arrays
  // 
  // TODO: I should update this to unique pointers to avoid this deallocation business
  
  deallocate_2d_array(x, sim.nsteps);
  deallocate_2d_array(y, sim.nsteps);
  
  return 0;
}  

//////////////////////////////////////////////////////////////////////////////////////////////////////////

