
/* calculate the area of the cells */

// COMPILATION AND RUN COMMANDS:
// g++ -O3 -Wl,-rpath=$HOME/hdf5/lib -L$HOME/hdf5/lib -I$HOME/hdf5/include ${spath}/calc_area_of_cells.cpp ${upath}/read_write.cpp -lhdf5 -lm -o calc_area_of_cells
// ./calc_area_of_cells out.h5 ${path}/Area_of_cells.txt

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

double calc_areas_per_frame (const double * const *x, const double * const *y,
			     vector<double> & area_per_frame,
			     const int nsteps, const int ncells, const int nbeads,
			     const vector<int> & nbpc) {
  /* calculate the area covered by cells per frame and return the averaged area over all the frames */
  
  double avg_area = 0.;
  
  for (int step = 0; step < nsteps; step++) {
    
    cout << "steps / nsteps : " << step << " / " << nsteps << endl;
    
    double areapf = 0.;
    int k = 0;
    
    for (int n = 0; n < ncells; n++) {
      double areapc = 0.;
      
      for (int j = 0; j < nbpc[n]-1; j++) {
        areapc += x[step][k]*y[step][k+1] - x[step][k+1]*y[step][k];
        k++;
	
      }		// nbpc
      
      areapc += x[step][k]*y[step][k-nbpc[n]+1] - x[step][k-nbpc[n]+1]*y[step][k];
      k++;
      areapf += areapc/2.;
      
    }	 	// ncells
    
    areapf /= ncells;
    area_per_frame[step] = areapf; 
    avg_area += areapf;
    
  }     	// nsteps
  
  avg_area /= nsteps;
  
  return avg_area;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

int main (int argc, char *argv[]) {

 // set data fetching properties and make allocations
 
 double **x;
 double **y;
 std::string cells_or_beads = "beads";	   // analysis will be done over cells
 bool get_image = false;			   // analysis will be done over images
					   // instead of unwrapped coords
					    
 // get the file name by parsing and load simulation data as well as positions
  
  string filename = argv[1];
  SimInfo sim(filename, cells_or_beads, get_image, x, y);  
  
  // print info
   
  cout << "Calculating areas from " << cells_or_beads << 
    " for the following file: \n" << filename << endl;
  cout << "nsteps = " << sim.nsteps << endl;
  cout << "ncells = " << sim.ncells << endl;
  cout << "nbeads = " << sim.nbeads << endl;

  // set variables related to the analysis
  
  vector<double> area_per_frame(sim.nsteps, 0.);
  
  // perform the analysis
  
  double avg_area = calc_areas_per_frame(x, y, area_per_frame, 
					 sim.nsteps, sim.ncells, sim.nbeads, sim.nbpc);
					   
  // write the computed data
  
  string outfilepath = argv[2];
  cout << "Writing average area of the cells to the following file: \n" << outfilepath << endl;  
  write_single_analysis_data(avg_area, outfilepath);
  string outfilepath_2 = argv[3];
  cout << "Writing area of the cells per frame to the following file: \n" << outfilepath_2 << endl;    
  write_1d_vec_analysis_data(area_per_frame, sim.nsteps, outfilepath_2);
  
  // deallocate the arrays
  
  for (int i = 0; i < sim.nsteps; i++) {
    delete [] x[i];
    delete [] y[i];
  }
  delete [] x;
  delete [] y;
  
  return 0;
}  
