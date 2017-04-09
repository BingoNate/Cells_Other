
/* calculate polar order parameter from the velocity vector */

// COMPILATION AND RUN COMMANDS:
// g++ -O3 -lm -Wl,-rpath=$HOME/hdf5/lib -L$HOME/hdf5/lib -I$HOME/hdf5/include ${spath}/calc_polar_order_intrinsic.cpp ${upath}/read_write.cpp -lhdf5 -o calc_polar_order_intrinsic
// ./calc_polar_order_intrinsic out.h5 ${path}/Polar_order_intrinsic.txt

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
#include "../Utility/basic.hpp"

#define pi M_PI

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////

double calc_polar_order_param (double **pol, 
			       const int nsteps, const int ncells) {
  /* calculate the polar order parameter */
    
  double polar_order_param = 0.;
  
  for (int step = 0; step < nsteps; step++) {
    
    double cost = 0.;
    double sint = 0.;
    
    for (int j = 0; j < ncells; j++) {
      
      double angle = pol[step][j]*pi/180;
      cost += cos(angle);
      sint += sin(angle);
      
    }	// cell loop
    
    polar_order_param += sqrt(cost*cost + sint*sint);
    
  }  	// timestep loop
  
  polar_order_param /= (nsteps*ncells);
    
  return polar_order_param;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////

int main (int argc, char *argv[]) {

  // get the file name by parsing
  
  string filename = argv[1];
  cout << "Calculating polar order parameter of the following file: \n" << filename << endl;

  // read in general simulation data
  
  int nsteps, nbeads, nsamp, ncells;
  nsteps = nbeads = nsamp = ncells = 0;
  double lx, ly, dt, eps, rho, fp, areak, bl, sigma;
  lx = ly = dt = eps = rho = fp = areak = bl = sigma = 0.;
  
  read_sim_data(filename, nsteps, nbeads, nsamp, ncells, lx, ly,
                dt, eps, rho, fp, areak, bl, sigma);
  
  // print simulation information
  
  cout << "nsteps = " << nsteps << endl;
  cout << "ncells = " << ncells << endl;
  
  // read in array data
  
  int nbpc[ncells];
  read_integer_array(filename, "/cells/nbpc", nbpc);
 
  // read in the cell position data all at once
  
  /* the position data is stored in the following format:
  (nsteps, 2, ncells)
  the data will be loaded as follows:
  (nsteps, ncells) in x and y separately 
  */
  
  double **x = new double*[nsteps];
  for (int i = 0; i < nsteps; i++) x[i] = new double[ncells];
  
  double **y = new double*[nsteps];
  for (int i = 0; i < nsteps; i++) y[i] = new double[ncells];
  
  double **pol = new double*[nsteps];
  for (int i = 0; i < nsteps; i++) pol[i] = new double[ncells];
  
  for (int i = 0; i < nsteps; i++) {
    for (int j = 0; j < ncells; j++) {
      x[i][j] = 0.;  y[i][j] = 0.;  pol[i][j] = 0.;
    }
  }
  
  read_all_pos_data(filename, x, y, nsteps, ncells, "/cells/comu");

  // read polarity data
  
  read_all_pol_data(filename, pol, nsteps, ncells, "/cells/pol"); 

  // calculate the velocity correlation in space
  
  double polar_order_param = calc_polar_order_param(pol, nsteps, ncells);

  // write the computed data
  
  string outfilepath = argv[2];
  cout << "Writing polar order parameter to the following file: \n" << outfilepath << endl;
  write_single_analysis_data(polar_order_param, outfilepath);
  
  // deallocate the arrays
  
  for (int i = 0; i < nsteps; i++) {
    delete [] x[i];
    delete [] y[i];
  }
  delete [] x;
  delete [] y;
  
  return 0;
}  

//////////////////////////////////////////////////////////////////////////////////////////////////////////
