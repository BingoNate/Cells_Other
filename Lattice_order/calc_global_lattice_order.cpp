
/* calculate lattice order (hex order, square lattice order, etc.) */

// COMPILATION AND RUN COMMANDS:
// g++ -Wl,-rpath=$HOME/hdf5/lib -L$HOME/hdf5/lib -I$HOME/hdf5/include ${spath}/calc_global_lattice_order.cpp ${upath}/read_write.cpp -lhdf5 -Wno-write-strings -o calc_glob_lat
// ./calc_glob_lat out.h5 ${path}/Glob_lattice_order.txt

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
#include "omp.h"
#include "../Utility/read_write.hpp"
#include "../Utility/kdtree.hpp"
#include "../Utility/cartesian.hpp"

#define pi M_PI

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void populate_extended_pos_array(vector<Point<2> > &points, const int npoints, const int ncells, const int step, double **x, double **y, double lx, double ly) {
  /* populate an extended array with all the image points */
    
  for (int i = 0; i < ncells; i++) {
    points[i].x[0] = x[step][i];
    points[i].x[1] = y[step][i];
    
    points[i+ncells].x[0] = x[step][i] - lx;
    points[i+ncells].x[1] = y[step][i] + ly;
    
    points[i+2*ncells].x[0] = x[step][i];
    points[i+2*ncells].x[1] = y[step][i] + ly;
    
    points[i+3*ncells].x[0] = x[step][i] + lx;
    points[i+3*ncells].x[1] = y[step][i] + ly; 
    
    points[i+4*ncells].x[0] = x[step][i] + lx;
    points[i+4*ncells].x[1] = y[step][i];
    
    points[i+5*ncells].x[0] = x[step][i] + lx;
    points[i+5*ncells].x[1] = y[step][i] - ly;  
    
    points[i+6*ncells].x[0] = x[step][i];
    points[i+6*ncells].x[1] = y[step][i] - ly; 
    
    points[i+7*ncells].x[0] = x[step][i] - lx;
    points[i+7*ncells].x[1] = y[step][i] - ly; 
    
    points[i+8*ncells].x[0] = x[step][i] - lx;
    points[i+8*ncells].x[1] = y[step][i];     
  }
  
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

double calc_glob_order (double **x, double **y, double lx, double ly, const int nn, const int ncells, const int nsteps) {
  /* calculate global lattice order for nn */

  const int npoints = 9*ncells;		// size of the extended array with all the periodic images
  double order_param = 0.; 
  
  for (int step = 0; step < nsteps; step++) {
  
    // populate a vector of points in an extended way with all the image particles around the central box per time step

    vector<Point<2> > points(npoints);
    populate_extended_pos_array(points, npoints, ncells, step, x, y, lx, ly);
    
    // build the kdtree 
    
    KDtree<2> kdtree(points);
    
    // calculate the global lattice order
    
    long double cost = 0.;
    long double sint = 0.;
    
    for (int j = 0; j < ncells; j++) {
      
      // find the nn number of nearest neighbors
      
      int neighbor_idx[nn];
      double neighbor_pos[nn];
      for (int ns = 0; ns < nn; ns++) {
	neighbor_idx[ns] = 0;   neighbor_pos[ns] = 0.;
      }
      kdtree.nnearest(j, neighbor_idx, neighbor_pos, nn);
      
      for (int k = 0; k < nn; k++) {
	
	double dx = dist_xy(points[neighbor_idx[k]], points[j], 0);
	double dy = dist_xy(points[neighbor_idx[k]], points[j], 1);

	double theta = nn*atan2(dy, dx);
	
	cost += cos(nn*theta);
	sint += sin(nn*theta);
	
      } // neighbors loop
      
    }  // cells loop
    
    cost /= (nn*ncells);
    sint /= (nn*ncells);
    
    order_param += cost*cost + sint*sint;
  
  }  // timestep loop

  order_param /= nsteps;
  
  return order_param;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

int main (int argc, char *argv[]) {

  // get the file name by parsing
  
  char *filename = argv[1];
  cout << "Calculating lattice order of the following file: \n" << filename << endl;

  // read in general simulation data
  
  int nsteps, nbeads, nsamp, ncells;
  nsteps = nbeads = nsamp = ncells = 0;
  double lx, ly, dt, eps, rho, fp, areak, kappa, bl, sigma;
  lx = ly = dt = eps = rho = fp = areak = kappa = bl = sigma = 0.;
  
  read_sim_data(filename, nsteps, nbeads, nsamp, ncells, lx, ly, dt, eps, rho, fp, areak, kappa, bl, sigma);
  
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
  
  for (int i = 0; i < nsteps; i++) {
    for (int j = 0; j < ncells; j++) {
      x[i][j] = 0.; y[i][j] = 0.;
    }
  }
  
  read_all_pos_data(filename, x, y, nsteps, ncells, "/cells/comu");
  
  // get the image positions in the central box for all time steps at once
  
  get_img_pos(x, y, nsteps, ncells, lx, ly);

  // set variables related to the analysis
  
  const int nn = 6;			// number of nearest neighbor points to investigate
  double order_param = calc_glob_order(x, y, lx, ly, nn, ncells, nsteps);
  
  // write the computed data
  
  char *outfilepath = argv[2];
  cout << "Writing lattice order parameter to the following file: \n" << outfilepath << endl;  
  write_single_analysis_data(order_param, outfilepath);
  
  // deallocate the arrays
  
  for (int i = 0; i < nsteps; i++) {
    delete [] x[i];
    delete [] y[i];
  }
  delete [] x;
  delete [] y;
  
  return 0;
}  
