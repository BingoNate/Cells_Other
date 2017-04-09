
/* calculate center of mass velocities */

// COMPILATION AND RUN COMMANDS:
// g++ -O3 -lm -Wl,-rpath=$HOME/hdf5/lib -L$HOME/hdf5/lib -I$HOME/hdf5/include ${spath}/calc_com_velocity.cpp ${upath}/read_write.cpp -lhdf5 -o calc_com_vel
// ./calc_com_vel out.h5 ${path}/Velocity_com.txt

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

double calc_velocity (vector<double> &vx, vector<double> &vy,
                    double **x, double **y,
                    const double lx, const double ly, const int delta,
                    const int ncells, const int nvels, const double dt) {
  /* calculate velocities with dt as delta */
  
  const double deltaDt = delta*dt;
  double comv = 0.;
  for (int i = 0; i < nvels; i++) {
    
    long double comvx = 0.;
    long double comvy = 0.;
    
    for (int j = 0; j < ncells; j++) {
      
      int curr_step = i*delta;
      int next_step = curr_step + delta;
      
      // note that UNWRAPPED COORDS ARE ASSUMED!
      
      double dx = x[next_step][j] - x[curr_step][j];
      vx[i*ncells+j] = dx/deltaDt;
      
      double dy = y[next_step][j] - y[curr_step][j];
      vy[i*ncells+j] = dy/deltaDt;
      
      comvx += vx[i*ncells+j];
      comvy += vy[i*ncells+j];
      
    }   // cell loop
    
    comvx /= ncells;
    comvy /= ncells;    
    
    for (int j = 0; j < ncells; j++) {
      vx[i*ncells+j] -= comvx;
      vy[i*ncells+j] -= comvy;
    } 	// cells loop
    
    comv += sqrt(comvx*comvx + comvy*comvy);
    
  }     // velocity step loop
  
  comv /= nvels;
  
  return comv;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

int main (int argc, char *argv[]) {

  // get the file name by parsing
  
  string filename = argv[1];
  cout << "Calculating c.o.m. velocity of the following file: \n" << filename << endl;

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
  
  for (int i = 0; i < nsteps; i++) {
    for (int j = 0; j < ncells; j++) {
      x[i][j] = 0.;  y[i][j] = 0.;
    }
  }
  
  read_all_pos_data(filename, x, y, nsteps, ncells, "/cells/comu");

  // set variables related to the analysis
  
  const int delta = 20;                 // number of data points between two steps
                                        // to calculate velocity
  const int nvels = nsteps/delta;       // number of data points in the velocity array
  
  // calculate the velocities
  
  vector<double> vx(nvels*ncells);
  vector<double> vy(nvels*ncells);
  double com_vel = calc_velocity(vx, vy, x, y, lx, ly, delta, ncells, nvels, dt);
 
  // write the computed data
  
  string outfilepath = argv[2];
  cout << "Writing c.o.m. velocity to the following file: \n" << outfilepath << endl;
  write_single_analysis_data(com_vel, outfilepath);
  
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
