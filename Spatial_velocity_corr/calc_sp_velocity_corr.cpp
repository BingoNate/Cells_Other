
/* calculate velocity correlation as a function of Euclidiean distance */

// COMPILATION AND RUN COMMANDS:
// g++ -O3 -lm -Wl,-rpath=$HOME/hdf5/lib -L$HOME/hdf5/lib -I$HOME/hdf5/include ${spath}/calc_sp_velocity_corr.cpp ${upath}/read_write.cpp -lhdf5 -o calc_sp_vel_corr
// ./calc_sp_vel_corr out.h5 ${path}/Sp_velocity_corr.txt

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
#include "../Utility/basic.hpp"
#include "../Utility/sim_info.hpp"

#define pi M_PI

// decide on which version to use
#define SUBTRACT_COM		// subtract center of mass velocities per frame

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void calc_velocity (vector<double> &vx, vector<double> &vy,
                    double **x, double **y,
                    const double lx, const double ly, const int delta,
                    const int ncells, const int nvels, const double dt) {
  /* calculate velocities with dt as delta */
  
  const double deltaDt = delta*dt;
  for (int i = 0; i < nvels; i++) {
    
    #ifdef SUBTRACT_COM
    long double comvx = 0.;
    long double comvy = 0.;
    #endif
    
    for (int j = 0; j < ncells; j++) {
      
      // note that UNWRAPPED COORDS ARE ASSUMED!
      
      double dx = x[i+delta][j] - x[i][j];
      vx[i*ncells+j] = dx/deltaDt;
      
      double dy = y[i+delta][j] - y[i][j];
      vy[i*ncells+j] = dy/deltaDt;
      
      #ifdef SUBTRACT_COM
      comvx += vx[i*ncells+j];
      comvy += vy[i*ncells+j];
      #endif
      
    }   // cell loop
    
    #ifdef SUBTRACT_COM
    comvx /= ncells;
    comvy /= ncells;    
    
    for (int j = 0; j < ncells; j++) {
      vx[i*ncells+j] -= comvx;
      vy[i*ncells+j] -= comvy;
    } 	// cells loop
    #endif
    
  }     // velocity step loop
  
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void calc_sp_vel_corr_v2 (double cvv[], 
                       const vector<double> &vx, const vector<double> &vy,
                       double **x, double **y,
                       double lx, double ly,
                       int ncells, int nsteps,
                       int delta, int ndata) {
  /* calculate spatial velocity correlation
  with the following definition (Wysocki, et. al.):
  Cvv(dr) = <v_i(r)*v_j(r+dr)>/(<v_i(r)^2/N>*(delta_ri*delta_rj))
  */
  
  // note that ndata here refers to the longest distance (longest_dist)
  // note that steps here refer to velocity data points
  // velocity time unit and position time unit are different
  
  double cnorm = 0.;
  
  for (int step = 0; step < nsteps; step++) {
        
    double cnorm_per_step = 0.;
    
    double cvv_per_step[ndata];
    int cnn_per_step[ndata];
    
    for (int i = 0; i < ndata; i++) {
      cvv_per_step[i] = 0.;  cnn_per_step[i] = 0;
    }
    
    for (int j1 = 0; j1 < ncells-1; j1++) {
      for (int j2 = j1+1; j2 < ncells; j2++) {
        
        // calculate the min. img. distance between the cells
        
        double dx = x[step*delta][j2] - x[step*delta][j1];
        dx = get_min_img_dist(dx, lx);
        double dy = y[step*delta][j2] - y[step*delta][j1];
        dy = get_min_img_dist(dy, ly);
        double dr = sqrt(dx*dx + dy*dy);
        int ibin = inearbyint(dr);
        
        cvv_per_step[ibin] += 2 * (vx[step*ncells+j1]*vx[step*ncells+j2] + vy[step*ncells+j1]*vy[step*ncells+j2]);
        cnn_per_step[ibin] += 2;

      } // inner cells loop
      
      cnorm_per_step += vx[step*ncells+j1]*vx[step*ncells+j1] + vy[step*ncells+j1]*vy[step*ncells+j1];
      
    }  // outer cells loop
    
    cnorm_per_step += vx[step*ncells+ncells-1]*vx[step*ncells+ncells-1] + vy[step*ncells+ncells-1]*vy[step*ncells+ncells-1]; 
    cnorm += cnorm_per_step;
    
    for (int i = 0; i < ndata; i++) {
      if (cnn_per_step[i] != 0) 
	cvv[i] += cvv_per_step[i]/cnn_per_step[i];
    }
    
  }  // timestep loop
  
  // normalize
  
  cnorm /= (nsteps*ncells);
  for (int j = 0; j < ndata; j++) {
    cvv[j] /= (nsteps*cnorm);
  }
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void calc_sp_vel_corr (double cvv[], 
                       const vector<double> &vx, const vector<double> &vy,
                       double **x, double **y,
                       double lx, double ly,
                       int ncells, int nsteps,
                       int delta, int ndata) {
  /* calculate spatial velocity correlation
  with the following definition (Wysocki, et. al.):
  Cvv(dr) = <v_i(r)*v_j(r+dr)>/(<v_i(r)^2/N>*(del_ri*del_rj))
  */
  
  // note that ndata here refers to the longest distance (longest_dist)
  // note that steps here refer to velocity data points
  // velocity time unit and position time unit are different
  
  
  for (int step = 0; step < nsteps; step++) {
        
    double cnorm_per_step[ndata];
    double cvv_per_step[ndata];
    int cnn_per_step[ndata];
    
    for (int i = 0; i < ndata; i++) {
      cvv_per_step[i] = 0.;  cnn_per_step[i] = 0;   cnorm_per_step[i] = 0.;
    }
    
    for (int j1 = 0; j1 < ncells-1; j1++) {
      for (int j2 = j1+1; j2 < ncells; j2++) {
        
        // calculate the min. img. distance between the cells
        
        double dx = x[step][j2] - x[step][j1];
        dx = get_min_img_dist(dx, lx);
        double dy = y[step][j2] - y[step][j1];
        dy = get_min_img_dist(dy, ly);
        double dr = sqrt(dx*dx + dy*dy);
        int ibin = inearbyint(dr);
        
        cvv_per_step[ibin] += 2*(vx[step*ncells+j1]*vx[step*ncells+j2] + vy[step*ncells+j1]*vy[step*ncells+j2]);
	cnorm_per_step[ibin] += vx[step*ncells+j1]*vx[step*ncells+j1] + vy[step*ncells+j1]*vy[step*ncells+j1] + vx[step*ncells+j2]*vx[step*ncells+j2] + vy[step*ncells+j2]*vy[step*ncells+j2];	
        cnn_per_step[ibin] += 1;

      } // inner cells loop
      
    }  // outer cells loop
    
    
    for (int i = 0; i < ndata; i++) {
      if (cnn_per_step[i] != 0 && cnorm_per_step[i] != 0.) {
	cvv[i] += cvv_per_step[i]/cnorm_per_step[i];
      }
    }
    
  }  // timestep loop
  
  // normalize
  
  for (int j = 0; j < ndata; j++) {
    cvv[j] /= nsteps;
  }
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

int main (int argc, char *argv[]) {

  // get the file name by parsing and load simulation data
  
  string filename = argv[1];
  cout << "Calculating spatial velocity correlation of the following file: \n" << filename << endl;
  SimInfo sim(filename);
  cout << "nsteps = " << sim.nsteps << endl;
  cout << "ncells = " << sim.ncells << endl;
 
  // read in the cell position data all at once
  
  /* the position data is stored in the following format:
  (nsteps, 2, ncells)
  the data will be loaded as follows:
  (nsteps, ncells) in x and y separately 
  */
  
  double **x = new double*[sim.nsteps];
  for (int i = 0; i < sim.nsteps; i++) x[i] = new double[sim.ncells];
  
  double **y = new double*[sim.nsteps];
  for (int i = 0; i < sim.nsteps; i++) y[i] = new double[sim.ncells];
  
  for (int i = 0; i < sim.nsteps; i++) {
    for (int j = 0; j < sim.ncells; j++) {
      x[i][j] = 0.; y[i][j] = 0.;
    }
  }
  
  read_all_pos_data(filename, x, y, sim.nsteps, sim.ncells, "/cells/comu");

  // set variables related to the analysis
  
  const int delta = 10;                       // number of data points between two steps
                                              // to calculate velocity
  const int nvels = sim.nsteps-delta;         // number of data points in the velocity array
  const int longest_dist = (int)(sim.lx+2);   // longest distance allowed by the sim. box
  double cvv[longest_dist];
  for (int i = 0; i < longest_dist; i++) 
    cvv[i] = 0.; 
  
  // calculate the velocities
  
  vector<double> vx(nvels*sim.ncells);
  vector<double> vy(nvels*sim.ncells);
  calc_velocity(vx, vy, x, y, sim.lx, sim.ly, delta, sim.ncells, nvels, sim.dt);
  
  // calculate the velocity correlation in space
  
  calc_sp_vel_corr (cvv, vx, vy, x, y, sim.lx, sim.ly, sim.ncells, nvels,
                    delta, longest_dist);
  

  // write the computed data
  
  string outfilepath = argv[2];
  cout << "Writing spatial velocity correlation to the following file: \n" << outfilepath << endl;
  write_1d_analysis_data(cvv, longest_dist, outfilepath);
  
  // deallocate the arrays
  
  for (int i = 0; i < nsteps; i++) {
    delete [] x[i];
    delete [] y[i];
  }
  delete [] x;
  delete [] y;
  
  return 0;
}  
