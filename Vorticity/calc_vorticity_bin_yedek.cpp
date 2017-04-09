
/* calculate vorticity */

// COMPILATION AND RUN COMMANDS:
// /usr/local/gcc/bin/g++ -O3 -lm -Wl,-rpath=$HOME/hdf5/lib -L$HOME/hdf5/lib -I$HOME/hdf5/include ${spath}/calc_vorticity.cpp ${upath}/read_write.cpp -lhdf5 -std=c++11 -o calc_vorticity
// ./calc_vorticity out.h5 ${path}/Vorticity_field.h5

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

#define pi M_PI

// decide on which version to use
//#define SUBTRACT_COM		// subtract center of mass velocities per frame

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
      
      int curr_step = i;
      int next_step = curr_step + delta;
      
      // note that UNWRAPPED COORDS ARE ASSUMED!
      
      double dx = x[next_step][j] - x[curr_step][j];
      vx[i*ncells+j] = dx/deltaDt;
      
      double dy = y[next_step][j] - y[curr_step][j];
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
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

inline int get_bin_number (double pos, double wbin) {
  /* get the bin number based on the position */
  
  return pos/wbin;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void calc_velocity_per_bin (vector<double> &vx_bin, vector<double> &vy_bin,
			    vector<vector<int> > &ncells_per_bin,
			    vector<int> &superbins,
			    //vector<vector<vector<double> > > &vx_bin, 
			    //vector<vector<vector<double> > > &vy_bin,
			    const vector<double> &vx, const vector<double> &vy, 
			    double **x, double **y, 
			    const double wbin, const int nbins, 
			    const double w_superbin, const int n_superbins,
			    const double lx, const double ly, 
			    const int ncells, const int nsteps) {
  /* calculate velocities for each bin for all timesteps */
  
  for (int step = 0; step < nsteps; step++) {
    
    vector<int> ncells_per_step(nbins*nbins, 0.0);

    // calculate the velocities 
    
    for (int j = 0; j < ncells; j++) {
      
      int x_bin = get_bin_number(x[step][j], wbin);
      int y_bin = get_bin_number(y[step][j], wbin); 
      
      int idx = step*nbins*nbins + x_bin*nbins + y_bin;
      vx_bin[idx] += vx[step*ncells+j];
      vy_bin[idx] += vy[step*ncells+j];
      ncells_per_step[x_bin*nbins+y_bin] += 1;
      
      int x_superbin = get_bin_number(x[step][j], w_superbin);
      int y_superbin = get_bin_number(y[step][j], w_superbin);        
      int idx_superbin = step*n_superbins*n_superbins + x_superbin*n_superbins + y_superbin;
      superbins[idx] = idx_superbin;
            
    }	  // cells
    
    // normalize the velocities per bin
    
    for (int i = 0; i < nbins; i++) {
      for (int j = 0; j < nbins; j++) {
	int idx = step*nbins*nbins + i*nbins + j;
	if (ncells_per_step[i*nbins+j] > 0) {
	  vx_bin[idx] /= ncells_per_step[i*nbins+j];
	  vy_bin[idx] /= ncells_per_step[i*nbins+j];	
	}
      }	 // ybins
    }	// xbins
    
    ncells_per_bin[step] = ncells_per_step;
    
  }	  // timesteps
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void calc_vorticity (vector<double> &wx_bin, vector<double> &wy_bin,
		    const vector<int> &superbins,
		    const vector<vector<int> > &ncells_per_bin,
		    const vector<double> &vx_bin, const vector<double> &vy_bin,
		    const double wbin, const int nbins, 
		    const double lx, const double ly, 
		    const int ncells, const int nsteps) {
  /* calculate the vorticity for all timesteps with finite difference differentiation */
  
  // histogram the vorticities from the binned velocities
  
  for (int step = 0; step < nsteps; step++) {
        
    for (int i = 0; i < nbins; i++) {
      for (int j = 0; j< nbins; j++) {
	
	if (ncells_per_bin[step][i*nbins+j] > 0) {
	  double rho = wbin/2./ncells_per_bin[step][i*nbins+j];
	  
	  int idx = step*nbins*nbins+i*nbins+j;
	  int idx_superbin = superbins[idx];
	  
	  int xfi = (i+1) % nbins;
	  int xbi = (i-1) % nbins;
	  wx_bin[idx_superbin] += vy_bin[step*nbins*nbins+xfi*nbins+j] - vy_bin[step*nbins*nbins+xbi*nbins+j];
	  
	  int yfi = (j+1) % nbins;
	  int ybi = (j-1) % nbins;
	  wy_bin[idx_superbin] += -vx_bin[step*nbins*nbins+i*nbins+yfi] + vx_bin[step*nbins*nbins+i*nbins+ybi];
	}
	
      }		// ybins
    }		// xbins
  }		// timesteps
  
  //
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

int main (int argc, char *argv[]) {

  // get the file name by parsing
  
  string filename = argv[1];
  cout << "Calculating vorticity of the following file: \n" << filename << endl;

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
  
  const int delta = 5;                   // number of data points between two steps
                                         // to calculate velocity
  const int nvels = nsteps-delta;        // number of data points in the velocity arrays
  
  // calculate the velocities per cell
  
  vector<double> vx(nvels*ncells, 0.0);
  vector<double> vy(nvels*ncells, 0.0);
  calc_velocity(vx, vy, x, y, lx, ly, delta, ncells, nvels, dt);
  
  // get all the image positions of cells
  
  get_img_pos(x, y, nsteps, ncells, lx, ly);
    
  // generate bins of cells to calculate velocities (square box and bin sizes are assumed!)
  
  double wbin = 4.; 			// bin width
  const int nbins = int(lx/wbin + 0.5);
  wbin = lx/nbins;
  vector<double> vx_bin(nvels*nbins*nbins, 0.0);
  vector<double> vy_bin(nvels*nbins*nbins, 0.0);
  vector<vector<int> > ncells_per_bin(nvels, vector<int>(nbins*nbins, 0.0));

  // calculate the velocoties per bin
  
  calc_velocity_per_bin(vx_bin, vy_bin, ncells_per_bin, 
			vx, vy, x, y, wbin, nbins, 
			lx, ly, ncells, nvels);
			
  // generate an outer bin/superbin of velocities to calculate vorticity
  
  double w_outer_bin = 5*w;
  const int n_outer_bins = int(lx/w_outer_bin + 0.5);
  w_outer_bin = lx/n_outer_bins;
  
  // calculate the vorticities per superbin
  
  vector<double> wx_bin(nvels*n_outer_bins*n_outer_bins, 0.0);
  vector<double> wy_bin(nvels*n_outer_bins*n_outer_bins, 0.0);
  vector<vector<int> > nvels_per_bin(nvels, vector<int>(n_outer_bins*n_outer_bins, 0.0));
  calc_vorticity(wx_bin, wy_bin, nvels_per_bin, 
		vx_bin, vy_bin, wbin, nbins, 
		lx, ly, ncells, nvels); 

//   // calculate the vorticities per bin
//   
//   vector<double> wx_bin(nvels*nbins*nbins, 0.0);
//   vector<double> wy_bin(nvels*nbins*nbins, 0.0);
//   calc_vorticity(wx_bin, wy_bin, ncells_per_bin, 
// 		vx_bin, vy_bin, wbin, nbins, 
// 		lx, ly, ncells, nvels);  
//   
  

  // write the computed data for vorticity field
  
  string outfilepath = argv[2];
  cout << "Writing vorticity field to the following file: \n" << outfilepath << endl;
  write_multid_analysis_data(wx_bin, wy_bin, vx_bin, vy_bin, nbins, nvels, outfilepath);
  
  // deallocate the arrays
  
  for (int i = 0; i < nsteps; i++) {
    delete [] x[i];
    delete [] y[i];
  }
  delete [] x;
  delete [] y;
  
  return 0;
}  
