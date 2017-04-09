
/* calculate vorticity */

// COMPILATION AND RUN COMMANDS:
// /usr/local/gcc/bin/g++ -std=c++11 -O3 -lm -Wl,-rpath=$HOME/hdf5/lib -L$HOME/hdf5/lib -I$HOME/hdf5/include ${spath}/calc_vorticity.cpp ${upath}/read_write.cpp -lhdf5 -std=c++11 -o calc_vorticity
// ./calc_vorticity out.h5 ${path}/Vorticity_field.txt ${path}/Enstrophy.txt

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

inline int get_bin_number (double pos, double wbin, int nbins) {
  /* get the bin number based on the position */
  
  int ix = static_cast<int>(pos/wbin);
  if (ix < 0)
    ix = nbins-1;
    
  return ix;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void add_velocity_to_bins (int ix_bin, int iy_bin, int step,
			   vector<double> &vx_bin, vector<double> &vy_bin,
			   vector<int> &ncells_per_bin_per_step,
			   double vxc, double vyc, 
			   const int nbins, const int ncells) {

  int idx = step*nbins*nbins + ix_bin*nbins + iy_bin;

  vx_bin[idx] += vxc;
  vy_bin[idx] += vyc;
  ncells_per_bin_per_step[ix_bin*nbins+iy_bin]++;
  
  return; 
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void calc_velocity_per_bin (vector<double> &vx_bin, vector<double> &vy_bin,
			    vector<vector<int> > &ncells_per_bin,
			    const vector<double> &vx, const vector<double> &vy, 
			    double **x, double **y, 
			    double wbin, const int nbins, const double woverlap,
			    const double lx, const double ly, 
			    const int ncells, const int nsteps) {
  /* calculate velocities for each bin for all timesteps */

  for (int step = 0; step < nsteps; step++) {
    
    vector<int> ncells_per_bin_per_step(nbins*nbins, 0.0);  

    // calculate the velocities 
    
    for (int j = 0; j < ncells; j++) {
      
      double vx_cell = vx[step*ncells+j];
      double vy_cell = vy[step*ncells+j];
      
      // add to the central bin
      
      int ix_bin = get_bin_number(x[step][j], wbin, nbins);
      int iy_bin = get_bin_number(y[step][j], wbin, nbins);

      add_velocity_to_bins(ix_bin, iy_bin, step, vx_bin, vy_bin,
			   ncells_per_bin_per_step, vx_cell, vy_cell, 
			   nbins, ncells);

      // add to the neighboring bins because of overlaps
      
      int xsub = fmod(x[step][j], wbin);
      int ysub = fmod(y[step][j], wbin);
      
      int ixneigh, iyneigh;
      
      if (xsub > woverlap) {
	ixneigh = (ix_bin+1) % nbins;
      }
      else {
	ixneigh = (ix_bin-1) % nbins;
      }
      if (ysub > woverlap) {
	iyneigh = (iy_bin+1) % nbins;
      }
      else {
	iyneigh = (iy_bin-1) % nbins;
      }
      
      if (ixneigh < 0) 
	ixneigh = nbins-1;
      if (iyneigh < 0)
	iyneigh = nbins-1;
            
      add_velocity_to_bins(ixneigh, iy_bin, step, vx_bin, vy_bin,
			   ncells_per_bin_per_step, vx_cell, vy_cell, nbins, ncells);      
      add_velocity_to_bins(ix_bin, iyneigh, step, vx_bin, vy_bin,
			   ncells_per_bin_per_step, vx_cell, vy_cell, nbins, ncells);      
      add_velocity_to_bins(ixneigh, iyneigh, step, vx_bin, vy_bin,
			   ncells_per_bin_per_step, vx_cell, vy_cell, nbins, ncells);      
            
    }	  // cells
            
    // normalize the velocities per bin
    
    for (int i = 0; i < nbins; i++) {
      for (int j = 0; j < nbins; j++) {
	int idx = step*nbins*nbins + i*nbins + j;
	if (ncells_per_bin_per_step[i*nbins+j] > 0) {
	  vx_bin[idx] /= ncells_per_bin_per_step[i*nbins+j];
	  vy_bin[idx] /= ncells_per_bin_per_step[i*nbins+j];	
	}
      }	 // ybins
    }	// xbins
    
    ncells_per_bin[step] = ncells_per_bin_per_step;
    
  }	  // timesteps
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

double calc_vorticity (vector<double> &w_bin,
		    const vector<vector<int> > &ncells_per_bin,
		    const vector<double> &vx_bin, const vector<double> &vy_bin,
		    const double wbin, const int nbins, 
		    const double lx, const double ly, 
		    const int ncells, const int nsteps) {
  /* calculate the vorticity for all timesteps with finite difference differentiation */
  
  double enstrophy = 0.0;
  
  for (int step = 0; step < nsteps; step++) {
    
    double ens_per_step = 0.0;
        
    for (int i = 0; i < nbins; i++) {
      for (int j = 0; j< nbins; j++) {

	if (ncells_per_bin[step][i*nbins+j] > 0) {
	  double rho = wbin/2./ncells_per_bin[step][i*nbins+j];
	  
	  int xfi = (i+1) % nbins;
	  int xbi = (i-1) % nbins;
	  if (xbi < 0) 
	    xbi = nbins-1;
	  double wx = vy_bin[step*nbins*nbins+xfi*nbins+j] - vy_bin[step*nbins*nbins+xbi*nbins+j];
	  
	  int yfi = (j+1) % nbins;
	  int ybi = (j-1) % nbins;
	  if (ybi < 0) 
	    ybi = nbins-1;
	  double wy = vx_bin[step*nbins*nbins+i*nbins+yfi] - vx_bin[step*nbins*nbins+i*nbins+ybi];
	  
	  w_bin[step*nbins*nbins+i*nbins+j] = (wx-wy)*rho;
	  
	  ens_per_step += (wx-wy)*(wx-wy)*rho*rho;
	}
	
      }		// ybins
    }		// xbins
    
    enstrophy += ens_per_step/nbins/nbins;
    
  }		// timesteps
  
  enstrophy /= (2.*nsteps);
  
  return enstrophy;
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
  
  const int delta = 4;                   // number of data points between two steps
                                         // to calculate velocity
  const int nvels = nsteps-delta;        // number of data points in the velocity arrays
  
  // calculate the velocities per cell
  
  vector<double> vx(nvels*ncells, 0.0);
  vector<double> vy(nvels*ncells, 0.0);
  calc_velocity(vx, vy, x, y, lx, ly, delta, ncells, nvels, dt);
  
  // get all the image positions of cells
  
  get_img_pos(x, y, nsteps, ncells, lx, ly);
  
  // calculate the velocities per bin (square box and bin sizes are assumed!)
  
  double wbin = 18.; 			// bin width
  double woverlap = wbin*3./4.;		// bin overlap distance (currently %75)
  const int nbins = int(lx/wbin + 0.5);
  wbin = lx/nbins;
  vector<double> vx_bin(nvels*nbins*nbins, 0.0);
  vector<double> vy_bin(nvels*nbins*nbins, 0.0);
  vector<vector<int> > ncells_per_bin(nvels, vector<int>(nbins*nbins, 0.0));
  calc_velocity_per_bin(vx_bin, vy_bin, ncells_per_bin, 
			vx, vy, x, y, wbin, nbins, woverlap,
			lx, ly, ncells, nvels);
			
  // calculate the vorticities per bin
  

  vector<double> w_bin(nvels*nbins*nbins, 0.0);
  double enstrophy = calc_vorticity(w_bin, ncells_per_bin, 
				    vx_bin, vy_bin, wbin, nbins, 
				    lx, ly, ncells, nvels);  
  
  // write the computed data for vorticity field
  
  string outfilepath = argv[2];
  cout << "Writing vorticity field to the following file: \n" << outfilepath << endl;
  write_multid_analysis_data(w_bin, vx_bin, vy_bin, nbins, nvels, outfilepath);

    // write the computed data for enstrophy
  
  string outfilepath2 = argv[3];
  cout << "Writing enstrophy to the following file: \n" << outfilepath2 << endl;
  write_single_analysis_data(enstrophy, outfilepath2);
  
  // deallocate the arrays
  
  for (int i = 0; i < nsteps; i++) {
    delete [] x[i];
    delete [] y[i];
  }
  delete [] x;
  delete [] y;
  
  return 0;
}  
