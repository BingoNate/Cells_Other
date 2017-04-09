
/* calculate velocity correlation length as a function of time separation */

// COMPILATION AND RUN COMMANDS:
// g++ -std=c++11 -O3 -lm -Wl,-rpath=$HOME/hdf5/lib -L$HOME/hdf5/lib -I$HOME/hdf5/include ${spath}/calc_vel_corr_length.cpp ${upath}/read_write.cpp -lhdf5 -o calc_vel_corr_length
// ./calc_vel_corr_length out.h5 ${path}/Sp_velocity_corr ${path}/Vel_corr_length.txt

// NOTE THAT there are 2 txt files given and the usual .txt is NOT written to the first one!!
// NOTE THAT unlike the usual case, this script requires c++11 because of some string opeartions

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
#define SUBTRACT_COM		// subtract center of mass velocities per frame

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void calc_velocity (vector<double> &vx, vector<double> &vy,
                    double **x, double **y,
                    double lx, double ly, int delta,
                    int ncells, int nvels, double dt) {
  /* calculate velocities with dt as delta */
  
  const double deltaDt = delta*dt;
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
        
        double dx = x[step*delta][j2] - x[step*delta][j1];
        dx = get_min_img_dist(dx, lx);
        double dy = y[step*delta][j2] - y[step*delta][j1];
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

double calc_corr_length (const double cvv[], const int longest_dist, 
			 const double lx, const double ly, const int ncells) {
  /* calculate correlation length from spatial velocity correlation
  for reference, see Wysocki, et. al., EPL */
  
  double corr_len = 0.;
  for (int len = 0; len < longest_dist/2; len++) {
    corr_len += 2*pi*len*cvv[len];
  }
  corr_len *= (ncells/(lx*ly));
  
  return corr_len;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

int main (int argc, char *argv[]) {

  // get the file name by parsing
  
  string filename = argv[1];
  cout << "Calculating velocity correlation length of the following file: \n" << filename << endl;

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
  // note that delta refers to time separation between two data points

  const int longest_dist = (int)lx+2;   // longest distance allowed by the sim. box

  const int ddelta = 5;			// separation between deltas
  const int ndelta = 30;		// total number of deltas
  double delta[ndelta];
  double corr_len[ndelta];
  for (int i = 0; i < ndelta; i++) {
    delta[i] = (i+1)*ddelta;
    corr_len[i] = 0.;
  }
  
  // run over deltas to calculate the spatial velocity correlation
  
  for (int i = 0; i < ndelta; i++) {
    
    int nvels = nsteps/delta[i];       // number of data points in the velocity array
    double cvv[longest_dist]; 
    for (int i = 0; i < longest_dist; i++) 
      cvv[i] = 0.;  
    
    // calculate the velocities with the given delta
    
    vector<double> vx(nvels*ncells);
    vector<double> vy(nvels*ncells);
    calc_velocity(vx, vy, x, y, lx, ly, delta[i], ncells, nvels, dt);
    
    // calculate the velocity correlation in space
    
    calc_sp_vel_corr (cvv, vx, vy, x, y, lx, ly, ncells, nvels,
		      (int)delta[i], longest_dist);

    // write the computed velocity correlation data for the given delta
    
    string outfilepath = argv[2];
    outfilepath += "_delta_" + to_string((int)delta[i]) + ".txt";
    cout << "Writing spatial velocity correlation to the following file: \n" << outfilepath << endl;
    write_1d_analysis_data(cvv, longest_dist, outfilepath);
    
    // calculate the correlation length for the given delta
    
    corr_len[i] = calc_corr_length(cvv, longest_dist, lx, ly, ncells); 
  
  } 	// delta loop
  
  // write the computed correlation length data 

  for (int i = 0; i < ndelta; i++) {
    delta[i] *= dt;	// convert deltas into time units
  }
  string outfilepath_2 = argv[3];
  cout << "Writing correlation length to the following file: \n" << outfilepath_2 << endl;  
  write_2d_analysis_data(delta, corr_len, (double)ndelta, outfilepath_2);
  
  // deallocate the arrays
  
  for (int i = 0; i < nsteps; i++) {
    delete [] x[i];
    delete [] y[i];
  }
  delete [] x;
  delete [] y;
  
  return 0;
}  
