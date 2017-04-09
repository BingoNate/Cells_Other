
/* calculate energy spectrum */

// COMPILATION AND RUN COMMANDS:
// g++ -O3 -lm -Wl,-rpath=$HOME/hdf5/lib -L$HOME/hdf5/lib -I$HOME/hdf5/include ${spath}/calc_energy_spectrum.cpp ${upath}/read_write.cpp -lhdf5 -o calc_energy_spectrum
// ./calc_energy_spectrum out.h5 ${path}/Energy_spectrum.txt

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

void calc_sp_vel_corr (double cvv[], 
		    const vector<double> &vx, const vector<double> &vy,
		    double **x, double **y,
		    double lx, double ly,
		    int ncells, int nsteps,
		    int delta, int ndata) {
  /* calculate spatial velocity correlation PER STEP with 2D bins */
  
  // note that ndata here refers to the longest distance (longest_dist)
  // note that steps here refer to velocity data points
  
  for (int step = 0; step < nsteps; step++) {
        
    double cvv_per_step[ndata*ndata];
    int cnn_per_step[ndata*ndata];
    
    for (int i = 0; i < ndata*ndata; i++) {
      cvv_per_step[i] = 0.;  cnn_per_step[i] = 0;  
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
        int ixbin = inearbyint(dx);
        int iybin = inearbyint(dy);
        
        cvv_per_step[ixbin*ndata+iybin] += vx[step*ncells+j1]*vx[step*ncells+j2] + vy[step*ncells+j1]*vy[step*ncells+j2];
        cnn_per_step[ixbin*ndata+iybin] += 1;

      } // inner cells loop
      
    }  // outer cells loop
    
    
    for (int i = 0; i < ndata; i++) {
      for (int j = 0; j < ndata; j++) {
        if (cnn_per_step[i*ndata+j] != 0) {
          cvv[step*ndata*ndata+i*ndata+j] = cvv_per_step[i*ndata+j]/cnn_per_step[i*ndata+j];
        }
      }
    }
    
  }  // timestep loop
  
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void calc_energy_spectrum (double ek[], double kvec[],
			  const double cvv[],
			  double lx, double ly,
			  int ncells, int nsteps,
			  int delta, int ndata, int Nmax) {
  /* calculate the energy spectrum */
  
  // note that ndata here refers to the longest distance (longest_dist)
  // note that steps here refer to velocity data points
  
  // set preliminary data
  
  double delk = 2*pi/lx;
  double lamda_min = 4.; 
  double lamda_max = lx/2.;
  double kmax = 2*pi/lamda_min;
  double kmin = 2*pi/lamda_max;
  int nks = 24;
   
  kmax = log(kmax);
  kmin = log(kmin);  
  double wbin = (kmax-kmin)/Nmax;

  cout << "delta k = " << delk << "\n" << "wbin = " << wbin << "\n" << "kmin = " << exp(kmin) << "\n" << "kmax = " << exp(kmax) << "\n" << "Nmax = " << Nmax << "\n" << endl;
  
  // populate log-spaced absolute value of the k vectors
  
  for (int j = 0; j < Nmax; j++) kvec[j] = kmin + j*wbin;
  for (int j = 0; j < Nmax; j++) kvec[j] = exp(kvec[j]);

  // populate circular orientation of the k vector
  
  double kxs[nks];
  double kys[nks];
  for (int j = 0; j < nks; j++) {
    kxs[j] = cos(2*pi*j/nks);
    kys[j] = sin(2*pi*j/nks);
  }

  for (int step = 0; step < nsteps; step++) {
    
    for (int k = 0; k < Nmax; k++) {
      
      for (int angle = 0; angle < nks; angle++) {
      
        double cost = 0.0;
        double sint = 0.0;
        for (int i = 0; i < ndata; i++) {
          for (int j = 0; j < ndata; j++) {
    
            double dotp = kvec[k]*kxs[angle]*i + kvec[k]*kys[angle]*j;
            cost += cos(dotp)*cvv[step*ndata*ndata+i*ndata+j];
            sint += sin(dotp)*cvv[step*ndata*ndata+i*ndata+j];
            
          }  // y distance loop
        } // x distance loop

        ek[k] += sqrt(cost*cost + sint*sint);
        
      }  // circular orientation of kvector loop
	
    } // absolute value of kvector loop
  
  } // timesteps
  
  // perform the normalization
  
  for (int j = 0; j < Nmax; j++)  { 
    ek[j] = ek[j]*kvec[j]/(nsteps*2*pi*nks);
  }  
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

int main (int argc, char *argv[]) {

  // get the file name by parsing
  
  string filename = argv[1];
  cout << "Calculating energy spectrum of the following file: \n" << filename << endl;

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
  const int nvels = nsteps-delta;       // number of data points in the velocity array
  const int longest_dist = (int)lx+2;   // longest distance allowed by the sim. box
  double cvv[nvels*longest_dist*longest_dist];
  for (int i = 0; i < nvels*longest_dist*longest_dist; i++)
    cvv[i] = 0.; 
  
  // calculate the velocities
  
  vector<double> vx(nvels*ncells);
  vector<double> vy(nvels*ncells);
  calc_velocity(vx, vy, x, y, lx, ly, delta, ncells, nvels, dt);
  
  // calculate the velocity correlation in space
  
  calc_sp_vel_corr (cvv, vx, vy, x, y, lx, ly, ncells, nvels,
                    delta, longest_dist);

  // set variables related to energy spectrum
  
  int Nmax = 100;
  double *kvec = new double[Nmax];
  double *ek = new double[Nmax];
  for (int j = 0; j < Nmax; j++) {
    ek[j] = 0.;  kvec[j] = 0.;
  }

  // calculate the energy spectrum in Fourier space
  
  calc_energy_spectrum (ek, kvec, cvv, lx, ly, ncells, nvels,
			delta, longest_dist, Nmax);
		    
  // write the computed data
  
  string outfilepath = argv[2];
  cout << "Writing energy spectrum n to the following file: \n" << outfilepath << endl;
  write_2d_analysis_data(kvec, ek, Nmax, outfilepath);
  
  // deallocate the arrays
  
  for (int i = 0; i < nsteps; i++) {
    delete [] x[i];
    delete [] y[i];
  }
  delete [] x;
  delete [] y;
  
  return 0;
}  
