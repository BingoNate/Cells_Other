
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
#include <fftw3.h>
#include <complex.h>

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
  
  cout << "Calculating the velocities" << endl;
  
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

void calc_energy_spectrum (double ek[], double kvec[],
			  const double cvv[],
			  double lx, double ly,
			  int ncells, int nsteps,
			  int delta, int ndata, int Nmax,
			  double kmax, double kmin, int njump, double delk) {
  /* calculate the energy spectrum */
  
  cout << "Calculating the energy spectrum" << endl;
  
  // note that ndata here refers to the longest distance (longest_dist)
  // note that steps here refer to velocity data points
  
  // set preliminary data

  int nks = 24;
  double wbin = (kmax-kmin)/Nmax;

  cout << "delta k = " << delk << "\n" << "wbin = " << wbin << "\n" << "kmin = " << exp(kmin) << "\n" << "kmax = " << exp(kmax) << "\n" << "Nmax = " << Nmax << "\n" << endl;
//   cout << "delta k = " << delk << "\n" << "wbin = " << wbin << "\n" << "kmin = " << kmin << "\n" << "kmax = " << kmax << "\n" << "Nmax = " << Nmax << "\n" << endl;
  
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

  for (int step = 0; step < nsteps; step+=njump) {
    
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
    ek[j] = ek[j]*kvec[j]*njump/(nsteps*2*pi*nks);
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
  
  const int delta = 10;                     // number of data points between two steps
                                            // to calculate velocity
  const int nvels = nsteps-delta;           // number of data points in the velocity array
  const int longest_dist = (int)(lx/2+1);   // longest distance allowed by the sim. box
  double *cvv = new double[nvels*longest_dist*longest_dist];
  for (int i = 0; i < nvels*longest_dist*longest_dist; i++)
    cvv[i] = 0.; 
  
  // calculate the velocities
  
  vector<double> vx(nvels*ncells);
  vector<double> vy(nvels*ncells);
  calc_velocity(vx, vy, x, y, lx, ly, delta, ncells, nvels, dt);
  
  // calculate the velocity autocorrelation in k space
  
  //calc_vel_autocorr_in_kspace();
  
  int arr_size = ncells/2 + 1;
  unsigned flags = 0;
  double *ek = new double[arr_size];
  double *kvec = new double[arr_size];
  for (int j = 0; j < arr_size; j++) {
    kvec[j] = 2.*pi*j/arr_size;
  }
  for (int tstep = 0; tstep < nsteps; tstep++) {
    
    double *v_per_step = (double *) fftw_malloc(sizeof(double)*2*ncells);
    for (int j = 0; j < ncells; j++) {
      v_per_step[j] = vx[j]; 
      v_per_step[ncells+j] = vy[j];
    }
    
    fftw_complex *vk_per_step = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*2*arr_size);
    fftw_plan p = fftw_plan_dft_r2c_2d(2, ncells, v_per_step, vk_per_step, FFTW_ESTIMATE);
    fftw_execute(p);
    
    for (int j = 0; j < arr_size; j++) {
      double real_part = vk_per_step[j][0] + vk_per_step[arr_size+j][0];
      double img_part = vk_per_step[j][1] + vk_per_step[arr_size+j][1];
      cout << vk_per_step[j][0] << "\t" << vk_per_step[arr_size+j][0] << endl;
      ek[j] += real_part*real_part + img_part*img_part;
    }
    
    fftw_destroy_plan(p);
    fftw_free(vk_per_step);
    fftw_free(v_per_step);
    
  }  // timestep loop
  
  for (int j = 0; j < arr_size; j++) {
    ek[j] = 4*pi*ek[j]*kvec[j]/nsteps;
  }

//   // set variables related to energy spectrum
// 
//   double delk = 2*pi/lx;
//   double lamda_min = 3.; 
//   double lamda_max = lx;
//   double kmax = 2*pi/lamda_min;
//   double kmin = 2*pi/lamda_max;
//   int nks = 24;
//   int njump = 20;
//   kmax = log(kmax);
//   kmin = log(kmin);  
//   int Nmax = static_cast<int>((kmax-kmin)/delk);
//   double wbin = (kmax-kmin)/Nmax;
//   double *kvec = new double[Nmax];
//   double *ek = new double[Nmax];
//   for (int j = 0; j < Nmax; j++) {
//     ek[j] = 0.;  kvec[j] = 0.;
//   }
// 
//   // calculate the energy spectrum in Fourier space
//   
//   calc_energy_spectrum (ek, kvec, cvv, lx, ly, ncells, nvels,
// 			delta, longest_dist, Nmax,
// 			kmax, kmin, njump, delk);
		    
  // write the computed data
  
  string outfilepath = argv[2];
  cout << "Writing energy spectrum to the following file: \n" << outfilepath << endl;
  write_2d_analysis_data(kvec, ek, arr_size, outfilepath);
  
  // deallocate the arrays
  
  delete [] cvv;
  for (int i = 0; i < nsteps; i++) {
    delete [] x[i];
    delete [] y[i];
  }
  delete [] x;
  delete [] y;
  
  return 0;
}  
