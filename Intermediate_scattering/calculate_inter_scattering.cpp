
/* calculate intermediate scattering function */

// COMPILATION AND RUN COMMANDS:
// g++ -Wl,-rpath=$HOME/hdf5/lib -L$HOME/hdf5/lib -I$HOME/hdf5/include ${spath}/calculate_inter_scattering.cpp ${upath}/read_write.cpp -lhdf5 -fopenmp -Wno-write-strings -o calc_inter_scattering
// ./calc_inter_scattering out.h5 ${path}/Inter_scattering.txt

//////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include "omp.h"
#include "../Utility/read_write.hpp"

#define pi M_PI

// decide on which version to use
#define SUBTRACT_COM		// subtract center of mass velocities per frame

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void compute_inter_scattering (double **x, double **y, double *Fs, double *delay, int ndelay, double qvector, int nsteps, int ncells, double dt) {
  /* calculate the self part of the intermedaite scattering function */
  
  // populate the circular orientation of the wavevector
  
  int nks = 24;
  double kxs[nks], kys[nks];
  for (int j = 0; j < nks; j++) { 
    kxs[j] = cos(2*pi*j/nks)*qvector;
    kys[j] = sin(2*pi*j/nks)*qvector;
  }

  // allocate the arrays
  
  double avg_over_qvector[ndelay];
  for (int j = 0; j < ndelay; j++) {
    avg_over_qvector[j] = 0.;
  }
  
  // calculate the intermediate scattering function
  
  delay[0] = 0.;
  Fs[0] = 1.;
  
  for (int d = 1; d < ndelay; d++) {
        
      // sum over all wavevectors with a circular discretization
      
      for (int j = 0; j < nks; j++) {
	
	double term_to_avg = 0.;
	
	// sum over different time origins
	
	for (int step = 0; step < nsteps-ndelay; step++) {
	  
	  double cost = 0.;
	  //double sint = 0.;

	  omp_set_num_threads(4);

	  #ifdef SUBTRACT_COM
	  // calculate center of mass drift
	  
	  double dcomx = 0.;
	  double dcomy = 0.;
	  #pragma omp parallel for reduction(+:dcomx,dcomy)
	  for (int p = 0; p < ncells; p++) {
	    double dx = x[d+step][p] - x[step][p];
	    double dy = y[d+step][p] - y[step][p];
	    dcomx += dx;
	    dcomy += dy;
	  }
	  
	  dcomx /= ncells;
	  dcomy /= ncells;
	  #endif
	  
	  // sum over particles
	  
	  //#pragma omp parallel for reduction(+:cost,sint)	
	  #pragma omp parallel for reduction(+:cost)	  
	  for (int p = 0; p < ncells; p++) {
	    
	    double dx = x[d+step][p] - x[step][p];
	    double dy = y[d+step][p] - y[step][p];
	    #ifdef SUBTRACT_COM
	    dx -= dcomx;
	    dy -= dcomy;
	    #endif
	    double dotp = kxs[j]*dx + kys[j]*dy;
	    cost += cos(dotp);
	    //sint += sin(dotp);
	    
	  } // particles
	  
	  //term_to_avg += sqrt(cost*cost + sint*sint);
	  term_to_avg += cost;
	
	} // time origins
	
	term_to_avg /= (ncells*(nsteps-ndelay));
	avg_over_qvector[d] += term_to_avg;
	
      } // wavectors
   
  } // delay
  
  // perform the normalizations
  
  for (int d = 1; d < ndelay; d++) {
    
    delay[d] = d*dt; 
    Fs[d] = avg_over_qvector[d]/nks;
    
  }

  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

int main (int argc, char *argv[]) {

  // get the file name by parsing
  
  char *filename = argv[1];
  cout << "Calculating intermediate scattering function of the following file: \n" << filename << endl;

  // read in general simulation data
  
  int nsteps, nbeads, nsamp, ncells;
  nsteps = nbeads = nsamp = ncells = 0;
  double lx, ly, dt, eps, rho, fp, areak, bl, sigma;
  lx = ly = dt = eps = rho = fp = areak = bl = sigma = 0.;
  
  read_sim_data(filename, nsteps, nbeads, nsamp, ncells, lx, ly, dt, eps, rho, fp, areak, bl, sigma);
  
  // print simulation information
  
  cout << "nsteps = " << nsteps << endl;
  cout << "ncells = " << ncells << endl;
  
  // read in array data
  
  int nbpc[ncells];
  for (int i = 0; i < ncells; i++) nbpc[i] = 0;
  //read_integer_array(filename, "/cell/nbpc", nbpc);
 
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
  
  // allocate the arrays
  
  int ndelay = nsteps/2;
  double delay[ndelay];
  double Fs[ndelay];
  for (int j = 0; j < ndelay; j++) {
    delay[j] = 0.; Fs[j] = 0.; 
  }  
  
  // calculate the intermediate scattering function
  
  double qvector = 2*pi/7.0;
  compute_inter_scattering(x, y, Fs, delay, ndelay, qvector, nsteps, ncells, dt);  
  
  // write the computed data
  
  char *outfilepath = argv[2];
  cout << "Writing intermediate scattering function to the following file: \n" << outfilepath << endl;  
  write_2d_analysis_data(delay, Fs, (double)ndelay, outfilepath);
  
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
