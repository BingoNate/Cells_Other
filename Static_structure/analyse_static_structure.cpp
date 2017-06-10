
/* calculate static structure factor */

// COMPILATION AND RUN COMMANDS:
// g++ -Wl,-rpath=$HOME/hdf5/lib -L$HOME/hdf5/lib -I$HOME/hdf5/include ${spath}/analyse_static_structure_log.cpp ${upath}/read_write.cpp -lhdf5 -fopenmp -o analyse_static_struct
// ./analyse_static_struct out.h5 ${path}/Static_struct.txt

//////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <vector>
#include "omp.h"
#include "../Utility/read_write.hpp"
#include "../Utility/basic.hpp"
#include "../Utility/sim_info.hpp"

#define pi M_PI

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////

void compute_static_structure (double **x, double **y, double *Sk, 
			       double *kxvec, double *kyvec, double *kvec,
			       int Nmax, int ndata, int natoms, int nsteps, double l, int njump) {
  /* calculate static structure per timestep with a running average */

  // set preliminary data
  
  vector<int> kcount(Nmax, 0);
  
  for (int step = 0; step < nsteps; step += njump) {
    
    cout << "step / nsteps: " << step << " / " << nsteps << endl;
    
    for (int kx = 0; kx < ndata; kx++) {
	double kxt = kxvec[kx]; 

	for (int ky = 0; ky < ndata; ky++) {
	  double costerm = 0.;
	  double sinterm = 0.;
	  double kyt = kyvec[ky];

	  omp_set_num_threads(4);
	  #pragma omp parallel for reduction(+:costerm,sinterm)
	  for (int j = 0; j < natoms; j++) {
	    double dotp = kxt*x[step][j] + kyt*y[step][j];
	    costerm += cos(dotp);
	    sinterm += sin(dotp);
    
	  } // particle loop

	  int knorm = (int)sqrt(kxt*kxt + kyt*kyt);
// 	  cout << "kx = " << kx << "\t ky = " << ky << endl;
// 	  cout << "knorm = " << knorm << "\t lx = " << l << "\t Nmax = " << Nmax << endl;
	  kcount[knorm]++;
	  Sk[knorm] += costerm*costerm + sinterm*sinterm;

	} // ky loop
		
    } // kx loop 
  
  } // timesteps
  
  // perform the normalization
  
  double totalData = nsteps/njump;
  for (int j = 0; j < Nmax; j++)  { 
    Sk[j] /= (natoms*kcount[j]); 
  }  
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

int main (int argc, char *argv[]) {
 
  // get the file name by parsing and load simulation data
  
  string filename = argv[1];
  cout << "Calculating static structure factor of the following file: \n" << filename << endl;
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
  
  // get the image positions in the central box
  
  get_img_pos(x, y, sim.nsteps, sim.ncells, sim.lx, sim.ly);
   
  // allocate the arrays and set preliminary information

  double delk = 2.*pi/sim.lx;
  int njump = 100;
  int ndata = (int)sim.lx;
  double *kxvec = new double[ndata];
  double *kyvec = new double[ndata];  
  for (int j = 0; j < ndata; j++) {
    kxvec[j] = delk*j;  kyvec[j] = delk*j; 
  }
  
  double mink = 0.;
  double maxk = ceil(sqrt(kxvec[ndata-1]*kxvec[ndata-1] + kyvec[ndata-1]*kyvec[ndata-1]));
  int Nmax = (int)maxk;
  double *Sk = new double[Nmax];
  double *kvec = new double[Nmax];
  for (int j = 0; j < Nmax; j++) {
    Sk[j] = 0.;  kvec[j] = j;
  }

  // calculate the static structure factor
  
  compute_static_structure(x, y, Sk, kxvec, kyvec, kvec, Nmax, ndata, sim.ncells, sim.nsteps, sim.lx, njump);  
  
  // write the computed data
  
  string outfilepath = argv[2];
  cout << "Writing static structure factor to the following file: \n" << outfilepath << endl;  
  for (int j = 0; j < Nmax; j++) {
    cout << kvec[j] << "\t" << Sk[j] << endl;
  }
  //write_2d_analysis_data(kvec, Sk, Nmax, outfilepath);
  
  // deallocate the arrays
  
  delete [] Sk;
  delete [] kvec;
  delete [] kxvec;
  delete [] kyvec;
  for (int i = 0; i < sim.nsteps; i++) {
    delete [] x[i];
    delete [] y[i];
  }
  delete [] x;
  delete [] y;
  
  return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
