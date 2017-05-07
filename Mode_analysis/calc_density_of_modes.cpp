/* calculate density of modes/states */

// COMPILATION AND RUN COMMANDS:
// g++ -O3 -lm -Wl,-rpath=$HOME/hdf5/lib -L$HOME/hdf5/lib -I$HOME/hdf5/include ${spath}/calc_density_of_modes.cpp ${upath}/read_write.cpp -lhdf5 -o calc_density_of_modes
// ./calc_density_of_modes out.h5 ${path}/Density_of_modes.txt

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
#include </usr/users/iff_th2/duman/Eigen/Eigen/Dense>
#include </usr/users/iff_th2/duman/Eigen/Eigen/Eigenvalues>
#include "../Utility/read_write.hpp"
#include "../Utility/basic.hpp"

#define pi M_PI

using namespace std;
using namespace Eigen;

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void calc_avg_pos (vector<double> &xavg, vector<double> &yavg, 
		   double **x, double **y, 
		   const int ncells, const int nsteps) {
  /* calculate average position of each cell in time */
  
  cout << "Calculating the average position of each cell in time" << endl; 
  
  for (int j = 0; j < ncells; j++) {
    
    for (int step = 0; step < nsteps; step++) {
      
      xavg[j] += x[step][j];
      yavg[j] += y[step][j];
      
    }	 // cells loop
    
    xavg[j] /= nsteps;
    yavg[j] /= nsteps;
    
  } 	// timesteps loop
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void calc_dynamical_matrix (MatrixXd &Dab, double **x, double **y,
			    const vector<double> &xavg, const vector<double> &yavg,
			    const int ncells, const int nsteps) {
  /* calculate the averaged dynamical matrix */

  cout << "Calculating and time averaging the dynamical matrix" << endl; 

  for (int step = 0; step < nsteps; step++) {
    
    int i_even = 0;
    int i_odd = 0;
    
    for (int i = 0; i < 2*ncells; i++) {
      
      double row = 0.;
      if (i % 2 == 0) {
        row = x[step][i_even] - xavg[i_even];
	i_even++;
      }
      else {
        row = y[step][i_odd] - yavg[i_odd];
	i_odd++;
      }
      
      int j_even = 0;
      int j_odd = 0;
      
      for (int j = 0; j < 2*ncells; j++) {
	
        double column = 0.;
        if (j % 2 == 0) {
          column = x[step][j_even] - xavg[j_even];
	  j_even++;
        }
        else {
          column = y[step][j_odd] - yavg[j_odd];
	  j_odd++;
        }
        
        Dab(i, j) += row*column;
	
      }	  // second cells loop (column)
    }  	  // first cells loop (row)
    
  }       // timesteps loop
  
  for (int i = 0; i < 2*ncells; i++) {
    for (int j = 0; j < 2*ncells; j++) {
      Dab(i, j) /= nsteps;
    }
  }
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void diag_dynamical_matrix(MatrixXd &Dab, double *evals, const int ncells) {
  /* diagonalize the dynamical matrix to calculate the eigenvalues */
  
  cout << "Diagonalizing the dynamical matrix" << endl;
  
  // diagonalize the dynamical matrix

  SelfAdjointEigenSolver<MatrixXd> es(Dab);
  if (es.info() != Success) abort();
  cout << es.eigenvalues() << endl;
  
  for (int j = 0; j < 2*ncells; j++)
    evals[j] = es.eigenvalues()[j];
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

int main (int argc, char *argv[]) {

  // get the file name by parsing
  
  string filename = argv[1];
  cout << "Calculating density of modes of the following file: \n" << filename << endl;

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
  
  // calculate the average position of each cell
  
  vector<double> xavg(ncells, 0.);
  vector<double> yavg(ncells, 0.);
  calc_avg_pos(xavg, yavg, x, y, ncells, nsteps);

  // set variables related to the analysis
  
  MatrixXd Dab(2*ncells, 2*ncells);
  Dab << MatrixXd::Zero(2*ncells, 2*ncells);
  double *evals = new double[2*ncells];
  for (int j = 0; j < 2*ncells; j++)
    evals[j] = 0.;
  
  // calculate the dynamical matrix
  
  calc_dynamical_matrix(Dab, x, y, xavg, yavg, ncells, nsteps);
  
  // diagonalize the dynamical matrix
  
  diag_dynamical_matrix(Dab, evals, ncells);
		    
  // write the computed data
  
  string outfilepath = argv[2];
  cout << "Writing density of modes to the following file: \n" << outfilepath << endl;
  write_1d_analysis_data(evals, 2*ncells, outfilepath);
  
  // deallocate the arrays
  
  delete [] evals;
  for (int i = 0; i < nsteps; i++) {
    delete [] x[i];
    delete [] y[i];
  }
  delete [] x;
  delete [] y;
  
  return 0;
}  
