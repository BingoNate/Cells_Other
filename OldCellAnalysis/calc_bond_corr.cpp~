

/* calculate bond correlation function */

// COMPILATION AND RUN COMMANDS:
// g++ -Wl,-rpath=$HOME/hdf5/lib -L$HOME/hdf5/lib -I$HOME/hdf5/include ${spath}/calculate_inter_scattering.cpp -lhdf5 -fopenmp -o calc_inter_scattering
// ./calc_inter_scattering out.h5 ${path}/Inter_scattering.txt

//////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <cmath>
#include "/usr/users/iff_th2/duman/hdf5/include/hdf5.h"
//#include "/homec/jiff26/jiff2610/hdf/include/hdf5.h"
#include "omp.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <map>
#include <algorithm>
#include <boost/multi_array.hpp>
#include "variables.h"
#include "basic.h"
#include "pointbox.h"
#include "kdtree.h"

#define pi M_PI

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////

int read_integer_data (hid_t file, char *path_in_file, int *buffer) {
  /* wrapper to read integer data from hdf5 file --note that buffer needs to be an array of size 1 for single entries-- */
    
  hid_t dataset = H5Dopen(file, path_in_file, H5P_DEFAULT);
  herr_t status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);

  H5Dclose(dataset);

  return buffer[0];
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

double read_double_data (hid_t file, char *path_in_file, double *buffer) {
  /* wrapper to read double data from hdf5 file --note that buffer needs to be an array of size 1 for single entries-- */
    
  hid_t dataset = H5Dopen(file, path_in_file, H5P_DEFAULT);
  herr_t status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);

  H5Dclose(dataset);

  return buffer[0];
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void read_integer_array (char *filename, char *path_in_file, int *buffer) {
  /* wrapper to read integer array data from hdf5 file --note that buffer needs to be the array size-- */

  // open the file pointer
  
  hid_t file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  
  // read the array
  
  hid_t dataset = H5Dopen(file, path_in_file, H5P_DEFAULT);
  herr_t status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);

  H5Dclose(dataset);
  H5Fclose(file);

  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void read_double_array (char *filename, char *path_in_file, double *buffer) {
  /* wrapper to read double array data from hdf5 file --note that buffer needs to be the array size-- */

  // open the file pointer
  
  hid_t file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  
  // read the array
  
  hid_t dataset = H5Dopen(file, path_in_file, H5P_DEFAULT);
  herr_t status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);

  H5Dclose(dataset);
  H5Fclose(file);

  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void read_sim_data (char *filename, int &nsteps, int &nbeads, int &nsamp, int &ncells, double &lx, double &ly, double &dt, double &eps, double &rho, double &fp, double &areak, double &bl, double &sigma) {
  /* read general simulation data in hdf5 format */
  
  // open the file pointer
  
  hid_t file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

  // create buffers to get single data --note that buffer needs to be an array of size 1 for single entries--
  
  int i_buffer[1];
  i_buffer[0] = 0;
  double d_buffer[1];
  d_buffer[0] = 0.;
    
  // read in the box info

  lx = read_double_data(file, "/info/box/x", d_buffer);
  ly = read_double_data(file, "/info/box/y", d_buffer);
  
  // read in the general simulation info

  dt = read_double_data(file, "/info/dt", d_buffer);
  nsteps = read_integer_data(file, "/info/nsteps", i_buffer);
  nbeads = read_integer_data(file, "/info/nbeads", i_buffer);
  nsamp = read_integer_data(file, "/info/nsamp", i_buffer);
  ncells = read_integer_data(file, "/info/ncells", i_buffer);  

  // read in the simulation parameters
  
  eps = read_double_data(file, "/param/eps", d_buffer);
  rho = read_double_data(file, "/param/rho", d_buffer);
  fp = read_double_data(file, "/param/fp", d_buffer);
  areak = read_double_data(file, "/param/areak", d_buffer);
  bl = read_double_data(file, "/param/bl", d_buffer);
  sigma = read_double_data(file, "/param/sigma", d_buffer);

  H5Fclose(file);
  
  // normalize some of the variables
  
  dt *= nsamp;
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void read_single_pos_data(int step, hid_t dataset, hid_t dataspace, double **x, double **y, int natoms) {
  /* read the position data at a single timestep */
  
  // read the data in the x direction
  
  // define the hyperslab in the dataset
  /* we are gonna reduce the data in the dataset from 3 dimensions to two 2 dimensional arrays */
    
  hsize_t offset[3];
  offset[0] = step; offset[1] = 0; offset[2] = 0;
  
  hsize_t count[3];
  count[0] = 1; count[1] = 1; count[2] = natoms;
  
  // select a 2D hyperslab from the original 3D dataset
  
  herr_t status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);

  // define memory dataspace
  
  hsize_t dimsm[3];	// dimensions and sizes in each dimension
  dimsm[0] = 1; dimsm[1] = 1; dimsm[2] = natoms;
  hid_t memspace = H5Screate_simple(3, dimsm, NULL);
  
  // define memory hyperslab
  
  hsize_t offset_out[3];
  offset_out[0] = 0; offset_out[1] = 0; offset_out[2] = 0;
  
  hsize_t count_out[3];
  count_out[0] = 1; count_out[1] = 1; count_out[2] = natoms;
  
  status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
  
  // read data from hyperslab in the file into the hyperslab in memory 
  
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, x[step]); 
      
  H5Sclose(memspace);
  
  // read the data in the y direction
  
  offset[0] = step; offset[1] = 1; offset[2] = 0;
  count[0] = 1; count[1] = 1; count[2] = natoms;
  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);   
  memspace = H5Screate_simple(3, dimsm, NULL);
  status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, y[step]);   
  H5Sclose(memspace);
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void read_all_pos_data(char *filename, double **x, double **y, int nsteps, int natoms) {
  /* read position data in hdf5 format all at once */
  
  // open the file pointer
  
  hid_t file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

  // get the dataset in the file
  /* DATASET (H5D) is the raw data (either singular or in arrays) */
  
  hid_t dataset = H5Dopen(file, "/beads/xu", H5P_DEFAULT);
  
  // get dataspace of the selected dataset
  /* DATASPACE (H5S) describes the number of dimensions and the size of the dataset in those dimension */
  
  hid_t dataspace = H5Dget_space(dataset);
 
  /* READ THE DATA
  load the data into the memory step by step 
  register the data at each step to the array
  */
  
  for (int step = 0; step < nsteps; step++) {
    read_single_pos_data(step, dataset, dataspace, x, y, natoms);
  }

  H5Sclose(dataspace);
  H5Dclose(dataset);
  H5Fclose(file);
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void write_data(double *delay, double *Fs, int ndelay, char *outpath) {
  /* write the analysis data to the outfile */
  
  ofstream fl(outpath);
  for (int j = 0; j < ndelay; j++) {
    fl << delay[j] << "\t\t" << Fs[j] << endl;
  }
  fl.close();
  
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

int main( int argc, char *argv[] ){
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++ Get the data
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    // Choose the data to be used
    const int nsamp = samplePostData;       // Sampling interval
    const double dtSamp = nsamp*dt;         // Time units between two data points
    
    // Instantiate the data structure
    vector<Cell> cells( M );
    
    // Set filenames
    char * x_input_file = argv[1];        // Filename for the x data
    char * y_input_file = argv[2];        // Filename for the y data
    
    // Read the data to the cells
    int sample_cnt = -1;
    int sample_data = 1;
    char getX = 'x';
    readData( cells, x_input_file, getX, sample_cnt, sample_data );
    sample_cnt = -1;
    char getY = 'y';
    readData( cells, y_input_file, getY, sample_cnt, sample_data );
    
    // Set general simulation variables
    Data simData;
    simData.setNumStep( cells[0].xpos.size() );
    simData.setNumDelay( sqrt( cells[0].xpos.size() ) );
    simData.setNumTotalDelay();
    
    const double T = simData.getNumStep();              // Total time
    const double D = simData.getNumDelay();             // Last delay time
    const double TD = simData.getNumTotalDelay();       // Total time - last delay time
    
    // Set the box
    SimBox box;
    box.setWidth( boxSize_x );
    box.setHeight( boxSize_y );
    
    const double Lx = box.getWidth();
    const double Ly = box.getHeight();
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++ Do the analysis
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++ Bond orientation
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    
    #ifdef BOND_ORIENT
    
    const int MM = 9*M;
    vector<Point<2> > cell_vec(MM);
    
    for(int m = 0; m < M; m++){
        
        cells[m].getImag(box, 0);
        
        cell_vec[m].x[0] = cells[m].xi;
        cell_vec[m].x[1] = cells[m].yi;
        
        cell_vec[M+m].x[0] = cells[m].xi - Lx;
        cell_vec[M+m].x[1] = cells[m].yi + Ly;
        
        cell_vec[2*M+m].x[0] = cells[m].xi;
        cell_vec[2*M+m].x[1] = cells[m].yi + Ly;
        
        cell_vec[3*M+m].x[0] = cells[m].xi + Lx;
        cell_vec[3*M+m].x[1] = cells[m].yi + Ly;
        
        cell_vec[4*M+m].x[0] = cells[m].xi + Lx;
        cell_vec[4*M+m].x[1] = cells[m].yi;
        
        cell_vec[5*M+m].x[0] = cells[m].xi + Lx;
        cell_vec[5*M+m].x[1] = cells[m].yi - Ly;
        
        cell_vec[6*M+m].x[0] = cells[m].xi;
        cell_vec[6*M+m].x[1] = cells[m].yi - Ly;
        
        cell_vec[7*M+m].x[0] = cells[m].xi - Lx;
        cell_vec[7*M+m].x[1] = cells[m].yi - Ly;
        
        cell_vec[8*M+m].x[0] = cells[m].xi - Lx;
        cell_vec[8*M+m].x[1] = cells[m].yi;
        
    }
    
    KDtree<2> kdtree(cell_vec);
    
    ofstream local_bond_out("local_bond.txt");
    vector<double> theta_real(M);
    vector<double> theta_img(M);
    vector<double> psi6(M);
    for(int m = 0; m < M; m++){
        theta_real[m] = 0;
        theta_img[m] = 0;
        psi6[m] = 0;
    }
    double glob_psi6 = 0;
    double glob_psi6_avg = 0;
    double psi6_bar_temp = 0;
    vector<double> psi6_bar;
    
    int n = 4;              // Number of neighbors
    
    int t_offset = 0;
    int t_offset_2 = T;
    int t_offset_total = t_offset_2 - t_offset;
    
    for( int i = t_offset; i < t_offset_2; i++){        // Run over time
        
        // Create the points that will make up the tree
        const int MM = 9*M;
        vector<Point<2> > cell_vec(MM);
        
        // Consider all the image points as well
        for(int m = 0; m < M; m++){
            
            cells[m].getImag(box, i);
            
            cell_vec[m].x[0] = cells[m].xi;
            cell_vec[m].x[1] = cells[m].yi;
            
            cell_vec[M+m].x[0] = cells[m].xi - Lx;
            cell_vec[M+m].x[1] = cells[m].yi + Ly;
            
            cell_vec[2*M+m].x[0] = cells[m].xi;
            cell_vec[2*M+m].x[1] = cells[m].yi + Ly;
            
            cell_vec[3*M+m].x[0] = cells[m].xi + Lx;
            cell_vec[3*M+m].x[1] = cells[m].yi + Ly;
            
            cell_vec[4*M+m].x[0] = cells[m].xi + Lx;
            cell_vec[4*M+m].x[1] = cells[m].yi;
            
            cell_vec[5*M+m].x[0] = cells[m].xi + Lx;
            cell_vec[5*M+m].x[1] = cells[m].yi - Ly;
            
            cell_vec[6*M+m].x[0] = cells[m].xi;
            cell_vec[6*M+m].x[1] = cells[m].yi - Ly;
            
            cell_vec[7*M+m].x[0] = cells[m].xi - Lx;
            cell_vec[7*M+m].x[1] = cells[m].yi - Ly;
            
            cell_vec[8*M+m].x[0] = cells[m].xi - Lx;
            cell_vec[8*M+m].x[1] = cells[m].yi;
            
        }
        
        // Build the tree from the points
        KDtree<2> kdtree(cell_vec);
        
        for(int m = 0; m < M; m++){                 // Run over cells
            
            int *nn = new int[n];
            double *dn = new double[n];
            kdtree.nnearest(m,nn,dn,n);
            
            for(int j = 0; j < n; j++){             // Run over neighbors
                double dy = dist_xy(cell_vec[nn[j]], cell_vec[m], 1);
                double dx = dist_xy(cell_vec[nn[j]], cell_vec[m], 0);
                double theta = 4*atan2( dy, dx );
                double cos_theta = cos( theta );
                double sin_theta = sin( theta );
                theta_real[m] += cos_theta;
                theta_img[m] += sin_theta;
            }
            
        }
        
        for(int m = 0; m < M; m++){
            double tmp_value = sqrt( theta_real[m]*theta_real[m] + theta_img[m]*theta_img[m] )/n;
            //double tmp_value = theta_real[m]/n;
            psi6[m] += tmp_value;
            glob_psi6 += tmp_value;
            psi6_bar_temp += tmp_value;
            theta_real[m] = 0;
            theta_img[m] = 0;
        }
        
        glob_psi6_avg += glob_psi6/M;
        glob_psi6 = 0;
        psi6_bar.push_back( psi6_bar_temp/M );
        psi6_bar_temp = 0;
        
    } // for i (time)
    
    for(int m = 0; m < M; m++){
        psi6[m] /= t_offset_total;
        local_bond_out << psi6[m] << "\n";
    }
    
    ofstream glob_bond_out;
    glob_bond_out.open( "glob_bond.txt", ios::out | ios::app );
    glob_bond_out << dens << "\t" << Fmc << "\t" << eps << "\t" << glob_psi6_avg/t_offset_total << "\n";
    glob_bond_out.close();
    
    double pre_sq = 0;
    double pre_avg = 0;
    for( int i = t_offset; i < t_offset_2; i++ ){
        pre_sq += psi6_bar[i]*psi6_bar[i];
        pre_avg += psi6_bar[i];
    }
    pre_sq /= t_offset_total;
    pre_avg /= t_offset_total;
    ofstream susc_out;
    susc_out.open( "susceptibility.txt", ios::out | ios::app );
    susc_out << dens << "\t" << Fmc << "\t" << eps << "\t" << M*(pre_sq - pre_avg*pre_avg) << "\n";
    susc_out.close();
    
    #endif // End of bond orientation
    
    #ifdef BOND_CORR
    
    double max_dist = 0.0;
    int sug_max_dist = (int)fmax( boxSize_x, boxSize_y );  // Suggest a sane approximation for maximum distance
    vector<double> cvv(sug_max_dist, 0.0);
    vector<double> cvv_n(sug_max_dist, 0.0);

    for( int i = t_offset; i < t_offset_2; i++ ){           // Run over all the timeseries
        
        for( int m1 = 0; m1 < M-1; m1++ ){                  // Run over cells to compute the distances
            for( int m2 = m1+1; m2 < M; m2++ ){
                
                // Calculate the distances between two cells with minimum image convention
                double dist_X = cells[m2].xpos[i] - cells[m1].xpos[i];
                dist_X = dist_X - Lx*dnearbyint( dist_X/Lx );
                
                double dist_Y = cells[m2].ypos[i] - cells[m1].ypos[i];
                dist_Y = dist_Y - Ly*dnearbyint( dist_Y/Ly );
                
                double dist_R = sqrt( dist_X*dist_X + dist_Y*dist_Y );
                
                // Get the bin corresponding to the distance
                int bin = inearbyint( dist_R );
                
                // Calculate the bond correlation function
                cvv[bin] += psi6[m2]*psi6[m1];
                cvv_n[bin] += 1;
                
                // Keep track of maximum distance between the cells
                if( dist_R > max_dist ){
                    max_dist = dist_R;
                }
                
            } // for m2
        } // for m1
        
    } // for running over time
    
    // Get the last bin
    int total_bins = (int)max_dist;    // Total number of bins
    
    ofstream bond_corr_out("bond_corr.txt");
    // Write the bond correlation function
    for( int i = 0; i < total_bins; i++ ){
        bond_corr_out << i/(2*R_avg) << "\t" << cvv[i]/cvv_n[i] << "\n";
    }
    
    #endif // End of bond correlation
    
} // End of the program
