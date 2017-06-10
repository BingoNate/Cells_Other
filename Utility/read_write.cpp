
/* helper function in data reading and writing */

// COMPILATION AND RUN COMMANDS:
// g++ -Wl,-rpath=$HOME/hdf5/lib -L$HOME/hdf5/lib -I$HOME/hdf5/include ${spath}/read_write.cpp -lhdf5 -o read_write
// -

//////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "read_write.hpp"

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////

int read_integer_data (hid_t file, string path_in_file, int *buffer) {
  /* wrapper to read integer data from hdf5 file --note that buffer needs to be an array of size 1 for single entries-- */
  
  const char *fl = path_in_file.c_str();	
  hid_t dataset = H5Dopen(file, fl, H5P_DEFAULT);
  herr_t status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);

  H5Dclose(dataset);

  return buffer[0];
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

double read_double_data (hid_t file, string path_in_file, double *buffer) {
  /* wrapper to read double data from hdf5 file --note that buffer needs to be an array of size 1 for single entries-- */
   
  const char *fl = path_in_file.c_str(); 
  hid_t dataset = H5Dopen(file, fl, H5P_DEFAULT);
  herr_t status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);

  H5Dclose(dataset);

  return buffer[0];
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void read_integer_array (string filename, string path_in_file, int *buffer) {
  /* wrapper to read integer array data from hdf5 file --note that buffer needs to be the array size-- */

  const char *fl1 = filename.c_str();
  const char *fl2 = path_in_file.c_str();

  // open the file pointer
  
  hid_t file = H5Fopen(fl1, H5F_ACC_RDONLY, H5P_DEFAULT);
  
  // read the array
  
  hid_t dataset = H5Dopen(file, fl2, H5P_DEFAULT);
  herr_t status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);

  H5Dclose(dataset);
  H5Fclose(file);

  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void read_double_array (string filename, string path_in_file, double *buffer) {
  /* wrapper to read double array data from hdf5 file --note that buffer needs to be the array size-- */

  const char *fl1 = filename.c_str();
  const char *fl2 = path_in_file.c_str();

  // open the file pointer
  
  hid_t file = H5Fopen(fl1, H5F_ACC_RDONLY, H5P_DEFAULT);
  
  // read the array
  
  hid_t dataset = H5Dopen(file, fl2, H5P_DEFAULT);
  herr_t status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);

  H5Dclose(dataset);
  H5Fclose(file);

  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void read_sim_data (string filename, int &nsteps, int &nbeads, int &nsamp, int &ncells, double &lx, double &ly, double &dt, double &eps, double &rho, double &fp, double &areak, double &kappa, double &bl, double &sigma) {
  /* read general simulation data in hdf5 format */
 
  const char *fl1 = filename.c_str();
 
  // open the file pointer
  
  hid_t file = H5Fopen(fl1, H5F_ACC_RDONLY, H5P_DEFAULT);

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
  kappa = read_double_data(file, "/param/kappa", d_buffer);
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

void read_all_pos_data(string filename, double **x, double **y, int nsteps, int natoms, string datapath) {
  /* read position data in hdf5 format all at once */
 
  const char *fl1 = filename.c_str();
  const char *fl2 = datapath.c_str();

  // open the file pointer
  
  hid_t file = H5Fopen(fl1, H5F_ACC_RDONLY, H5P_DEFAULT);

  // get the dataset in the file
  /* DATASET (H5D) is the raw data (either singular or in arrays) */
  
  hid_t dataset = H5Dopen(file, fl2, H5P_DEFAULT);
  
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

void read_single_pol_data(int step, hid_t dataset, hid_t dataspace, double **pol, int natoms) {
  /* read the polarity data at a single timestep */
 
  // read the data in the x direction
  
  // define the hyperslab in the dataset
  /* we are gonna reduce the data in the dataset from 3 dimensions to two 2 dimensional arrays */
    
  hsize_t offset[2];
  offset[0] = step; offset[1] = 0; 
  
  hsize_t count[2];
  count[0] = 1; count[1] = natoms;
  
  // select a 1D hyperslab from the original 2D dataset
  
  herr_t status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);

  // define memory dataspace
  
  hsize_t dimsm[2];	// dimensions and sizes in each dimension
  dimsm[0] = 1; dimsm[1] = natoms; 
  hid_t memspace = H5Screate_simple(2, dimsm, NULL);
  
  // define memory hyperslab
  
  hsize_t offset_out[2];
  offset_out[0] = 0; offset_out[1] = 0; 
  
  hsize_t count_out[2];
  count_out[0] = 1; count_out[1] = natoms; 
  
  status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
  
  // read data from hyperslab in the file into the hyperslab in memory 
  
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, pol[step]); 
      
  H5Sclose(memspace);
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void read_all_pol_data(string filename, double **pol, int nsteps, int natoms, string datapath) {
  /* read polarity data in hdf5 format all at once */
 
  const char *fl1 = filename.c_str();
  const char *fl2 = datapath.c_str();

  // open the file pointer
  
  hid_t file = H5Fopen(fl1, H5F_ACC_RDONLY, H5P_DEFAULT);

  // get the dataset in the file
  /* DATASET (H5D) is the raw data (either singular or in arrays) */
  
  hid_t dataset = H5Dopen(file, fl2, H5P_DEFAULT);
  
  // get dataspace of the selected dataset
  /* DATASPACE (H5S) describes the number of dimensions and the size of the dataset in those dimension */
  
  hid_t dataspace = H5Dget_space(dataset);
 
  /* READ THE DATA
  load the data into the memory step by step 
  register the data at each step to the array
  */
  
  for (int step = 0; step < nsteps; step++) {
    read_single_pol_data(step, dataset, dataspace, pol, natoms);
  }

  H5Sclose(dataspace);
  H5Dclose(dataset);
  H5Fclose(file);
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void write_single_analysis_data(double data, string outpath) {
  /* write analysis data with single parameter to the outfile */
 
  const char *outfl = outpath.c_str();
  ofstream fl(outfl);
  fl << data << endl;
  fl.close();
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void write_1d_analysis_data(double *x, int ndata, string outpath) {
  /* write the 1d analysis data to the outfile */
  
  const char *outfl = outpath.c_str();
  ofstream fl(outfl);
  for (int j = 0; j < ndata; j++) {
    fl << j << "\t\t" << x[j] << endl;
  }
  fl.close();
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void write_2d_analysis_data(double *x, double *y, double ndata, string outpath) {
  /* write the 2d analysis data to the outfile */
  
  const char *outfl = outpath.c_str();
  ofstream fl(outfl);
  for (int j = 0; j < ndata; j++) {
    fl << x[j] << "\t\t" << y[j] << endl;
  }
  fl.close();
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void write_multid_analysis_data(const vector<double> &w, 
			    const vector<double> &vx, const vector<double> &vy, 
			    const int ndata, const int nsteps, string outpath) {
  /* write multid analysis data to the outfile */

  const char *outfl = outpath.c_str();
  ofstream fl(outfl);
  for (int step = 0; step < nsteps; step++) {
    for (int i = 0; i < ndata; i++) {
      for (int j = 0; j < ndata; j++) {
	int idx = step*ndata*ndata + ndata*i + j;
	fl << step << "\t" << i << "\t" << j << "\t" << w[idx] << "\t" << vx[idx] << "\t" << vy[idx] << "\t" << sqrt(vx[idx]*vx[idx]+vy[idx]*vy[idx]) << endl;
      }
    }
  }	
    
  fl.close();
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
