
#pragma once

#include <string>
#include <iostream>
#include <vector>

class SimInfo {
  /* container for general simulation data */
  public:
    
    SimInfo(std::string);
    int nsteps, nbeads, nsamp, ncells;
    double lx, ly, dt, eps, rho, fp, areak, kappa, bl, sigma;
    std::vector<int> nbpc;
  
};

SimInfo::SimInfo(std::string filename) {
  
  // read in general simulation data

  nsteps = nbeads = nsamp = ncells = 0;
  lx = ly = dt = eps = rho = fp = areak = kappa = bl = sigma = 0.;

  read_sim_data(filename, nsteps, nbeads, nsamp, ncells, lx, ly, dt, eps, rho, fp, areak, kappa, bl, sigma);

  int nbpc_buffer[ncells];
  for (int i = 0; i < ncells; i++) nbpc_buffer[i] = 0;
  read_integer_array(filename, "/cells/nbpc", nbpc_buffer);  
  for (int i = 0; i < ncells; i++) nbpc.push_back(nbpc_buffer[i]);

  return;
}