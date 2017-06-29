
/* basic calculations */

// COMPILATION AND RUN COMMANDS:
// -
// -

//////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <iostream>
#include <vector>
#include <cmath>

//////////////////////////////////////////////////////////////////////////////////////////////////////////

inline int ifloor (double doubval) {
    return ( doubval >= 0. ? (int)(doubval) : ((int)doubval)-1 );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

inline long double dnearbyint (double doubval) {
    return ( doubval >= 0. ? (long)(doubval + 0.5) : (long)(doubval - 0.5) );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

inline int inearbyint (double doubval) {
    return ( doubval >= 0. ? (int)(doubval + 0.5) : (int)(doubval - 0.5) );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class T>
inline void SWAP (T &a, T &b) {
    T dum=a;
    a=b;
    b=dum;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

inline double get_single_img_pos (double x, double l) {
  /* calculate the image position in the central box */
  
  return x-ifloor(x/l)*l;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void get_img_pos (double **x, double **y, int nsteps, int natoms, double lx, double ly) {
  /* calculate the image positions in the central box */
  
  for (int i = 0; i < nsteps; i++) {
    for (int j = 0; j < natoms; j++) {
      
      x[i][j] -= ifloor(x[i][j]/lx)*lx;
      y[i][j] -= ifloor(y[i][j]/ly)*ly;
      
    }
  }
 
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

inline double get_min_img_dist (double dx, double lx) {
  /* calculate the minimum image distance */
  
  return dx-lx*dnearbyint(dx/lx);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class T>
inline void deallocate_1d_array (T *&arr) {
  /* deallocate 1D array */
  
  delete [] arr;
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class T>
inline void deallocate_2d_array (T **&arr, int ndata) {
  /* deallocate 2D array */
  
  for (int i = 0; i < ndata; i++) {
    delete [] arr[i];
  }
  delete [] arr;
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

inline int get_bin_number (double x, double w, int nbins) {
  /* calculate the bin of the coordinate assumming it is image */
  
  int bin = static_cast<int>(x/w);
  if (bin < 0) {
    bin += nbins;
  }
  else if (bin >= nbins) {
    bin -= nbins;
  }
  
  return bin;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
