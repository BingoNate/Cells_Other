
#pragma once

// Define boundary conditions
#define PERIODIC

// Define starting from a backup
#define INIT

// Define parallel programming
//#define PARALLEL
#ifdef PARALLEL
#include "omp.h"
#endif

// Define debug
//#define DEB

// Include packages
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <iomanip>

// Include dependent header files
#include "basic.h"

// Pi number
#define pi M_PI


//+++++++++++++++++++
//++ Variable declarations
//+++++++++++++++++++

//+++++++++++++++++++
//++ Core variables
//+++++++++++++++++++

// Simulation variables
extern const double dt;                 // Time step size in MD simulation
extern const int sampleFreq;            // Sampling frequency in MD simulation

// Cell variables
extern const int M;                     // Number of cells
extern int L;                           // Maximum total number of particles in the simulation
extern int* Npart;                      // Number of particles on each cell
extern const int N_avg;                 // Average number of particles on each cell
extern const double B;                  // Bond length between particles of cells
extern double* R;                       // Radii of cells
extern const double R_avg;              // Mean radius of a cell
extern const double polyd;              // Polydispersity in radius of cells (percentage)
extern double* A;                       // Initial area of cells
extern const double mass;               // Mass -- ( Viscous time scale is mass/gamma )
extern double* theta;                   // Angle between the particles in a cell
extern double* phi;                     // Polarity of the cells


// Initial grid
extern const double interDist_x;    // Interparticle distance
extern const double interDist_y;
extern const double initOffset_x;   // Offset value for the initial positions
extern const double initOffset_y;
extern const int gridDim_x;         // Number of particles on the grid
extern const int gridDim_y;
extern const double gridSize_x;     // Sizes of the grid
extern const double gridSize_y;

// Boundary conditions
extern const double packing;
extern const double boxSize_x;
extern const double boxSize_y;

//+++++++++++++++++++
//++ Force variables
//+++++++++++++++++++

// Spring force
extern const double kcons;          // Spring constant

// Bending force
extern const double Kbend;          // Bending rigidity

// Area constraint force
extern const double A0_;            // 'Target'/Equilibrium area
extern double* A0;
extern const double Acons;          // Membrane elasticity

// Drag force
extern const double gam;            // Translational friction coefficient

// Diffusion
extern const double kB;             // Boltzmann constant
extern const double T;              // Temperature
extern const double Dt;             // Translational diffusion coefficient (kB*T/gam)

// Polarity rotation
extern const double Dr;             // Rotational diffusion coefficient

// Motility force
extern double* Fm;                  // Motility force constant matrix
extern double* Fmsin;
extern double* Fmcos;
extern const double Fmc;            // Motility force constant

// Lennard Jones interaction
extern const double eps;            // Depth (minimum) of the LJ potential
extern const double totalStep;      // Total number of MD steps
extern const double sig;            // Distance where the LJ potential becomes zero
extern const double Rcut;           // Cutoff distance for the LJ potential


//++++++++++++++++
//++ Verlet listing
//++++++++++++++++

// Verlet listing
extern int* NV;                     // Number of particles in the Verlet list
extern const double Rver;           // Verlet radius
extern const double skin;           // Skin for updating the Verlet list
extern int vlistsize;
extern int* verList;                // Verlet list of particle
extern int* vlist_offset;
extern int* Nhead;                  // Mark particles with cells they belong

// Boxing
extern int kx;                      // Box index
extern int ky;
extern int kbox;                    // Serial box index
extern int nbox_x;                  // Neighbor box index
extern int nbox_y;
extern const double Rbox_x;         // Box size for boxing
extern const double Rbox_y;
extern const int KX;                // Number of boxes
extern const int KY;
extern const int K2;
extern int* hoc;                    // Head of chain
extern int* link_list;              // Linked list

//+++++++++++++++++
//++ Other variables
//+++++++++++++++++

// Arrays
extern double* xc;                  // x center of the cell
extern double* yc;                  // y center of the cell
extern double* x;                   // x position of the cell
extern double* y;                   // y position of the cell
extern double* vx;                  // x velocity of the cell
extern double* vy;                  // y velocity of the cel
extern double* Fx;                  // x force on particle
extern double* Fy;                  // y force on particle
extern double* xb;                  // x position in the central box
extern double* yb;                  // y position in the central box
extern double* xv;                  // Old position of the particles in x
extern double* yv;                  // Old position of the particles in y

// Useful calculation variables
extern const double sampleDt;
extern const int totalSample;
extern const double Rcut2;
extern const double Rver2;
extern const double skin2;


