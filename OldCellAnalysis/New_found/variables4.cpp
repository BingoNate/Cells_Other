
//+++++++++++++++++++
//++ Variable definitions
//+++++++++++++++++++

#include "variables.h"

using namespace std;

#define pi M_PI


// Simulation variables
const double dt = 1e-3;                 // Time step size in MD simulation
const int sampleFreq = 50000;           // Sampling frequency in MD simulation

// Cell variables
const int M = 600;                      // Number of cells
int L;                                  // Maximum total number of particles in the simulation
int* Npart = new int[M];                // Number of particles on each cell
const int N_avg = 40;                   // Average number of particles on each cell
const double B = 1;                     // Bond length between particles of cells
double* R = new double[M];              // Radii of cells
const double R_avg = N_avg*B/(2*pi);    // Mean radius of a cell
const double polyd = 0.08;              // Polydispersity in radius of cells (percentage)
double* A = new double[M];              // Initial area of cells
const double mass = 1;                  // Mass -- ( Viscous time scale is mass/gamma )
double* theta = new double[M];          // Angle between the particles in a cell
double* phi = new double[M];            // Polarity of the cells

// Initial grid
const double interDist_x = 3*R_avg;                                   // Interparticle distance
const double interDist_y = 3*R_avg;
const double initOffset_x = interDist_x;                                // Offset value for the initial positions
const double initOffset_y = interDist_y;
const int gridDim_x = (int)sqrt(M);                                     // Number of particles on the grid
const int gridDim_y = (int)sqrt(M);
const double gridSize_x = initOffset_x + interDist_x*gridDim_x;         // Size of the grid
const double gridSize_y = initOffset_y + interDist_y*gridDim_y;

// Boundary conditions
const double packing = 0.3;
const double cell_area = 77312.6;
const double boxSize_x = sqrt(cell_area/packing);
const double boxSize_y = sqrt(cell_area/packing);

//+++++++++++++++++++
//++ Force variables
//+++++++++++++++++++

// Spring force
const double kcons = 1e+4;          // Spring constant

// Bending force
const double Kbend = 100;           // Bending rigidity

// Area constraint force
const double A0_ = 0.9;             // 'Target'/Equilibrium area
double* A0 = new double[M];
const double Acons = 100;           // Membrane elasticity

// Drag force
const double gam = 10;              // Translational friction coefficient

// Diffusion
const double kB = 1;                // Boltzmann constant
const double T = 0.1;               // Temperature
const double Dt = kB*T/gam;         // Translational diffusion coefficient (kB*T/gam)

// Polarity rotation
const double Dr = 0.01;             // Rotational diffusion coefficient

// Motility force
double* Fm = new double[M];         // Motility force constant array
double* Fmsin = new double[M];
double* Fmcos = new double[M];
const double Fmc = 1;               // Motility force constant

// Lennard Jones interaction
const double eps = 10;                              // Depth (minimum) of the LJ potential
const double totalStep = 1e+8;                      // Total number of MD steps
const double sig = 2;                               // Distance where the LJ potential becomes zero
const double Rcut = 1.5*sig;                        // Cutoff distance for the LJ potential


//++++++++++++++++
//++ Verlet listing
//++++++++++++++++

// Verlet listing
const double Rver = 1.7*sig;                        // Verlet radius
const double skin = Rver-Rcut;                      // Skin for updating the Verlet list
int vlistsize;                                      // make a sane estimation for the vlist entry

// Boxing
const double Rbox_x = (double)ceil( boxSize_x/ifloor( boxSize_x/Rver ) );       // Box size for boxing
const double Rbox_y = (double)ceil( boxSize_y/ifloor( boxSize_y/Rver ) );
const int KX = ifloor( boxSize_x/Rbox_x );                                      // Number of boxes
const int KY = ifloor( boxSize_y/Rbox_y );
const int K2 = KX*KY;
int* hoc = new int[K2];                                                         // Head of chain


//+++++++++++++++++
//++ Other variables
//+++++++++++++++++

// Arrays
double* xc = new double[M];         // Initial centers of the cells
double* yc = new double[M];
double *x, *y, *vx, *vy, *Fx, *Fy, *xb, *yb, *xv, *yv;
int *link_list, *vlist_offset, *Nhead, *verList, *NV;

// Useful calculation variables
const double sampleDt = dt*sampleFreq;
const int totalSample = (int)totalStep/sampleFreq;
const double Rcut2 = Rcut*Rcut;
const double Rver2 = Rver*Rver;
const double skin2 = skin*skin;

