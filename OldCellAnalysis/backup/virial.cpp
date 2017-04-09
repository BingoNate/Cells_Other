
//+++++++++++++++++
//++ Data analysis with postprocessing for virial stress calculation
//+++++++++++++++++

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <map>
#include <algorithm>
#include "variables.h"
#include "basic.h"
#include "pointbox.h"
#include "kdtree.h"
#include "data_struct.cpp"

using namespace std;

#define pi M_PI


double sig2 = sig*sig;
double eps48 = 48*eps;

// Calculate virial stress
double calcVirStress( double dx, double dy ){
    
    double d2 = dx*dx + dy*dy;
    double Fx = 0.;
    double Fy = 0.;
    if( d2 < Rcut2 ){
        
        double p3 = sig2/d2;
        p3 = p3*p3*p3;
        double tempFx = eps48*( p3*p3 - p3/2. )*dx/d2;
        double tempFy = eps48*( p3*p3 - p3/2. )*dy/d2;
        Fx = tempFx;
        Fy = tempFy;
    }
    
    return dx*Fx + dy*Fy;

}

// ++++++



int main( int argc, char *argv[] ){
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++ Get the data
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    // Choose the data to be used
    const int nsamp = sampleFreq;           // Sampling interval
    const double dtSamp = nsamp*dt;         // Time units between two data points
    
    // Set general simulation variables
    Data simData;
    simData.setNumStep( totalStep/sampleFreq );
    simData.setNumDelay( sqrt( totalStep/sampleFreq ) );
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

    // Instantiate the data structure
    vector<Particle> particles;
    for(int i = 0; i < L; i++){
        Particle p(L, T);
        particles.push_back( p );
    }
    vector<Cell> cells;
    for(int m = 0; m < M; m++){
        Cell c(M, N);
        cells.push_back( c );
    }
    
    // Set filenames
    char * input_file = argv[1];        // Filename for the particle position data in xyz file format
    
    // Read the data to the cells
    int sample_cnt = -1;
    int sample_data = 1;
    readPosData( particles, input_file, sample_cnt, sample_data );
    
    // Box link lists and Verlet list
    CellBox cellbox(Rver, L);
    cellbox.setBoxProp( box );
    vector<int> numVerList(L, 0);
    int verListSize = L*1024;
    vector<int> verList(verListSize, -1);
    vector<int> verListOffset(L, 0);
    verletList( particles, cellbox, Lx, Ly, numVerList, verList, verListOffset, verListSize, box, 0 );

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++ Do the analysis
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    
    // GLOBAL VIRIAL STRESS CALCULATION
    
    double vir_avg = 0.;
    
    for(int ti = 0; ti < T; ti++){  // Time loop
        
        updateList(particles, cellbox, Lx, Ly, numVerList, verList, verListOffset, verListSize, box, ti);
        
        double vir_ens = 0.;
        double vir_cnt = 0.;
        for(int i = 0; i < L; i++){ // Particle loop
            
            int jmax = verListOffset[i] + numVerList[i];
            for(int jj = verListOffset[i]; jj < jmax; jj++){
                
                double dx = particles[i].xpos[ti] - particles[verList[jj]].xpos[ti];
                dx -= Lx*dnearbyint( dx/Lx );
                double dy = particles[i].ypos[ti] - particles[verList[jj]].ypos[ti];
                dy -= Ly*dnearbyint( dy/Ly );
                
                vir_ens += calcVirStress(dx, dy);
                vir_cnt++;
            }
            
        }   // Particle loop
        
        // Take the ensemble average
        if( vir_cnt != 0 ){
            vir_avg += vir_ens/L;   // Should I normalise by the total number of particles or the counted particles that are actually interacting?
        }
        
    }   // Time loop
    
    // Take the time average
    vir_avg /= T;

    // Save the data
    string glob_vir_path = "glob_vir_stress.txt";
    ofstream glob_vir_out;
    glob_vir_out.open( glob_vir_path, ios::out | ios::app );
    glob_vir_out << dens << "\t" << Fmc << "\t" << eps << "\t" << vir_avg << endl;
    glob_vir_out.close();

}
