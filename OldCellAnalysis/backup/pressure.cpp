
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
    
    // Radii information
    ifstream r_in;
    r_in.open( "radii.dat", ios::binary | ios::in );
    r_in.read( (char *)R, M*sizeof(double) );
    r_in.close();
    for(int m = 0; m < M; m++){
        B[m] = pi*R[m]/20;
        A[m] = pi*R[m]*R[m];
        A0[m] = A0_*A[m];
    }
    
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

    
    // BOXED VIRIAL STRESS CALCULATION
    
    // Create the box structure
    double box_grid = 5.;
    double rsep = Lx/box_grid;
    CellBox virbox(rsep, L);
    virbox.setBoxProp( box );
    
    // Print box info
    ofstream box_info_out("img/dat/box_info.txt");
    box_info_out << virbox.Bx << "\t" << virbox.By << "\t" << virbox.Rx << "\t" << virbox.Ry << endl;

//    string boxvir_out_path = "img/dat/local_vir_stress_0.txt";
//    ofstream boxvir_out( boxvir_out_path );
//    boxvir_out << 0 << "\t" << 0 << "\t" << 0 << endl;
//    string pres_out_path = "img/dat/local_pressure_0.txt";
//    ofstream pres_out( pres_out_path );
//    pres_out << 0 << "\t" << 0 << endl;
    
    vector<double> vir_avg(virbox.BN, 0.);;
//    vector<double> P(M, 0.);

    for(int ti = 1; ti < T; ti++){      // Time loop
        
        initList( particles, virbox, box, ti );

        for(int i = 0; i < L; i++){     // Particle loop
            
            particles[i].getImag(box, ti);
            int kx = boxNum_periodic( particles[i].xi, virbox.Rx, virbox.Bx );
            int ky = boxNum_periodic( particles[i].yi, virbox.Ry, virbox.By );
            
            int kbox = kx*virbox.By + ky;
            int jp = virbox.head_box[kbox];
            
            vector<double> vir_ens(virbox.BN, 0.);
            vector<double> vir_cnt(virbox.BN, 0.);
            
            while( jp != -1 ){          // Link list
                
                if( i != jp && i/N != jp/N ){
                    
                    double dx = particles[i].xpos[ti] - particles[jp].xpos[ti];
                    double dy = particles[i].ypos[ti] - particles[jp].ypos[ti];
                    double d2 = dx*dx + dy*dy;
                    
                    if( d2 < Rcut2 ){
                        vir_ens[kbox] += calcVirStress(dx, dy);
                        vir_cnt[kbox]++;
                    }
                }
                
                jp = virbox.link_box[jp];
                
            }   // Link list
            
            if( vir_cnt[kbox] != 0 ){
                vir_avg[kbox] += vir_ens[kbox]/vir_cnt[kbox];
            }
            
        }       // Particle loop
        
        string boxvir_out_path = "img/dat/local_vir_stress_" + to_string(ti) + ".txt";
        ofstream boxvir_out( boxvir_out_path );
        
        for(int j = 0; j < virbox.BN; j++){
            int jx = j/virbox.By;
            int jy = j % virbox.By;
            boxvir_out << jx << "\t" << jy << "\t" << vir_avg[j]/ti << endl;
        }
        
        // PRESSURE CALCULATION
       
	vector<double> P(M, 0.); 
        string pres_out_path = "img/dat/local_pressure_" + to_string(ti) + ".txt";
        ofstream pres_out( pres_out_path );
        
        for(int m = 0; m < M; m++){
            
            double area_of_cell = 0.;
            int i;
            for(int n = 0; n < N-1; n++){
                i = m*N + n;
                area_of_cell += particles[i].xpos[ti]*particles[i+1].ypos[ti] - particles[i+1].xpos[ti]*particles[i].ypos[ti];
                
            }
            i = m*N + (N-1);
            area_of_cell += particles[i].xpos[ti]*particles[0].ypos[ti] - particles[0].xpos[ti]*particles[i].ypos[ti];
            area_of_cell /= 2; 
            P[m] += 0.5*Acons*(area_of_cell - A0[m])*(area_of_cell - A0[m]);
            
            pres_out << m << "\t" << P[m] << endl;
        }
        
    }           // Time loop

    
    // ++++++


}
