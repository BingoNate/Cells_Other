
//+++++++++++++++++
//++ Data analysis with postprocessing for virial stress calculation
//+++++++++++++++++

#include "data_struct.cpp"

using namespace std;

#define pi M_PI


double sig2 = sig*sig;
double eps48 = 48*eps;


int main( int argc, char *argv[] ){
    
    // +++++++++++
    // Get the simulation data
    
    const int nsamp = sampleFreq;                       // Sampling interval
    const double dtSamp = nsamp*dt;                     // Time steps between two data points
    Data simData;
    simData.setNumStep( totalStep/sampleFreq );
    simData.setNumDelay( sqrt( totalStep/sampleFreq ) );
    simData.setNumTotalDelay();
    double T = simData.getNumStep();              // Total time
    double D = simData.getNumDelay();             // Last delay time
    double TD = simData.getNumTotalDelay();       // Total time - last delay time
    
    SimBox box;
    box.setWidth( boxSize_x );
    box.setHeight( boxSize_y );
    const double Lx = box.getWidth();
    const double Ly = box.getHeight();
    
    char * input_file = argv[1];                        // Filename for the particle position data in xyz file format
    
    readRadii();
    
    vector<Cell> cells;
    for(int m = 0; m < M; m++){
        Cell c(m, T, Npart[m], M);
        cells.push_back( c );
    }
    
    vector<Particle> particles;
    int nidx = -1;
    for(int m = 0; m < M; m++){
        for(int n = 0; n < Npart[m]; n++){
            ++nidx;
            Particle p(T, m);
            particles.push_back( p );
            cells[m].p_idx.push_back(nidx);
        }
    }
    
    int sample_cnt = -1;
    int sample_data = 1;
    int totalSamp = readPosData( particles, input_file, sample_cnt, sample_data );
    if( T != totalSamp ){
        simData.setNumStep( totalSamp );
        simData.setNumDelay( sqrt( totalSamp ) );
        simData.setNumTotalDelay();
        T = simData.getNumStep();              // Total time
        D = simData.getNumDelay();             // Last delay time
        TD = simData.getNumTotalDelay();       // Total time - last delay time
    }
    
    CellBox cellbox(Rver, L);
    cellbox.setBoxProp( box );
    int verListSize = L*1024;
    vector<int> numVerList(L, 0);
    vector<int> verList(verListSize, -1);
    vector<int> verListOffset(L, 0);
    verletList( particles, cellbox, Lx, Ly, numVerList, verList, verListOffset, verListSize, box, 0 );
    
    // +++++++++++
    
    
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
    
    vector<double> vir_avg(virbox.BN, 0.);
    
    double sig2 = sig*sig;
    double eps48 = 48*eps;
    
    
    // AT EACH TIME INSTANT ...
    for(int ti = 0; ti < T; ti++){                  // Time loop
        
        nidx = -1;
        for(int m = 0; m < M; m++){                 // Cell loop
            
            double xcom = 0.;
            double ycom = 0.;
            for(int n = 0; n < Npart[m]; n++){      // Particle loop
                ++nidx;
                xcom += particles[nidx].xpos[ti];
                ycom += particles[nidx].ypos[ti];
            }       // Particle loop
            cells[m].xpos[ti] = xcom/Npart[m];
            cells[m].ypos[ti] = ycom/Npart[m];
            
            // Set image positions
            cells[m].getImag(box, ti);
            
    
        }           // Cell loop
        
    }

    
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
                
                if( i != jp && particles[i].cell_idx != particles[jp].cell_idx ){
                    
                    double dx = particles[jp].xpos[ti] - particles[i].xpos[ti];
                    double dy = particles[jp].ypos[ti] - particles[i].ypos[ti];
                    double d2 = dx*dx + dy*dy;
                    
                    if( d2 < Rcut2 ){
                        vir_ens[kbox] += calcVirStress(dx, dy, eps48, sig2);
                        vir_cnt[kbox]++;
                    }
                }
                
                jp = virbox.link_box[jp];
                
            }   // Link list
            
            if( vir_cnt[kbox] != 0 ){
                vir_avg[kbox] += vir_ens[kbox]/vir_cnt[kbox];
            }
            
        }       // Particle loop
        
        
        char buffer[30];
        sprintf(buffer, "%d", ti);
        char boxvir_out_path[200] = "img/dat/local_vir_stress_";
        strcat(boxvir_out_path, buffer);
        char buffer2[30] = ".txt";
        strcat(boxvir_out_path, buffer2);
        ofstream boxvir_out( boxvir_out_path );
        
        for(int j = 0; j < virbox.BN; j++){
            int jx = j/virbox.By;
            int jy = j % virbox.By;
            boxvir_out << jx << "\t" << jy << "\t" << vir_avg[j]/ti << endl;
        }
        
        // PRESSURE CALCULATION
        vector<double> P(M, 0.);
        char pres_out_path[200] = "img/dat/local_pressure_";
        strcat(pres_out_path, buffer);
        strcat(pres_out_path, buffer2);
        ofstream pres_out( pres_out_path );
        
        nidx = -1;
        for(int m = 0; m < M; m++){
            
            double area_of_cell = 0.;
            for(int n = 0; n < Npart[m]-1; n++){
                ++nidx;
                area_of_cell += particles[nidx].xpos[ti]*particles[nidx+1].ypos[ti] - particles[nidx+1].xpos[ti]*particles[nidx].ypos[ti];
                
            }
            ++nidx;
            area_of_cell += particles[nidx].xpos[ti]*particles[nidx-Npart[m]+1].ypos[ti] - particles[nidx-Npart[m]+1].xpos[ti]*particles[nidx].ypos[ti];
            area_of_cell /= 2; 
            P[m] += 0.5*Acons*(area_of_cell - A0[m])*(area_of_cell - A0[m]);
            
            pres_out << m << "\t" << P[m] << endl;
        }
        
    }           // Time loop

    
    // ++++++


}
