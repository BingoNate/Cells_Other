
//+++++++++++++++++
//++ Data analysis with postprocessing for virial stress calculation
//+++++++++++++++++

#include "data_struct.cpp"


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
    
    
    double vir_time_avg = 0.;
    double sig2 = sig*sig;
    double eps48 = 48*eps;
    
    
    // AT EACH TIME INSTANT ...
    for(int ti = 0; ti < T; ti++){                  // Time loop
        
        
        // Calculate positions of cells
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
        
        
        // Interparticle virial stress
        updateList(particles, cellbox, Lx, Ly, numVerList, verList, verListOffset, verListSize, box, ti);
        
        double vir_ens_avg = 0.;
        
        for(int i = 0; i < L; i++){     // Particle loop
            
            double vir_neigh_avg = 0.;
            double vir_cnt = 0.;
            
            int jmax = verListOffset[i] + numVerList[i];
            for(int jj = verListOffset[i]; jj < jmax; jj++){
                
                double dx = particles[verList[jj]].xpos[ti] - particles[i].xpos[ti];
                dx -= Lx*dnearbyint( dx/Lx );
                double dy = particles[verList[jj]].ypos[ti] - particles[i].ypos[ti];
                dy -= Ly*dnearbyint( dy/Ly );
                
                vir_neigh_avg += calcVirStress(dx, dy, eps48, sig2);
                vir_cnt++;
                
            }
            
            if( vir_cnt != 0 ){
                vir_ens_avg += vir_neigh_avg/vir_cnt;
            }
            
        }       // Particle loop
        
        vir_time_avg += vir_ens_avg / L;

        
    }       // Time loop
    
    ofstream glob_vir_out;
    glob_vir_out.open( "glob_vir_stress.txt", ios::out | ios::app );
    glob_vir_out << packing << "\t" << Fmc << "\t" << eps << "\t" << vir_time_avg/T << endl;
    glob_vir_out.close();
    
}
