
//+++++++++++++++
//++ Data analysis with postprocessing for mean square displacement
//+++++++++++++++


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
    
    double cell_diam_avg = 4*R_avg*R_avg;
    
    // +++++++++++
        
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++ Do the analysis
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    // Mean positions in each time frame
    vector<double> x_mean(T,0);
    vector<double> y_mean(T,0);
    
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
            
            
            // Mean positions in a frame
            x_mean[ti] += cells[m].xpos[ti];
            y_mean[ti] += cells[m].ypos[ti];
            
        }           // Cell loop

        
        // Mean positions in a frame
        x_mean[ti] /= M;
        y_mean[ti] /= M;
    
        
    }               // Time loop
        

    // DYNAMIC CORRELATIONS
    ofstream msd_out("msd_post.txt");
    
    for(int di = 0; di < D; di++){              // Delay loop
        
        double msd_ens_avg = 0.;
        
        for(int m = 0; m < M; m++){             // Cell loop
            
            double msd_time_avg = 0.;
            
            for(int ti = 0; ti < TD; ti++){        // Different time origins
                
                // Mean square displacement
                double x_1 = cells[m].xpos[ti] - x_mean[ti];
                double x_2 = cells[m].xpos[ti+di] - x_mean[ti];
                double y_1 = cells[m].ypos[ti] - y_mean[ti];
                double y_2 = cells[m].ypos[ti+di] - y_mean[ti];
                
                double delta_x = x_2 - x_1;
                double delta_y = y_2 - y_1;
                
                msd_time_avg += (delta_x*delta_x + delta_y*delta_y)/cell_diam_avg;
                
            }       // Different time origins
            
            msd_ens_avg += msd_time_avg/TD;
            
        }           // Cell loop
        
        double time_stamp = dtSamp*di;
    
        // Mean square displacement
        msd_out << time_stamp << "\t" << msd_ens_avg/M << endl;
        
    }               // Delay loop

    
    // ++++++++++++++
    
}
