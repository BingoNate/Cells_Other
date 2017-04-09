
//+++++++++++++++++
//++ Data analysis with postprocessing for density correlators
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
    
    
    // Spatial correlations
    double max_dist = 0.;
    int sug_max_dist = (int)fmax(Lx, Ly);
    vector<double> cnn(sug_max_dist, 0);

    // Static correlations
    int maxk = (int)(fmax(Lx,Ly)/2);
    vector<double> S_x_time_avg(maxk, 0);
    vector<double> S_y_time_avg(maxk, 0);
    
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
        
        
        // CORRELATION IN SPACE ...
        // Spatial bond order correlation,
        // pair correlation, and
        // spatial velocity correlation
        
        vector<double> cnn_avg(sug_max_dist);
        vector<double> c_cnt(sug_max_dist);
        for(int i = 0; i < sug_max_dist; i++){
            cnn_avg[i] = 0.;
            c_cnt[i] = 0.;
        }

        
        for(int ki = 0; ki < maxk; ki++){               // Wavevector loop
            
            double k = 2*pi*ki/maxk;
            
            double S_cos_x = 0.;
            double S_sin_x = 0.;
            double S_cos_y = 0.;
            double S_sin_y = 0.;
            double S_cnt = 0.;
            
            for(int m1 = 0; m1 < M; m1++){              // First cell loop
                for(int m2 = 0; m2 < M; m2++){          // Second cell loop
                    
                    double dx = cells[m2].xpos[ti] - cells[m1].xpos[ti];
                    dx -= Lx*dnearbyint( dx/Lx );
                    double dy = cells[m2].ypos[ti] - cells[m1].ypos[ti];
                    dy -= Ly*dnearbyint( dy/Ly );
                    double dr = sqrt( dx*dx + dy*dy );
                    
                    if(ki == 0){
                        int bin = inearbyint( dr );
                        c_cnt[bin]++;
                        cnn_avg[bin] += 1;
                        
                        // Keep track of maximum distance between the cells
                        if( dr > max_dist ){
                            max_dist = dr;
                        }
                    }
                    else{
                        
                        S_cnt++;
                        double tut = k*dx;
                        S_cos_x += cos(tut);
                        S_sin_x += sin(tut);
                        
                        tut = k*dy;
                        S_cos_y += cos(tut);
                        S_sin_y += sin(tut);
                        
                        
                    }
                    
                }      // First cell loop
            }          // Second cell loop
            
            S_x_time_avg[ki] += (S_cos_x*S_cos_x + S_sin_x*S_sin_x)/S_cnt;
            S_y_time_avg[ki] += (S_cos_y*S_cos_y + S_sin_y*S_sin_y)/S_cnt;
            
        }               // Wavevector loop
        
        
        // Ensemble averages for correlation in space
        for(int i = 0; i < sug_max_dist; i++){
            cnn[i] += cnn_avg[i];
        }
        
        
    }           // Time loop
    
    
    
    // Spatial correlations
    int total_bins = (int)max_dist;
    ofstream pair_out("pair_corr.txt");
    double rho = M/(Lx*Ly);
    for(int i = 1; i < total_bins; i++){
            double n_ideal = 2*pi*i*rho;
            pair_out << i/(2*R_avg) << "\t" << cnn[i]/(n_ideal*T*M) << endl;
        
    }
    
    
    // Static correlations
    ofstream static_struct_out("static_struct.txt");
    for(int ki = 1; ki < maxk; ki++){
        double k = 2*pi*ki/maxk;
        S_x_time_avg[ki] /= T;
        S_y_time_avg[ki] /= T;
        static_struct_out << k*2*R_avg << "\t" << (S_x_time_avg[ki]+S_y_time_avg[ki])/(2*M) << endl;
    }
    
    
    
    
    // DYNAMIC CORRELATIONS ...
    // Intermediate scattering function
    // Dynamic structure factor
    
    ofstream inter_scattering_out("inter_scattering.txt");
    double k_int_scatter = 0.5;
    vector<double> Fr(D);
    vector<double> Fi(D);
    int maxw = (int)D/2;
    
    for(int di = 0; di < D; di++){              // Delay loop
        
        double cv_ens_avg = 0.;
        
        double Fr_x = 0.;
        double Fr_y = 0.;
        double Fi_x = 0.;
        double Fi_y = 0.;
        
        double msd_ens_avg = 0.;
        
        for(int m = 0; m < M; m++){             // Cell loop
            
            double F_cos_x = 0.;
            double F_sin_x = 0.;
            double F_cos_y = 0.;
            double F_sin_y = 0.;
            
            for(int ti = 0; ti < TD; ti++){        // Different time origins
                
                
                // Density autocorrelation in Fourier space (intermediate scattering function)
                double dx = cells[m].xpos[ti+di] - cells[m].xpos[ti];
                double dy = cells[m].ypos[ti+di] - cells[m].ypos[ti];
                
                double tut = k_int_scatter*dx;
                F_cos_x += cos( tut );
                F_sin_x += sin( tut );
                
                tut = k_int_scatter*dy;
                F_cos_y += cos( tut );
                F_sin_y += sin( tut );
                
            }       // Different time origins
            
            Fr_x += F_cos_x/TD;
            Fr_y += F_cos_y/TD;
            Fi_x += F_sin_x/TD;
            Fi_y += F_sin_y/TD;
            
        }           // Cell loop
        
        double time_stamp = dtSamp*di;
        
        // Intermediate scattering
        Fr[di] = (Fr_x+Fr_y)/(2*M);
        Fi[di] = (Fi_x+Fi_y)/(2*M);
        
        if( di != 0 ){
            inter_scattering_out << time_stamp << "\t" << Fr[di] << endl;
        }
        else{
            inter_scattering_out << 0 << "\t" << 1 << endl;
        }
        
    }               // Delay loop
    
    
    // Dynamic structure factor
    ofstream dync_struct_out("dynamic_struct.txt");
    for(int di = 0; di < D; di++){
        
        double sumreal = 0.;
        double sumimag = 0.;
        double jj = di-maxw;
        double phase = 2*pi*jj/D;
        
        for(int ii = 0; ii < D; ii++){      // Frequency loop
            
            sumreal += Fr[ii]*cos(phase*ii) + Fi[ii]*sin(phase*ii);
            sumimag += Fi[ii]*cos(phase*ii) - Fr[ii]*sin(phase*ii);
            
        }   // Frequency loop
        
        sumreal /= D;
        sumimag /= D;
        
        dync_struct_out << phase << "\t" << (sumreal*sumreal + sumimag*sumimag)/(2*pi) << endl;
    }
    
} // end of the program
