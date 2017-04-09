
//+++++++++++++++
//++ Data analysis with postprocessing for velocity correlators
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
    
    double th = 2*R_avg + Rcut;
    
    // +++++++++++
    
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++ Do the analysis
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    
    // Average speed per cell calculation
    double speed_per_cell_time_avg = 0.;
    
    
    // Mean velocities in each time frame
    vector<double> vx_mean(T-1,0);
    vector<double> vy_mean(T-1,0);
    
    
    // Spatial correlations
    double max_dist = 0.;
    int sug_max_dist = (int)fmax( Lx, Ly );
    double denom_time_avg = 0.;
    vector<double> cvv(sug_max_dist);
    for(int i = 0; i < sug_max_dist; i++){
        cvv[i] = 0.;
    }
    
    
    int num_neigh = 12;                         // Number of neighbours to check to make sure they are touching

    
    // Static correlations
    double max_cluster_size_avg = 0.;
    
    
    // Average pressure per cell calculation
    double avg_pres_time = 0.;
    
    // Moments of area
    double avg_mean_area = 0.;
    double avg_var_area = 0.;
    double avg_skew_area = 0.;
    double avg_kurt_area = 0.;

    
    
    // AT EACH TIME INSTANT ...
    for(int ti = 0; ti < T; ti++){                  // Time loop
        
        
        // Build KD tree
        const int MM = 9*M;
        vector<Point<2> > cell_vec(MM);

        
        // Average speed per cell calculation
        vector<double> speed(M);
        double speed_per_cell_ens_avg = 0.;
        
        
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
            cells[m].setImag(box, cell_vec);
            
            
            // Calculate velocities of cells and average speed of a cell
            if( ti > 0 ){
                cells[m].xvel[ti-1] = (cells[m].xpos[ti] - cells[m].xpos[ti-1])/dtSamp;
                cells[m].yvel[ti-1] = (cells[m].ypos[ti] - cells[m].ypos[ti-1])/dtSamp;
                double tmp_value = sqrt( cells[m].xvel[ti-1]*cells[m].xvel[ti-1] + cells[m].yvel[ti-1]*cells[m].yvel[ti-1] );
                speed[m] = tmp_value;
                speed_per_cell_ens_avg += tmp_value;
                
                // Mean velocities in a frame
                vx_mean[ti-1] += cells[m].xvel[ti-1];
                vy_mean[ti-1] += cells[m].yvel[ti-1];
            }
            
        }           // Cell loop
        
        
        // Build the KD tree
        KDtree<2> kdtree(cell_vec);
        
        
        // Neighbourhood check to be absolutely sure that cells are touching each other
        for(int m = 0; m < M; m++){
            
            vector<int> nn(num_neigh);              // Nearest neighbour vector
            vector<double> dn(num_neigh);           // Nearest neighbour distances
            
            // Create the neighbor list
            kdtree.nnearest(m,nn,dn,num_neigh);
            
            for (int i = 0; i < num_neigh; i++) {
                int out_flag = 0;
                int j = nn[i];                      // cell_vec index
                int jj = j % M;                     // cell index
                
                // Run over the particles of the original and the neighbor cell to make sure if they are indeed neighbors or not
                for(int n1 = 0; n1 < Npart[m]; n1++){       // n1
                    
                    // Break the for loop as the cells are already found to be neighbors
                    if (out_flag == 1){
                        break;
                    }
                    
                    for(int n2 = 0; n2 < Npart[jj]; n2++){   // n2
                        
                        double dx = particles[cells[jj].p_idx[n2]].xpos[ti] - particles[cells[m].p_idx[n1]].xpos[ti];
                        dx -= Lx*dnearbyint( dx/Lx );
                        double dy = particles[cells[jj].p_idx[n2]].ypos[ti] - particles[cells[m].p_idx[n1]].ypos[ti];
                        dy -= Ly*dnearbyint( dy/Ly );
                        
                        double dr2 = dx*dx + dy*dy;
                        
                        // If particles are within LJ cutoff distance, call them neighbors and break the for loop
                        if (dr2 < Rcut2) {
                            out_flag = 1;
                            cells[m].neighs.insert( jj );
                            break;
                        }
                        
                    }             // n2
                    
                }                 // n1
            }                     // it
        }                         // m
        
        
        //
        // By now, I have all the touching neighbors in the Cell object as set<int> neighs
        //

        
        
        // Mean velocities in a frame
        if( ti > 0 ){
            vx_mean[ti-1] /= M;
            vy_mean[ti-1] /= M;
        }
        
        
        // Average speed per cell calculation
        speed_per_cell_time_avg += speed_per_cell_ens_avg / M;
        
        
        // Pressure calculation
        nidx = -1;
        double mean_area = 0.;
        double avg_pres_ens = 0.;
        for (int m = 0; m < M; m++) {
            
            double area_of_cell = 0.;
            for (int n = 0; n < Npart[m]-1; n++) {
                ++nidx;
                area_of_cell += particles[nidx].xpos[ti]*particles[nidx+1].ypos[ti] - particles[nidx+1].xpos[ti]*particles[nidx].ypos[ti];
                
            }
            ++nidx;
            area_of_cell += particles[nidx].xpos[ti]*particles[nidx-Npart[m]+1].ypos[ti] - particles[nidx-Npart[m]+1].xpos[ti]*particles[nidx].ypos[ti];
            area_of_cell /= 2;
            cells[m].area = area_of_cell - A0[m];
            mean_area += cells[m].area;
            avg_pres_ens += 0.5*Acons*(area_of_cell - A0[m])*(area_of_cell - A0[m]);
            
        }
        
        avg_pres_time += avg_pres_ens/M;
        mean_area /= M;
        avg_mean_area += mean_area;
        
        // Moments of area
        
        // Second moment of area (standard deviation)
        double var = 0.;
        for(int m = 0; m < M; m++){
            double tmp_area = cells[m].area - mean_area;
            var += tmp_area*tmp_area;
        }
        var /= M;       // Not M-1, because we know the mean a priori
        double std_dev = sqrt(var);
        avg_var_area += var;
        
        // Third moment (skewness) and fourth moment (kurtosis)
        double skew = 0.;
        double kurt = 0.;
        for(int m = 0; m < M; m++){
            double tmp_area = (cells[m].area - mean_area)/std_dev;
            double tmp_area_3 = tmp_area*tmp_area*tmp_area;
            skew += tmp_area_3;
            kurt += tmp_area_3*tmp_area;
        }
        skew /= M;
        avg_skew_area += skew;
        kurt /= M;
        kurt -= 3;
        avg_kurt_area += kurt;
        
        
        
        
        
        // CORRELATION IN SPACE ...
        // spatial velocity correlation
        
        vector<double> cvv_avg(sug_max_dist);
        vector<double> cvv_cnt(sug_max_dist);
        for(int i = 0; i < sug_max_dist; i++){
            cvv_avg[i] = 0.;
            cvv_cnt[i] = 0.;
        }
        double denom_ens_avg = 0.;
        double denom_cnt = 0.;
        

        for(int m1 = 0; m1 < M; m1++){              // First cell loop
            for(int m2 = 0; m2 < M; m2++){          // Second cell loop
                
                double dx = cells[m2].xpos[ti] - cells[m1].xpos[ti];
                dx -= Lx*dnearbyint( dx/Lx );
                double dy = cells[m2].ypos[ti] - cells[m1].ypos[ti];
                dy -= Ly*dnearbyint( dy/Ly );
                double dr = sqrt( dx*dx + dy*dy );
                
                int bin = inearbyint( dr );
                
                // Spatial velocity correlation
                if( ti > 0 ){
                    double vx_1 = cells[m1].xvel[ti-1] - vx_mean[ti-1];
                    double vy_1 = cells[m1].yvel[ti-1] - vy_mean[ti-1];
                    double vx_2 = cells[m2].xvel[ti-1] - vx_mean[ti-1];
                    double vy_2 = cells[m2].yvel[ti-1] - vy_mean[ti-1];
                    
                    if( m1 == m2 ){
                        denom_ens_avg += vx_1*vx_2 + vy_1*vy_2;
                        denom_cnt++;
                    }
                    
                    cvv_avg[bin] += vx_1*vx_2 + vy_1*vy_2;
                    cvv_cnt[bin]++;
                }
                
                // Keep track of maximum distance between the cells
                if( dr > max_dist ){
                    max_dist = dr;
                }
            
        
            }      // First cell loop
        }          // Second cell loop
            
   
        // Ensemble averages for correlation in space
        for(int i = 0; i < sug_max_dist; i++){
            if( ti > 0 && cvv_cnt[i] != 0 ){
                cvv[i] += cvv_avg[i]/cvv_cnt[i];
            }
        }
        
        if( ti > 0 && denom_cnt != 0 ){
            denom_time_avg += denom_ens_avg/denom_cnt;
        }

        
        
        // Maximum cluster size
        // Set up the graph
        Graph graph( (int)M/5, th );
        int up_limit = graph.numNodes;
        double threshold2 = graph.threshold*graph.threshold;
        
        // Choose the fastest cells
        int break_cnt = 0;
        set<int> idx;
        vector<int> a = sort_indexes(speed);
        for(vector<int>::iterator it = a.begin(); it != a.end(); it++){
            idx.insert( *it );
            break_cnt++;
            if( break_cnt == up_limit ){
                break;
            }
        }
        
        // CHANGE THIS PART MAYBE??
        // Build the graph
        for (set<int>::iterator it = idx.begin(); it != idx.end(); it++) {  // loop over i
            
            set<int> destination = cells[*it].neighs;
            graph.addEdge( *it, destination );
            
        }  // loop over i
        
        // Search the fastest cells to find the connected components
        vector<int> max_cluster_size;
        vector<bool> visited(M);
        for (set<int>::iterator it = idx.begin(); it != idx.end(); it++) { visited[*it] = false; }
        stack<int> nodes;
        
        // For all the fastest nodes ...
        for (set<int>::iterator it = idx.begin(); it != idx.end(); it++) {
            
            // If they are not visited ...
            if ( !visited[*it] ) {
                max_cluster_size.push_back( (dfsFromNode(graph, *it, nodes, visited))/2 );
            }
        }
        
        max_cluster_size_avg += *max_element( max_cluster_size.begin(), max_cluster_size.end() );
        
        
        // Clear the neighbor lists, don't do neighbor analysis after this point
        for (int m = 0; m < M; m++) {
            cells[m].neighs.clear();
        }
        
        
    }           // Time loop
    
    
    // Average speed per cell calculation
    ofstream speed_per_cell_out;
    speed_per_cell_out.open( "avg_speed_per_cell.txt", ios::out | ios::app );
    speed_per_cell_out << packing << "\t" << Fmc << "\t" << eps << "\t" << speed_per_cell_time_avg / (T-1) << endl;
    speed_per_cell_out.close();
    
    
    // Average pressure per cell calculation
    ofstream pres_per_cell_out;
    pres_per_cell_out.open( "avg_pressure_per_cell.txt", ios::out | ios::app );
    pres_per_cell_out << packing << "\t" << Fmc << "\t" << eps << "\t" << avg_pres_time/T << endl;
    pres_per_cell_out.close();
    
    
    // Moments of area calculations
    
    ofstream mean_out;
    mean_out.open( "avg_mean_area_per_cell.txt", ios::out | ios::app );
    mean_out << packing << "\t" << Fmc << "\t" << eps << "\t" << avg_mean_area/T << endl;
    mean_out.close();
    
    ofstream var_out;
    var_out.open( "avg_var_area_per_cell.txt", ios::out | ios::app );
    var_out << packing << "\t" << Fmc << "\t" << eps << "\t" << avg_var_area/T << endl;
    var_out.close();
    
    ofstream skew_out;
    skew_out.open( "avg_skew_area_per_cell.txt", ios::out | ios::app );
    skew_out << packing << "\t" << Fmc << "\t" << eps << "\t" << avg_skew_area/T << endl;
    skew_out.close();
    
    ofstream kurt_out;
    kurt_out.open( "avg_kurt_area_per_cell.txt", ios::out | ios::app );
    kurt_out << packing << "\t" << Fmc << "\t" << eps << "\t" << avg_kurt_area/T << endl;
    kurt_out.close();
    
    
    // Spatial correlations
    int total_bins = (int)max_dist;
    ofstream sp_vacf_out("sp_vacf.txt");
    denom_time_avg /= (T-1);
    for(int i = 0; i < total_bins; i++){
        sp_vacf_out << i/(2*R_avg) << "\t" << cvv[i]/(denom_time_avg*(T-1)) << endl;
    }
    
    
    // Write the maximum cluster size average
    ofstream max_cluster_size_out;
    max_cluster_size_out.open( "max_cluster_size.txt", ios::out | ios::app );
    max_cluster_size_out << packing << "\t" << Fmc << "\t" << eps << "\t" << max_cluster_size_avg/T << endl;
    max_cluster_size_out.close();


    
    // DYNAMIC CORRELATIONS ...
    // Velocity autocorrelation
   
    ofstream vacf_out("vacf_cm.txt");
    
    // Normalization
    double denom_ens_avg = 0.;
    for(int m = 0; m < M; m++){
        double denom_time_avg = 0.;
        for(int ti = 0; ti < TD; ti++){
            double vx_1 = cells[m].xvel[ti];
            double vy_1 = cells[m].yvel[ti];
            
            denom_time_avg += vx_1*vx_1 + vy_1*vy_1;
        }
        denom_ens_avg += denom_time_avg/TD;
        
    }
    denom_ens_avg /= M;
    
    for(int di = 0; di < D; di++){              // Delay loop
        
        double cv_ens_avg = 0.;
        
        for(int m = 0; m < M; m++){             // Cell loop
            
            double cv_time_avg = 0.;
            
            for(int ti = 0; ti < TD; ti++){        // Different time origins
                
                // Velocity autocorrelation
                double vx_1 = cells[m].xvel[ti];
                double vx_2 = cells[m].xvel[ti+di];
                double vy_1 = cells[m].yvel[ti];
                double vy_2 = cells[m].yvel[ti+di];
                
                cv_time_avg += vx_1*vx_2 + vy_1*vy_2;
                
            }       // Different time origins
            
            cv_ens_avg += cv_time_avg/TD;
            
        }           // Cell loop
        
        // Velocity correlation
        double time_stamp = dtSamp*di;
        vacf_out << time_stamp << "\t" << cv_ens_avg/(M*denom_ens_avg) << endl;
        
    }               // Delay loop

    
} // end of the program
