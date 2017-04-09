
//+++++++++++++++
//++ Data analysis with postprocessing for tensile or compressive detection
//+++++++++++++++

#include "data_struct.cpp"


int main( int argc, char *argv[] ){
    
    // +++++++++++
    // Get the simulation data
    
    const int nsamp = sampleFreq;                       // Sampling interval
    const double dtSamp = nsamp*dt;                     // Time stDr between two data points
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
    
    
    int num_neigh = 12;                         // Number of neighbours to check to make sure they are touching
    int boundary_identifier_value = 5;          // Identify the bulk and the boundaries of the cluster
    
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

        
        // Pressure calculation
        nidx = -1;
        double mean_area = 0.;
        double avg_pres_ens = 0.;
        double cnt_m = 0.;
        for (int m = 0; m < M; m++) {
            
            double area_of_cell = 0.;
            for (int n = 0; n < Npart[m]-1; n++) {
                ++nidx;
                area_of_cell += particles[nidx].xpos[ti]*particles[nidx+1].ypos[ti] - particles[nidx+1].xpos[ti]*particles[nidx].ypos[ti];
                
            }
            ++nidx;
            area_of_cell += particles[nidx].xpos[ti]*particles[nidx-Npart[m]+1].ypos[ti] - particles[nidx-Npart[m]+1].xpos[ti]*particles[nidx].ypos[ti];
            area_of_cell /= 2;
            
            set<int> neighs = cells[m].neighs;
            int size_value = neighs.size();
            if (size_value < boundary_identifier_value) {
                cells[m].area = area_of_cell - A0[m];
                mean_area += cells[m].area;
                avg_pres_ens += 0.5*Acons*(area_of_cell - A0[m])*(area_of_cell - A0[m]);
                cnt_m++;
            }
            
        }
        if (cnt_m != 0) {
            avg_pres_time += avg_pres_ens/cnt_m;
            mean_area /= cnt_m;
            avg_mean_area += mean_area;
        }
        
        // Moments of area
        
        // Second moment of area (standard deviation)
        double var = 0.;
        for(int m = 0; m < M; m++){
            double tmp_area = cells[m].area - mean_area;
            set<int> neighs = cells[m].neighs;
            int size_value = neighs.size();
            if (size_value < boundary_identifier_value) {
                var += tmp_area*tmp_area;
            }
        }
        
        double std_dev = 1e-6;
        if (cnt_m != 0) {
            var /= cnt_m;       // Not M-1, because we know the mean a priori
            std_dev = sqrt(var);
            avg_var_area += var;
        }
        
        // Third moment (skewness) and fourth moment (kurtosis)
        double skew = 0.;
        double kurt = 0.;
        for(int m = 0; m < M; m++){
            double tmp_area = (cells[m].area - mean_area)/std_dev;
            double tmp_area_3 = tmp_area*tmp_area*tmp_area;
            set<int> neighs = cells[m].neighs;
            int size_value = neighs.size();
            if (size_value < boundary_identifier_value) {
                skew += tmp_area_3;
                kurt += tmp_area_3*tmp_area;
            }
        }
        if (cnt_m != 0) {
            skew /= cnt_m;
            avg_skew_area += skew;
            kurt /= cnt_m;
            kurt -= 3;
            avg_kurt_area += kurt;
        }
        
        
    }           // Time loop
    
    
    // Average pressure per cell calculation
    ofstream pres_per_cell_out;
    pres_per_cell_out.open( "avg_pressure_per_cell.txt", ios::out | ios::app );
    pres_per_cell_out << Kbend << "\t" << Acons << "\t" << Dr << "\t" << avg_pres_time/T << endl;
    pres_per_cell_out.close();
    
    
    // Moments of area calculations
    
    ofstream mean_out;
    mean_out.open( "avg_mean_area_per_cell.txt", ios::out | ios::app );
    mean_out << Kbend << "\t" << Acons << "\t" << Dr << "\t" << avg_mean_area/T << endl;
    mean_out.close();
    
    ofstream var_out;
    var_out.open( "avg_var_area_per_cell.txt", ios::out | ios::app );
    var_out << Kbend << "\t" << Acons << "\t" << Dr << "\t" << avg_var_area/T << endl;
    var_out.close();
    
    ofstream skew_out;
    skew_out.open( "avg_skew_area_per_cell.txt", ios::out | ios::app );
    skew_out << Kbend << "\t" << Acons << "\t" << Dr << "\t" << avg_skew_area/T << endl;
    skew_out.close();
    
    ofstream kurt_out;
    kurt_out.open( "avg_kurt_area_per_cell.txt", ios::out | ios::app );
    kurt_out << Kbend << "\t" << Acons << "\t" << Dr << "\t" << avg_kurt_area/T << endl;
    kurt_out.close();
  
    
} // end of the program


