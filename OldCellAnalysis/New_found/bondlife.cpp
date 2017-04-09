
//+++++++++++++++
//++ Data analysis with postprocessing for bondlife correlators
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
    double T = simData.getNumStep();                    // Total time
    double D = simData.getNumDelay();                   // Last delay time
    double TD = simData.getNumTotalDelay();             // Total time - last delay time
    
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
    int boundary_identifier_value = 4;          // If the cell has less than or equal to neighbors w.r.t. this value, it is marked as a boundary cell
    
    // Hash function variables for the time for the bond
    const int MM = 9*M;
    const int TD2 = T/2;
    const int D2 = sqrt(TD2);
    int keysize = M*M*TD2;
    Hash<Ullong, int, Nullhash> *bondtime = new Hash<Ullong, int, Nullhash>(keysize, keysize);
    Ranhash hashfunc;
    
    
    int t_offset = 0;
    int t_total = T-t_offset;
    for (int ti = t_offset; ti < T; ti++) {            // Time loop
        
        // Build KD tree
        vector<Point<2> > cell_vec(MM);
        
        
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
            
            vector<int> nn(num_neigh);              // Nearest neighbour indices vector
            vector<double> dn(num_neigh);           // Nearest neighbour distances vector
            
            // Create the neighbor list
            kdtree.nnearest(m,nn,dn,num_neigh);
            
            for (int i = 0; i < num_neigh; i++) {
                int out_flag = 0;
                int j = nn[i];                      // cell_vec index
                int jj = j % M;                     // cell index
                
                // Run over the particles of the original and the neighbor cell to make sure if they are indeed neighbors or not
                for(int n1 = 0; n1 < Npart[m]; n1++){       // n1
                    
                    // Break the for loop as the cells are already found to be neighbors
                    if (out_flag == 1) {
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
        
        if (ti < TD2) {
        
        
            // Store bonds between bulk cells at each time step
            for(int m = 0; m < M; m++){       // Cell loop
                
                // Fill out the hash table with zeroes and later identify the neighbors as ones
                for (int m2 = 0; m2 < M; m2++) {
                    Ullong key = hashfunc.int64(m) ^ hashfunc.int64(m2) ^ hashfunc.int64(ti);
                    bondtime -> set(key, 0);
                    
                    if (m == m2) {
                        bondtime -> erase(key);
                        bondtime -> set(key, 1);
                    }
                    
                }
                
                // Find the neighbors of the bulk cell currently being iterated
                set<int> neighs = cells[m].neighs;
                int size_value = neighs.size();
                if (size_value > boundary_identifier_value) {
                    
                    
                    // Add the neighbors to the bondtime hash table
                    for (set<int>::iterator it = neighs.begin(); it != neighs.end(); it++) {
                        Ullong key = hashfunc.int64(m) ^ hashfunc.int64(*it) ^ hashfunc.int64(ti);
                        bondtime -> erase(key);
                        bondtime -> set(key, 1);
                    }
                    
                }
                
            }           // Cell loop
            
        }
        
        
        // DELETE THE NEIGHBOR LIST
        for(int m = 0; m < M; m++){
            cells[m].neighs.clear();
        }
        
        
        //
        // NO FURTHER ANALYSIS ON NEIGHBORS BEYOND THIS POINT
        //
        

        
        
    }       // End of time loop
    
    
    // Calculate autocorrelation of bond life time
    ofstream cbl_out("bondlife_corr.txt");
    for (int di = 0; di < D2; di++) {                    // Delay loop
        
        double cbl_time_avg = 0.;
        
        for (int ti = t_offset; ti < TD2; ti++) {        // Time origin loop
            
            // Calculate the neighbors
            
            // Build KD tree
            vector<Point<2> > cell_vec(MM);
            
            
            // Calculate images
            for(int m = 0; m < M; m++){                 // Cell loop
                
                // Set image positions
                cells[m].getImag(box, ti);
                cells[m].setImag(box, cell_vec);
                
                
            }           // Cell loop
            
            
            // Build the KD tree
            KDtree<2> kdtree(cell_vec);
            
            
            // Neighbourhood check to be absolutely sure that cells are touching each other
            for(int m = 0; m < M; m++){
                
                vector<int> nn(num_neigh);              // Nearest neighbour indices vector
                vector<double> dn(num_neigh);           // Nearest neighbour distances vector
                
                // Create the neighbor list
                kdtree.nnearest(m,nn,dn,num_neigh);
                
                for (int i = 0; i < num_neigh; i++) {
                    int out_flag = 0;
                    int j = nn[i];                      // cell_vec index
                    int jj = j % M;                     // cell index
                    
                    // Run over the particles of the original and the neighbor cell to make sure if they are indeed neighbors or not
                    for(int n1 = 0; n1 < Npart[m]; n1++){       // n1
                        
                        // Break the for loop as the cells are already found to be neighbors
                        if (out_flag == 1) {
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

            
            // Calculate the correlation function
            
            double cbl_ens_avg = 0.;
            double cnt_ens = 0.;
            
            for (int m = 0; m < M; m++) {               // Run over cells
                
                // Check if the cell is in the bulk of the cluster
                set<int> neighs = cells[m].neighs;
                int size_value = neighs.size();
                if (size_value > boundary_identifier_value) {
                    
                    for (set<int>::iterator it = neighs.begin(); it != neighs.end(); it++) {        // Neighbor loop
                        
                        int j1, j2;
                        Ullong key1 = hashfunc.int64(m) ^ hashfunc.int64(*it) ^ hashfunc.int64(ti);
                        bondtime -> get(key1, j1);
                        Ullong key2 = hashfunc.int64(m) ^ hashfunc.int64(*it) ^ hashfunc.int64(ti+di);
                        bondtime -> get(key2, j2);
                        
                        cbl_ens_avg += j1*j2;
                        cnt_ens++;
                        
                    }   // Neighbor loop
                
                }       // if cell is in the bulk of the cluster
                
            }           // Run over cells
            
            if (cnt_ens != 0) {
                cbl_time_avg += cbl_ens_avg/cnt_ens;
            }
            
            // DELETE THE NEIGHBOR LIST
            for(int m = 0; m < M; m++){
                cells[m].neighs.clear();
            }
            
            
            //
            // NO FURTHER ANALYSIS ON NEIGHBORS BEYOND THIS POINT
            //
 
            
        }       // Time orisgin loop
        
        cbl_time_avg /= (TD2-t_offset);
        cbl_out << di*dtSamp << "\t" << cbl_time_avg << endl;
        
    }           // Delay loop






}


