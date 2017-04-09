

//+++++++++++++++
//++ Data analysis with postprocessing for velocity correlators
//+++++++++++++++

#include "data_struct_debug.cpp"


int main( int argc, char *argv[] ){
    
    // +++++++++++
    // Get the simulation data
    
    const int nsamp = sampleFreq;                       // Sampling interval
    const double dtSamp = nsamp*dt;                     // Time steps between two data points
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
    
    
    // Bond order parameter
    int psin = 4;
    double psin_time_avg = 0.;
    double susc_sq_avg = 0.;
    
    
    int num_neigh = 12;                         // Number of neighbours to check to make sure they are touching
    int boundary_identifier_value = 4;          // If the cell has less than or equal to neighbors w.r.t. this value, it is marked as a boundary cell
    
    
    // Average neighbor change
    double neigh_chng_cnt_avg = 0.;
    double neigh_chng_cnt_avg_in_bulk = 0.;
//    vector<set<int> > prev_neighs(M);
    
    
    // Cluster size distribution
    vector<double> cluster_size_dist(M+1,0);

    
    // AT EACH TIME INSTANT ...
    for(int ti = 0; ti < T; ti++){                  // Time loop
        
        
        // Build KD tree
        const int MM = 9*M;
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
            
            vector<int> nn(num_neigh);              // Nearest neighbour vector
            vector<double> dn(num_neigh);           // Nearest neighbour distances
            
            // Create the neighbor list
            kdtree.nnearest(m,nn,dn,num_neigh);
            
            for (int i = 0; i < num_neigh; i++) {
                int out_flag = 0;
                int j = nn[i];                      // cell_vec index
                int jj = j % M;                     // cell index (for the neighbour cell)
                
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

        
        // CLUSTER SIZE DISTRIBUTION ANALYSIS
        
        // Set up the graph
        Graph graph( M, th );

        
        // Build the graph
        for (int m = 0; m < M; m++) {                   // Cell loop
            
            // Add the neighbors of the cell as edges
            // So nodes are the cells and edges are the neighbors, pretty straightforward
            set<int> destination = cells[m].neighs;
            graph.addEdge( m, destination );
            
        }                                               // Cell loop
        
        
        // Search the edges to find the connected components (or clusters, if you like)
        vector<bool> visited(M);
        for (int m = 0; m < M; m++) { visited[m] = false; }
        stack<int> nodes;
        
        // For all the nodes ...
        char buffer[30];
        sprintf(buffer, "%d", ti);
        char connected_comp_list_out_path[200] = "img/dat/cluster_list_";
        strcat(connected_comp_list_out_path, buffer);
        char buffer_txt[30] = ".txt";
        strcat(connected_comp_list_out_path, buffer_txt);
        
        ofstream connected_comp_list_out(connected_comp_list_out_path);
        
        char connected_comp_size_out_path[200] = "img/dat/cluster_size_";
        strcat(connected_comp_size_out_path, buffer);
        strcat(connected_comp_size_out_path, buffer_txt);
        
        ofstream connected_comp_size_out(connected_comp_size_out_path);
        
        for (int m = 0; m < M; m++) {
            
            // ... if they are not visited ...
            if ( !visited[m] ) {
                //connected_comp_list_out << m << "\t";
                int size_value = dfsFromNode(graph, m, nodes, visited, connected_comp_list_out);
                connected_comp_list_out << endl;
                cluster_size_dist[size_value]++;
                connected_comp_size_out << size_value << endl;
            }
        }
        
//        
//        // Analyze the number of neighbor changes
//        
//        if (ti > 0) {
//            double neigh_chng_cnt_in_bulk = 0.;
//            double neigh_chng_cnt = 0.;
//            for(int m = 0; m < M; m++){
//                set<int> neighs_m = cells[m].neighs;
//                int size_value_m = neighs_m.size();
//                set<int> neighs_m_prev = prev_neighs[m];
//                int size_value_m_prev = neighs_m_prev.size();
//                vector<int> pd;
//                set_difference( neighs_m.begin(), neighs_m.end(), neighs_m_prev.begin(), neighs_m_prev.end(), back_inserter(pd) );
//                if (size_value_m > boundary_identifier_value && size_value_m_prev > boundary_identifier_value) {
//                    neigh_chng_cnt_in_bulk += pd.size();
//                }
//                else{
//                    neigh_chng_cnt += pd.size();
//                }
//            }
//            neigh_chng_cnt_avg_in_bulk += neigh_chng_cnt_in_bulk;
//            neigh_chng_cnt_avg += neigh_chng_cnt;
//        }
        
        
        // Save the previous time neighbour values
        // and DELETE THE NEIGHBOR LIST
        for(int m = 0; m < M; m++){
//            prev_neighs[m] = cells[m].neighs;
            cells[m].neighs.clear();
        }
        

        //
        // NO FURTHER ANALYSIS ON NEIGHBORS BEYOND THIS POINT
        //
        
        
    }           // Time loop
    
    
    // Save the data
    

    
    // Cluster size distribution
    ofstream cluster_size_dist_out("cluster_size_dist.txt");
    for (int m = 0; m < M+1; m++) {
        cluster_size_dist_out << m << "\t" << cluster_size_dist[m]/T << endl;
    }
    
    
//    // Average number of neighbour changes in the bulk
//    ofstream neigh_chng_in_bulk_out;
//    neigh_chng_in_bulk_out.open( "avg_neigh_change_in_bulk.txt", ios::out | ios::app );
//    neigh_chng_in_bulk_out << packing << "\t" << Fmc << "\t" << eps << "\t" << neigh_chng_cnt_avg_in_bulk/(2*(T-1)) << endl;
//    neigh_chng_in_bulk_out.close();
//    
//    
//    // Average number of neighbour changes
//    ofstream neigh_chng_out;
//    neigh_chng_out.open( "avg_neigh_change.txt", ios::out | ios::app );
//    neigh_chng_out << packing << "\t" << Fmc << "\t" << eps << "\t" << neigh_chng_cnt_avg/(2*(T-1)) << endl;
//    neigh_chng_out.close();


    
} // End of the program
