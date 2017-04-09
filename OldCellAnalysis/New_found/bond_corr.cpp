

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
    
    
    // Spatial correlations
    double denom_time_avg = 0.;
    double max_dist = 0.;
    int sug_max_dist = (int)fmax( Lx, Ly );
    vector<double> cbb(sug_max_dist);
    for(int i = 0; i < sug_max_dist; i++){
        cbb[i] = 0.;
    }
    
    
    int num_neigh = 12;                         // Number of neighbours to check to make sure they are touching
    int boundary_identifier_value = 4;          // If the cell has less than or equal to neighbors w.r.t. this value, it is marked as a boundary cell
    
    
    // Average number of neighbors per cell
    double avg_neigh_time = 0.;
    
    
    // Average number of neighbors per time frame
    ofstream avg_num_neigh_per_frame_out("num_neigh_per_frame.txt");
    
    
    // Average neighbor change
    double neigh_chng_cnt_avg = 0.;
    double neigh_chng_cnt_avg_in_bulk = 0.;
    vector<set<int> > prev_neighs(M);
    
    
    // Average number of T1 transitions
    double t1_cnt_avg = 0.;
    vector<T1set> t1set_prev;
    
    
    // Average number of neighbors histogram
    int max_num_neigh_size = 11;
    vector<double> avg_num_neigh(max_num_neigh_size,0);
    
    
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
        
        
        
        // T1 transition analysis
        // depends on the candidate list from the previous time step
        
        if (ti > 0) {
            
            // For each pair forming a T1 set ...
            for(int i = 0; i < t1set_prev.size(); i++){
                
                // Make sure previously neighbored points are separated (they shouldn't be neighbors!)
                bool flag = true;
                set<int> neighs_n1 = cells[t1set_prev[i].n1].neighs;
                for(set<int>::iterator it = neighs_n1.begin(); it != neighs_n1.end(); it++){
                    if (*it == t1set_prev[i].n2) {   // If they are neighbors, raise the flag to false and break the for loop
                        flag = false;
                        break;
                    }
                }
                
                // If n1 and n2 are separated
                if (flag == true) {
                    
                    // Check if the previously separate points are now neighbors or not (they should be!)
                    bool flag2 = false;
                    set<int> neighs_non1 = cells[t1set_prev[i].non1].neighs;
                    for(set<int>::iterator it = neighs_non1.begin(); it != neighs_non1.end(); it++){
                        if (*it == t1set_prev[i].non2) {    // If they are indeed neighbors, raise the flag to true
                            flag2 = true;
                        }
                    }
                    
                    // And, if their neighborhood contains 2 identical points or not
                    if (flag2 == true) {
                        vector<int> intersect;
                        set<int> neighs_non2 = cells[t1set_prev[i].non2].neighs;
                        set_intersection( neighs_non1.begin(), neighs_non1.end(), neighs_non2.begin(), neighs_non2.end(), back_inserter(intersect) );
                        
                        // It is a T1 transition!
                        if (intersect.size() == 2) {
                            t1_cnt_avg++;
                        }
                    }
                    
                    
                }
                
                
            }   // for each pair in T1
        }
        
        
        // Build the possible T1 candidate list for further analysis in the next step
        vector<T1set> t1set;
        
        for(int m = 0; m < M; m++){                                                             // Vertices
            
            set<int> neighs_m = cells[m].neighs;                                                // Neighbors of the vertex
            
            for(set<int>::iterator it = neighs_m.begin(); it != neighs_m.end(); it++) {         // Run over the neighbors of the vertex
                
                T1set ta;
                
                vector<int> intersect;                                                          // Intersection of neighbors of the vertex with neighbors of the neighbor
                set<int> neighs_i = cells[*it].neighs;                                          // Neighbors of the neighbor
                
                // Calculate the intersection
                set_intersection( neighs_m.begin(), neighs_m.end(), neighs_i.begin(), neighs_i.end(), back_inserter(intersect) );
                
                // Add the intersecting points to the T1 set
                // They are forming a T1 set as for every point, each has a neighborhood not containing the other
                if (intersect.size() == 2) {
                    int icnt = 0;
                    for(vector<int>::iterator jt = intersect.begin(); jt != intersect.end(); jt++){
                        icnt++;
                        ta.addPair(*jt,icnt);
                    }         // jt -- intersection with neighbors and neighbors of neighbor
                    
                    // Make sure that these identical neighbors don't contain each other in their respective neighborhood at this point in time
                    bool flag = true;
                    set<int> neighs_1 = cells[ta.non1].neighs;                                   // Neighbor set
                    for(set<int>::iterator jt = neighs_1.begin(); jt != neighs_1.end(); jt++){
                        if (*jt == ta.non2) {           // If neighbor equals to the other point, they are not a candidate for T1 transition
                            flag = false;
                        }
                    }
                    
                    // If the points are qualified as a T1 set, store them for further checking in the next step in time
                    if (flag == true) {
                        ta.addElem(m, *it);
                        t1set.push_back( ta );
                    }
                    
                }
                
            }               // it -- neighbors of vertex
            
        }                   // m -- vertex
        
        // Save the candidate list
        t1set_prev = t1set;

        
        // BOND ORDER PARAMETER and NEIGHBORHOOD ANALYSIS
        
        // Parameters for bond order
        double real_angle = 0.;
        double img_angle = 0.;
        vector<double> psin_per_cell(M, 0);
        
        // Parameters for average number of neighbors
        double avg_neigh_ens = 0.;
        
        for(int m = 0; m < M; m++){       // Cell loop
            
            set<int> neighs = cells[m].neighs;
            int size_value = neighs.size();
            if (size_value > boundary_identifier_value) {
                
                // Calculate average number of neighbors per cell
                avg_neigh_ens += size_value;
                
            }
            
            // Calculate average number of neighbor histogram
            avg_num_neigh[size_value]++;
            
            // Calculate bond order parameter
            double real_angle_per_cell = 0.;
            double img_angle_per_cell = 0.;
            
            vector<int> nn(psin);
            vector<double> dn(psin);
            kdtree.nnearest(m,nn,dn,psin);
            
            for(int j = 0; j < psin; j++){
                double dx = dist_xy(cell_vec[nn[j]], cell_vec[m], 0);
                double dy = dist_xy(cell_vec[nn[j]], cell_vec[m], 1);
                double angle = psin*atan2( dy, dx );
                double cos_angle = cos(angle);
                double sin_angle = sin(angle);
                
                real_angle += cos_angle;
                img_angle += sin_angle;
                
                real_angle_per_cell += cos_angle;
                img_angle_per_cell += sin_angle;
            }
            
            psin_per_cell[m] = sqrt( real_angle_per_cell*real_angle_per_cell + img_angle_per_cell*img_angle_per_cell ) / psin;
            
        }           // Cell loop
        
        // Calculate bond order parameter and susceptibility
        double tmp_value = sqrt( real_angle*real_angle + img_angle*img_angle )/(M*psin);
        psin_time_avg += tmp_value;
        susc_sq_avg += tmp_value*tmp_value;
        
        // Calculate average number of neighbors
        tmp_value = avg_neigh_ens/M;
        avg_neigh_time += tmp_value;
        avg_num_neigh_per_frame_out << ti << "\t" << tmp_value << endl;
        
        
        
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
        for (int m = 0; m < M; m++) {
            
            // ... if they are not visited ...
            if ( !visited[m] ) {
                int size_value = dfsFromNode(graph, m, nodes, visited);
                cluster_size_dist[size_value]++;
            }
        }
        
        
        
        
        // CORRELATION IN SPACE ...
        // Spatial bond order correlation
        
        // Spatial bond correlation parameters
        vector<double> cbb_avg(sug_max_dist);
        vector<double> c_cnt(sug_max_dist);
        for(int i = 0; i < sug_max_dist; i++){
            cbb_avg[i] = 0.;
            c_cnt[i] = 0.;
        }
        
        for(int m1 = 0; m1 < M; m1++){              // First cell loop
            for(int m2 = 0; m2 < M; m2++){          // Second cell loop
                
                double dx = cells[m2].xpos[ti] - cells[m1].xpos[ti];
                dx -= Lx*dnearbyint( dx/Lx );
                double dy = cells[m2].ypos[ti] - cells[m1].ypos[ti];
                dy -= Ly*dnearbyint( dy/Ly );
                double dr = sqrt( dx*dx + dy*dy );
                
                    int bin = inearbyint( dr );
                    c_cnt[bin]++;
                    
                    // Spatial bond correlation
                    cbb_avg[bin] += psin_per_cell[m2]*psin_per_cell[m1];
                    
                    // Keep track of maximum distance between the cells
                    if( dr > max_dist ){
                        max_dist = dr;
                    }

                
            }       // First cell loop
        }          // Second cell loop
        
        // Ensemble averages for correlation in space
        for(int i = 0; i < sug_max_dist; i++){
            if( c_cnt[i] != 0 ){
                cbb[i] += cbb_avg[i]/c_cnt[i];
            }
        }

        
        
        // Analyze the number of neighbor changes
        
        if (ti > 0) {
            double neigh_chng_cnt_in_bulk = 0.;
            double neigh_chng_cnt = 0.;
            for(int m = 0; m < M; m++){
                set<int> neighs_m = cells[m].neighs;
                int size_value_m = neighs_m.size();
                set<int> neighs_m_prev = prev_neighs[m];
                int size_value_m_prev = neighs_m_prev.size();
                vector<int> pd;
                set_difference( neighs_m.begin(), neighs_m.end(), neighs_m_prev.begin(), neighs_m_prev.end(), back_inserter(pd) );
                if (size_value_m > boundary_identifier_value && size_value_m_prev > boundary_identifier_value) {
                    neigh_chng_cnt_in_bulk += pd.size();
                }
                else{
                    neigh_chng_cnt += pd.size();
                }
            }
            neigh_chng_cnt_avg_in_bulk += neigh_chng_cnt_in_bulk;
            neigh_chng_cnt_avg += neigh_chng_cnt;
        }
        
        
        // Save the previous time neighbour values
        // and DELETE THE NEIGHBOR LIST
        for(int m = 0; m < M; m++){
            prev_neighs[m] = cells[m].neighs;
            cells[m].neighs.clear();
        }
        

        //
        // NO FURTHER ANALYSIS ON NEIGHBORS BEYOND THIS POINT
        //
        
        
    }           // Time loop
    
    
    // Save the data
    
    
    // Bond order parameter
    psin_time_avg /= T;
    ofstream psin_out;
    psin_out.open( "glob_bond.txt", ios::out | ios::app );
    psin_out << packing << "\t" << Fmc << "\t" << eps << "\t" << psin_time_avg << endl;
    
    
    // Bond order susceptibility
    susc_sq_avg /= T;
    ofstream susc_out;
    susc_out.open( "susceptibility.txt", ios::out | ios::app );
    susc_out << packing << "\t" << Fmc << "\t" << eps << "\t" <<  M*(susc_sq_avg - psin_time_avg*psin_time_avg) << endl;
    
    
    // Average number of neighbours per cell
    ofstream avg_neigh_out;
    avg_neigh_out.open( "avg_neigh_per_cell.txt", ios::out | ios::app );
    avg_neigh_out << packing << "\t" << Fmc << "\t" << eps << "\t" << avg_neigh_time/(T-1) << endl;
    avg_neigh_out.close();
    
    
    // Cluster size distribution
    ofstream cluster_size_dist_out("cluster_size_dist.txt");
    for (int m = 0; m < M+1; m++) {
        cluster_size_dist_out << m << "\t" << cluster_size_dist[m]/T << endl;
    }
    
    
    // Average number of neighbour histogram
    ofstream avg_num_neigh_out("num_neigh_dist.txt");
    for (int i = 0; i < max_num_neigh_size; i++) {
        avg_num_neigh_out << i << "\t" << avg_num_neigh[i]/T << endl;
    }
    
    
    // Average number of neighbour changes in the bulk
    ofstream neigh_chng_in_bulk_out;
    neigh_chng_in_bulk_out.open( "avg_neigh_change_in_bulk.txt", ios::out | ios::app );
    neigh_chng_in_bulk_out << packing << "\t" << Fmc << "\t" << eps << "\t" << neigh_chng_cnt_avg_in_bulk/(2*(T-1)) << endl;
    neigh_chng_in_bulk_out.close();
    
    
    // Average number of neighbour changes
    ofstream neigh_chng_out;
    neigh_chng_out.open( "avg_neigh_change.txt", ios::out | ios::app );
    neigh_chng_out << packing << "\t" << Fmc << "\t" << eps << "\t" << neigh_chng_cnt_avg/(2*(T-1)) << endl;
    neigh_chng_out.close();

    
    
    // Average number of T1 transitions
    ofstream t1_out;
    t1_out.open( "avg_t1_transition.txt", ios::out | ios::app);
    t1_out << packing << "\t" << Fmc << "\t" << eps << "\t" << t1_cnt_avg/(T-1) << endl;
    t1_out.close();
    
    
    // Spatial correlations
    int total_bins = (int)max_dist;
    cout << "Maximum distance observed is " << max_dist << endl;
    ofstream bond_corr_out("bond_corr.txt");
    for(int i = 0; i < total_bins; i++){
        bond_corr_out << i/(2*R_avg) << "\t" << cbb[i]/T << endl;
    }
    
    
} // End of the program
