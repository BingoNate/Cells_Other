
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
    
    char * input_file = argv[1];                  // Filename for the particle position data in xyz file format
    
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
    
    
    // Average speed per cell calculation
    double speed_per_cell_time_avg = 0.;
    
    
    // Mean positions in each time frame
    vector<double> x_mean(T,0);
    vector<double> y_mean(T,0);
    
    
    // Mean velocities in each time frame
    vector<double> vx_mean(T-1,0);
    vector<double> vy_mean(T-1,0);
    
    
    // Bond order parameter
    int psin = 4;
    double psin_time_avg = 0.;
    double susc_sq_avg = 0.;
    
    
    // Average number of neighbors per cell
    double avg_neigh_time = 0.;
    double neigh_chng_cnt_avg = 0.;
    vector<set<int> > prev_neighs(M);
    int num_neigh = 12;                     // Number of neighbours to check
    double t1_cnt_avg = 0.;                 // Average number of T1 transitions
    vector<T1set> t1set_prev;
    
    
    // Spatial correlations
    double denom_time_avg = 0.;
    double max_dist = 0.;
    int sug_max_dist = (int)fmax( Lx, Ly );
    vector<double> cbb(sug_max_dist);
    vector<double> cnn(sug_max_dist);
    vector<double> cvv(sug_max_dist);
    for(int i = 0; i < sug_max_dist; i++){
        cbb[i] = 0.;
        cvv[i] = 0.;
        cnn[i] = 0.;
    }
    
    
    // Static correlations
    int maxk = (int)(fmax(Lx,Ly)/2);
    vector<double> S_x_time_avg(maxk, 0);
    vector<double> S_y_time_avg(maxk, 0);
    double vir_time_avg = 0.;
    double sig2 = sig*sig;
    double eps48 = 48*eps;
    double max_cluster_size_avg = 0.;
    double th = 2*R_avg + Rcut;

    
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
            
            // Mean positions in a frame
            x_mean[ti] += cells[m].xpos[ti];
            y_mean[ti] += cells[m].ypos[ti];
            
        }           // Cell loop
    
        
        // Build the KD tree
        KDtree<2> kdtree(cell_vec);
        
        // Neighbourhood check to be absolutely sure
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
        
        
        // T1 transition analysis
        
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
        
        t1set_prev = t1set;


        
        // Mean positions in a frame
        x_mean[ti] /= M;
        y_mean[ti] /= M;
        
        // Mean velocities in a frame
        if( ti > 0 ){
            vx_mean[ti-1] /= M;
            vy_mean[ti-1] /= M;
        }
        
        
        // Average speed per cell calculation
        speed_per_cell_time_avg += speed_per_cell_ens_avg / M;
        
        
        // BOND ORDER PARAMETER and NEIGHBORHOOD ANALYSIS
        double real_angle = 0.;
        double img_angle = 0.;
        vector<double> psin_per_cell(M, 0);
        
        double avg_neigh_ens = 0.;
        
        for(int m = 0; m < M; m++){       // Cell loop
            
            set<int> neighs = cells[m].neighs;
            avg_neigh_ens += neighs.size();
            
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
        
        double tmp_value = sqrt( real_angle*real_angle + img_angle*img_angle )/(M*psin);
        psin_time_avg += tmp_value;
        susc_sq_avg += tmp_value*tmp_value;
        
        avg_neigh_time += avg_neigh_ens/M;
        
        // CORRELATION IN SPACE ...
        // Spatial bond order correlation,
        // pair correlation, and
        // spatial velocity correlation
        
        vector<double> cbb_avg(sug_max_dist);
        vector<double> cnn_avg(sug_max_dist);
        vector<double> cvv_avg(sug_max_dist);
        vector<double> c_cnt(sug_max_dist);
        vector<double> cvv_cnt(sug_max_dist);
        for(int i = 0; i < sug_max_dist; i++){
            cbb_avg[i] = 0.;
            cvv_avg[i] = 0.;
            cnn_avg[i] = 0.;
            c_cnt[i] = 0.;
            cvv_cnt[i] = 0.;
        }
        double denom_ens_avg = 0.;
        double denom_cnt = 0.;
        
        /// +++++++++++

        
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
                        
                        // Spatial bond correlation
                        cbb_avg[bin] += psin_per_cell[m2]*psin_per_cell[m1];
                        
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
                    
                }       // First cell loop
            }          // Second cell loop

            S_x_time_avg[ki] += (S_cos_x*S_cos_x + S_sin_x*S_sin_x)/S_cnt;
            S_y_time_avg[ki] += (S_cos_y*S_cos_y + S_sin_y*S_sin_y)/S_cnt;
            
        }               // Wavevector loop
        
        /// +++++++++++
        


        // Ensemble averages for correlation in space
        for(int i = 0; i < sug_max_dist; i++){
            if( c_cnt[i] != 0 ){
                cbb[i] += cbb_avg[i]/c_cnt[i];
            }
            if( ti > 0 && cvv_cnt[i] != 0 ){
                cvv[i] += cvv_avg[i]/cvv_cnt[i];
            }
            cnn[i] += cnn_avg[i];
        }
        
        if( ti > 0 && denom_cnt != 0 ){
            denom_time_avg += denom_ens_avg/denom_cnt;
        }
        
        
        // STATIC CORRELATIONS ...
        // Interparticle virial stress
        // Maximum cluster size
        

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
        
        
        // Maximum cluster size
        // Set up the graph
        Graph graph((int)M/5, th);
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
        
        for(set<int>::iterator it = idx.begin(); it != idx.end(); it++){  // loop over i
            
            set<int> destination;
            for(set<int>::iterator jt = idx.begin(); jt != idx.end(); jt++){  // loop over j
                if( *it != *jt ){
                    double dx = cells[*jt].xpos[ti] - cells[*it].xpos[ti];
                    dx -= Lx*dnearbyint( dx/Lx );
                    double dy = cells[*jt].ypos[ti] - cells[*it].ypos[ti];
                    dy -= Ly*dnearbyint( dy/Ly );
                    double d2 = dx*dx + dy*dy;
                    if( d2 < threshold2 ) { destination.insert( *jt ); }
                }
            }  // loop over j
            
            graph.addEdge( *it, destination );
            
        }  // loop over i

        // Search the fastest cells to find the connected components
        vector<int> max_cluster_size;
        vector<bool> visited(M);
        for(set<int>::iterator it = idx.begin(); it != idx.end(); it++){ visited[*it] = false; }
        stack<int> nodes;
        
        // For all the fastest nodes ...
        for(set<int>::iterator it = idx.begin(); it != idx.end(); it++){
            
            // If they are not visited ...
            if( !visited[*it] ){
                max_cluster_size.push_back( dfsFromNode(graph, *it, nodes, visited) );
            }
        }
        
        max_cluster_size_avg += *max_element( max_cluster_size.begin(), max_cluster_size.end() );
        
        // Analyze the number of neighbor changes
        if (ti > 0) {
            double neigh_chng_cnt = 0.;
            for(int m = 0; m < M; m++){
                set<int> neighs_m = cells[m].neighs;
                set<int> neighs_m_prev = prev_neighs[m];
                vector<int> pd;
                set_difference( neighs_m.begin(), neighs_m.end(), neighs_m_prev.begin(), neighs_m_prev.end(), back_inserter(pd) );
                neigh_chng_cnt += pd.size();
            }
            neigh_chng_cnt_avg += neigh_chng_cnt;
        }
        
        
        // Save the previous time neighbour values
        for(int m = 0; m < M; m++){
            prev_neighs[m] = cells[m].neighs;
            cells[m].neighs.clear();
        }
        
                
    }               // Time loop
    
    
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
    
    
    // Average speed per cell calculation
    ofstream speed_per_cell_out;
    speed_per_cell_out.open( "avg_speed_per_cell.txt", ios::out | ios::app );
    speed_per_cell_out << packing << "\t" << Fmc << "\t" << eps << "\t" << speed_per_cell_time_avg / (T-1) << endl;
    speed_per_cell_out.close();
    
    
    // Spatial correlations
    int total_bins = (int)max_dist;
    cout << "Maximum distance observed is " << max_dist << endl;
    ofstream bond_corr_out("bond_corr.txt");
    ofstream pair_out("pair_corr.txt");
    ofstream sp_vacf_out("sp_vacf.txt");
    double rho = M/(Lx*Ly);
    denom_time_avg /= (T-1);
    for(int i = 0; i < total_bins; i++){
        bond_corr_out << i/(2*R_avg) << "\t" << cbb[i]/T << endl;
        
        if( i != 0 ){
            double n_ideal = 2*pi*i*rho;
            pair_out << i/(2*R_avg) << "\t" << cnn[i]/(n_ideal*T*M) << endl;
        }
        
        sp_vacf_out << i/(2*R_avg) << "\t" << cvv[i]/(denom_time_avg*(T-1)) << endl;
    }
    
    
    // Static structure factor
    ofstream static_struct_out("static_struct.txt");
    for(int ki = 1; ki < maxk; ki++){
        double k = 2*pi*ki/maxk;
        S_x_time_avg[ki] /= T;
        S_y_time_avg[ki] /= T;
        static_struct_out << k*2*R_avg << "\t" << (S_x_time_avg[ki]+S_y_time_avg[ki])/(2*M) << endl;
    }
    
    
    // Interparticle virial stress
    ofstream glob_vir_out;
    glob_vir_out.open( "glob_vir_stress.txt", ios::out | ios::app );
    glob_vir_out << packing << "\t" << Fmc << "\t" << eps << "\t" << vir_time_avg/T << endl;
    glob_vir_out.close();
    
    
    // Maximum cluster size
    ofstream max_cluster_size_out;
    max_cluster_size_out.open( "max_cluster_size.txt", ios::out | ios::app );
    max_cluster_size_out << packing << "\t" << Fmc << "\t" << eps << "\t" << max_cluster_size_avg/T << endl;
    max_cluster_size_out.close();
    
    
    // DYNAMIC CORRELATIONS ...
    // Velocity autocorrelation
    // Intermediate scattering function
    // Dynamic structure factor
    // Mean square displacement
    
    ofstream vacf_out("vacf_cm.txt");
    ofstream inter_scattering_out("inter_scattering.txt");
    double k_int_scatter = 0.5;
    vector<double> Fr(D);
    vector<double> Fi(D);
    int maxw = (int)D/2;
    ofstream msd_out("msd_post.txt");
    double cell_diam_avg = 4*R_avg*R_avg;
    
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
        
        double Fr_x = 0.;
        double Fr_y = 0.;
        double Fi_x = 0.;
        double Fi_y = 0.;
        
        double msd_ens_avg = 0.;
        
        for(int m = 0; m < M; m++){             // Cell loop
            
            double cv_time_avg = 0.;
            
            double F_cos_x = 0.;
            double F_sin_x = 0.;
            double F_cos_y = 0.;
            double F_sin_y = 0.;
            
            double msd_time_avg = 0.;
            
            for(int ti = 0; ti < TD; ti++){        // Different time origins
                
                // Velocity autocorrelation
                double vx_1 = cells[m].xvel[ti];
                double vx_2 = cells[m].xvel[ti+di];
                double vy_1 = cells[m].yvel[ti];
                double vy_2 = cells[m].yvel[ti+di];
                
                cv_time_avg += vx_1*vx_2 + vy_1*vy_2;
                
                // Density autocorrelation in Fourier space (intermediate scattering function)
                double dx = cells[m].xpos[ti+di] - cells[m].xpos[ti];
                double dy = cells[m].ypos[ti+di] - cells[m].ypos[ti];
                
                double tut = k_int_scatter*dx;
                F_cos_x += cos( tut );
                F_sin_x += sin( tut );
                
                tut = k_int_scatter*dy;
                F_cos_y += cos( tut );
                F_sin_y += sin( tut );
                
                // Mean square displacement
                double x_1 = cells[m].xpos[ti] - x_mean[ti];
                double x_2 = cells[m].xpos[ti+di] - x_mean[ti];
                double y_1 = cells[m].ypos[ti] - y_mean[ti];
                double y_2 = cells[m].ypos[ti+di] - y_mean[ti];
                
                double delta_x = x_2 - x_1;
                double delta_y = y_2 - y_1;
                
                msd_time_avg += (delta_x*delta_x + delta_y*delta_y)/cell_diam_avg;
                
            }       // Different time origins
            
            cv_ens_avg += cv_time_avg/TD;
            
            Fr_x += F_cos_x/TD;
            Fr_y += F_cos_y/TD;
            Fi_x += F_sin_x/TD;
            Fi_y += F_sin_y/TD;
            
            msd_ens_avg += msd_time_avg/TD;
            
        }           // Cell loop
        
        // Velocity correlation
        double time_stamp = dtSamp*di;
        vacf_out << time_stamp << "\t" << cv_ens_avg/(M*denom_ens_avg) << endl;
        
        // Intermediate scattering
        Fr[di] = (Fr_x+Fr_y)/(2*M);
        Fi[di] = (Fi_x+Fi_y)/(2*M);
        
        if( di != 0 ){
            inter_scattering_out << time_stamp << "\t" << Fr[di] << endl;
        }
        else{
            inter_scattering_out << 0 << "\t" << 1 << endl;
        }
        
        // Mean square displacement
        msd_out << time_stamp << "\t" << msd_ens_avg/M << endl;
        
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
    
}




