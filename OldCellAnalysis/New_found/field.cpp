
//+++++++++++++++
//++ Data analysis with postprocessing for velocity field calculation
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
    
    // +++++++++++
    
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++ Do the analysis
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    
    int psin = 4;
    double cell_diameter_2 = 4*R_avg*R_avg;
    char buffer_txt[30] = ".txt";
    int num_neigh = 12;
    int up_limit = (int)M/5;

    
    // AT EACH TIME INSTANT ...
    for(int ti = 0; ti < T; ti++){                  // Time loop
        
        
        // Build KD tree
        const int MM = 9*M;
        vector<Point<2> > cell_vec(MM);
        
        
        // Average speed per cell calculation
        vector<double> speed(M);
        
        
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
            }
            
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
        
        
        // Paths
        char buffer[30];
        sprintf(buffer, "%d", ti);
        
        char msd_out_path[200] = "img/dat/msd_field_";
        strcat(msd_out_path, buffer);
        strcat(msd_out_path, buffer_txt);
        ofstream msd_out( msd_out_path );
        
        char bond_out_path[200] = "img/dat/bond_field_";
        strcat(bond_out_path, buffer);
        strcat(bond_out_path, buffer_txt);
        ofstream bond_out( bond_out_path );
        
        char neigh_out_path[200] = "img/dat/neighs_";
        strcat(neigh_out_path, buffer);
        strcat(neigh_out_path, buffer_txt);
        ofstream neigh_out( neigh_out_path );
    

        for(int m = 0; m < M; m++){       // Cell loop
            
            // NEIGHBOUR ANALYSIS
            set<int> neighs = cells[m].neighs;
            neigh_out << m << "\t" << neighs.size() << endl;
            
            // LOCAL MSD
            double delta_x = cells[m].xpos[ti] - cells[m].xpos[0];
            double delta_y = cells[m].ypos[ti] - cells[m].ypos[0];
            msd_out << m << "\t" << (delta_x*delta_x + delta_y*delta_y)/cell_diameter_2 << endl;
            
            double real_angle_per_cell = 0.;
            double img_angle_per_cell = 0.;
            
            vector<int> nn(psin);
            vector<double> dn(psin);
            kdtree.nnearest(m,nn,dn,psin);
            
            // LOCAL BOND ORDER PARAMETER
            for(int j = 0; j < psin; j++){
                double dx = dist_xy(cell_vec[nn[j]], cell_vec[m], 0);
                double dy = dist_xy(cell_vec[nn[j]], cell_vec[m], 1);
                double angle = psin*atan2( dy, dx );
                double cos_angle = cos(angle);
                double sin_angle = sin(angle);
                
                real_angle_per_cell += cos_angle;
                img_angle_per_cell += sin_angle;
            }
            
            bond_out << m << "\t" << sqrt( real_angle_per_cell*real_angle_per_cell + img_angle_per_cell*img_angle_per_cell ) / psin << endl;
            
        }           // Cell loop
        
        for (int m = 0; m < M; m++) {
            cells[m].neighs.clear();
        }
        
        
        // FASTEST CELLS
    
        if (ti > 0) {
            
            char fastest_cells_out_path[200] = "img/dat/fastest_cells_";
            strcat(fastest_cells_out_path, buffer);
            strcat(fastest_cells_out_path, buffer_txt);
            ofstream fastest_cells_out( fastest_cells_out_path );
            
            // Choose the fastest cells
            int break_cnt = 0;
            set<int> idx;
            vector<int> a = sort_indexes(speed);
            for(vector<int>::iterator it = a.begin(); it != a.end(); it++){
                fastest_cells_out << ti*dtSamp << "\t" << *it << "\t" << speed[*it] << endl;
                idx.insert( *it );
                break_cnt++;
                if( break_cnt == up_limit ){
                    break;
                }
            }
            
            char all_cells_out_path[200] = "img/dat/cells_";
            strcat(all_cells_out_path, buffer);
            strcat(all_cells_out_path,buffer_txt);
            ofstream all_cells_out( all_cells_out_path );
            
            for (int m = 0; m < M; m++) {
                all_cells_out << ti*dtSamp << "\t" << m << "\t" << speed[m] << "\t" << cells[m].xi << "\t" << cells[m].yi << "\t" << cells[m].xvel[ti-1] << "\t" << cells[m].yvel[ti-1] << "\t" << cells[m].xpos[ti] << "\t" << cells[m].ypos[ti] << "\n";
            }

        }
        
    }           // Time loop
    
} // End of the program
