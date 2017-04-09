
//+++++++++++++++
//++ Data analysis with postprocessing for velocity correlators
//+++++++++++++++

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <map>
#include <set>
#include <stack>
#include <algorithm>
#include <boost/multi_array.hpp>
#include "variables.h"
#include "basic.h"

using namespace std;

#define pi M_PI

#define VEL
#define VACF
#define SP_VACF
#define AVG_VEL
#define VEL_DIR
#define VEL_SEG

// +++++++
// DATA TYPES
// +++++++

// +++++
// Cell
// +++++

class Cell
{
    
public:
    vector<double> xpos;
    vector<double> ypos;
    vector<double> xvel;
    vector<double> yvel;
    int numCell;
    
    Cell() {}
    ~Cell() {}
    
};

// +++++


// +++++
// The periodic box
// +++++

class Box
{
    
public:
    Box() {}
    ~Box() {}
    void setWidth( double wdth );
    void setHeight( double hght );
    double getWidth() { return width; }
    double getHeight() { return height; }
    double getArea() { return width*height; }
    
private:
    double width;
    double height;
    
};

void Box::setWidth( double wdth ){
    width = wdth;
}

void Box::setHeight( double hght ){
    height = hght;
}

// +++++


// +++++
// Simulation details
// +++++

class Data
{
    
public:
    Data() {}
    ~Data() {}
    void setNumStep( double nst );
    void setNumDelay( double nd );
    void setNumTotalDelay();
    double getNumStep() { return numStep; }
    double getNumDelay() { return numDelay; }
    double getNumTotalDelay() { return numTotalDelay; }
    
private:
    double numStep;
    double numDelay;
    double numTotalDelay;
    
};

void Data::setNumStep( double nst ){
    numStep = nst;
}

void Data::setNumDelay( double nd ){
    numDelay = nd;
}

void Data::setNumTotalDelay(){
    numTotalDelay = numStep - numDelay;
}

// +++++


// +++++
// Graph
// +++++

// ( Note that edges are lines connecting nodes or vertices )

class Graph
{
    
public:
    Graph() {}
    ~Graph() {}
    void setNumNodes( int nds );
    int getNumNodes() { return numNodes; }
    void addEdge( int nd_idx, set<int> dest );
    map<int,set<int> > edges; // [node, destination]
    void setThreshold( double cutoff, double carpan );
    double getThreshold() { return threshold; }
    
private:
    int numNodes;
    double threshold;
    
};

void Graph::setNumNodes( int nds ){
    numNodes = nds;
}

void Graph::addEdge( int nd_idx, set<int> dest ){
    edges.insert( pair<int,set<int> >( nd_idx, dest ) );
}

void Graph::setThreshold( double cutoff, double carpan ){
    threshold = 2*R_avg + carpan*cutoff;
}


// +++++


// Data type for holding multidimensional data
typedef boost::multi_array<double,2> array_type;


// Function for reading the positions
void readData( vector<Cell> & cells, char * file_name, char xory, int& sample_cnt, int sample_data ){
    
    // Inputs
    int cell_cnt = -1;
    ifstream inputData;
    inputData.open( file_name );
    
    // Holder for each line
    string line_holder;
    
    // Yield error in case of a problem
    if( inputData.fail() ){
        cerr << "The file cannot be opened!\n";
        exit(1);
    }
    
    // Decide whether to read the x value or the value
    switch (xory) {
            
            // For x
        case 'x':
            
            // Read the lines
            while( inputData.good() ){ // If everything is good ...
                
                while( getline(inputData, line_holder) )
                { // Read the lines one by one
                    
                    // Data sampling -optional-
                    sample_cnt++;
                    if( sample_cnt % sample_data == 0 ){
                        
                        // Holder for each element in the read line
                        istringstream stream_holder( line_holder );
                        
                        // Yet another holder for data type of each element
                        double double_holder;
                        
                        // Now read the values in a line, finally ...
                        while ( stream_holder >> double_holder )
                        {
                            
                            // Each element in the line consists of center of mass of cells at that given time instant
                            cell_cnt++;
                            cells[cell_cnt].xpos.push_back( double_holder );
                            
                        } // Read a particular line
                        
                        cell_cnt = -1;
                        
                    } // read line by line
                } // data sampling
            } // while data is good
            
            break;
            
            // For y
        case 'y':
            
            // Read the lines
            while( inputData.good() ){ // If everything is good ...
                
                while( getline(inputData, line_holder) )
                { // Read the lines one by one
                    
                    // Data sampling -optional-
                    sample_cnt++;
                    if( sample_cnt % sample_data == 0 ){
                        
                        // Holder for each elements in the read line
                        istringstream stream_holder( line_holder );
                        
                        // Yet another holder for data type of each element
                        double double_holder;
                        
                        // Now read the values in a line, finally ...
                        while ( stream_holder >> double_holder )
                        {
                            
                            // Each element in the line consists of center of mass of cells at that given time instant
                            cell_cnt++;
                            cells[cell_cnt].ypos.push_back( double_holder );
                            
                        } // Read a particular line
                        
                        cell_cnt = -1;
                        
                    } // read line by line
                } // data sampling
            } // while data is good
            
            break;
            
    } // switch end
    
    // Good, you're done!!
    inputData.close();
    
}


// Function to sort a vector by value while keeping the indexes
template <typename T>
vector<int> sort_indexes( const vector<T> &v ) {

	// Initialize original index locations
	vector<int> idx( v.size() );
	for( int i = 0; i != idx.size(); ++i ) idx[i] = i;

	// Sort indexes by comparing values in v
	sort( idx.begin(), idx.end(), [&v](int i1, int i2) { return v[i1] > v[i2]; } );
	
	return idx;
}


// Function for depth-first search from a starting node
int dfsFromNode( Graph& graph, int i, stack<int>& S, vector<bool>& visited ){
    
    int connected_comp = 0;
    
    // Add the node to the stack
    S.push( i );
    
    // While there are nodes that are not visited
    // (so, as long as stack is not empty ..)
    while( !S.empty() ){
        
        // Remove the top of the stack (backtracking process)
        i = S.top();
        S.pop();
        if( !visited[i] ){
            visited[i] = true;
//            cout << i << "\t";
            connected_comp++;   // Add a connected component when you visit a new node
            set<int> neighbors;
            neighbors = graph.edges[i];
            for( auto j: neighbors ){
                i = j;
                S.push( i );
            }
        } // if the node is visited before, get out
        
    } // while loop to check if the stack is empty or not
    
    return connected_comp;

}

// +++++++


int main( int argc, char *argv[] ){
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++ Get the data
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    // Choose the data to be used
    const int nsamp = samplePostData;                // Sampling interval
    const double dtSamp = nsamp*dt;         // Time units between two data points
    
    // Instantiate the data structure
    vector<Cell> cells( M );
    
    // Set filenames
    char * x_input_file = argv[1];        // Filename for the x data
    char * y_input_file = argv[2];        // Filename for the y data

    // Read the data to the cells
    int sample_cnt = -1;
    int sample_data = 1;
    char getX = 'x';
    readData( cells, x_input_file, getX, sample_cnt, sample_data );
    sample_cnt = -1;
    char getY = 'y';
    readData( cells, y_input_file, getY, sample_cnt, sample_data );
    
    // Set general simulation variables
    Data simData;
    simData.setNumStep( cells[0].xpos.size() );
    simData.setNumDelay( sqrt( cells[0].xpos.size() ) );
    simData.setNumTotalDelay();
    
    const double T = simData.getNumStep();              // Total time
    const double D = simData.getNumDelay();             // Last delay time
    const double TD = simData.getNumTotalDelay();       // Total time - last delay time
    
    // Set the box
    Box box;
    box.setWidth( boxSize_x );
    box.setHeight( boxSize_y );
    
    const double Lx = box.getWidth();
    const double Ly = box.getHeight();
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++ Do the analysis
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++ Velocity correlators
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    #ifdef VEL

    // Calculate the velocities
    for(int m = 0; m < M; m++){
        for(int i = 1; i < T; i++){ // Careful with the starting index!
            cells[m].xvel.push_back( (cells[m].xpos[i]-cells[m].xpos[i-1])/dtSamp );
            cells[m].yvel.push_back( (cells[m].ypos[i]-cells[m].ypos[i-1])/dtSamp );
        }
    }
 
    const int TV = cells[0].xvel.size();
    const int DV = (int)sqrt(TV);
    const int TVD = TV - DV;
 
    // Calculate the average velocity magnitudes at different sampling intervals
    #ifdef AVG_VEL
    double vsq = 0;                                 // Velocity magnitude
    double vavg = 0;                                // Average velocity magnitude
    ofstream avgvel_out;
    avgvel_out.open( "avgvel.txt", ios::out | ios::app );
    for(int m = 0; m < M; m++){
        for(int i = 0; i < TV; i++){
            vsq += sqrt( cells[m].xvel[i]*cells[m].xvel[i] + cells[m].yvel[i]*cells[m].yvel[i] );
        }
        vavg += vsq/TV;
        vsq = 0;
    }
    avgvel_out << dens << "\t" << Fmc << "\t" << eps << "\t" << vavg/M << endl;
    avgvel_out.close();
    #endif // end of AVG_VEL
    
    
    // ++++++++++
    // Velocity autocorrelation
    // ++++++++++
    
    #ifdef VACF
    
    double tmpVACF = 0;
    double vacf = 0;
    ofstream vacf_out("vacf_cm.txt");
    
    // Calculate normalization coefficient
    double origin_avg = 0;
    for(int m = 0; m < M; m++){ // Run over all cells
        for(int i = 0; i < TVD; i++){ // Run over different time origins
            
            double vx_origin = cells[m].xvel[i];
            double vy_origin = cells[m].yvel[i];
            
            tmpVACF += vx_origin*vx_origin + vy_origin*vy_origin;
            
        }
        origin_avg += tmpVACF/TVD; // Average over time
        tmpVACF = 0;
    }
    origin_avg /= M; // Average over cells

    for(int d = 0; d < DV; d++){ // Run over all delay times
        for(int m = 0; m < M; m++){ // Run over all cells
            for(int i = 0; i < TVD; i++){ // Run over different time origins
                
                double vx_origin = cells[m].xvel[i];
                double vy_origin = cells[m].yvel[i];
                double vx_delay = cells[m].xvel[i+d];
                double vy_delay = cells[m].yvel[i+d];
                
                tmpVACF += vx_delay*vx_origin + vy_delay*vy_origin;
                
            }
            
            vacf += tmpVACF/TVD; // Average over time
            tmpVACF = 0;
        }
        vacf_out << dtSamp*d << "\t" << vacf/(M*origin_avg) << endl; // Average over cells and normalize
        vacf = 0;
    }
    
    #endif // end of VACF
    
    
    // ++++++++
    // Spatial velocity correlation
    // ++++++++
    
    #ifdef SP_VACF
    
    double* vx_mean = new double[TV];
    double* vy_mean = new double[TV];
    for( int i = 0; i < TV; i++ ){
        vx_mean[i] = 0;
        vy_mean[i] = 0;
    }
   
    double denom_sum = 0;
    double denom_avg = 0;
    int t_offset = 0;
    int t_offset_2 = T;
    int t_offset_total = t_offset_2 - t_offset;
    
    // Calculate the normalization factor (the denominator in the equation)
    for( int i = t_offset; i < t_offset_2; i++){
        
        // Calculate the mean velocity vector in each time frame
        for(int m = 0; m < M; m++){
            vx_mean[i] += cells[m].xvel[i];
            vy_mean[i] += cells[m].yvel[i];
        }
        vx_mean[i] /= M;
        vy_mean[i] /= M;
        
        // Calculate the normalization factor
        for(int m = 0; m < M; m++){
            double vx_1 = cells[m].xvel[i] - vx_mean[i];
            double vy_1 = cells[m].yvel[i] - vy_mean[i];
            
            denom_sum += vx_1*vx_1 + vy_1*vy_1;
        }
        denom_avg += denom_sum/M;
        denom_sum = 0;
    }
    
    // Normalize the normalization factor
    denom_avg /= t_offset_total;
    
    double max_dist = 0.0;
    int sug_max_dist = (int)fmax( boxSize_x, boxSize_y );  // Suggest a sane approximation for maximum distance
    vector<double> cvv(sug_max_dist);
    vector<double> cvv_n(sug_max_dist);
    #ifdef VEL_DIR
    vector<double> cvv_par(sug_max_dist);
    vector<double> cvv_perp(sug_max_dist);
    #endif
    for( int i = 0; i < sug_max_dist; i++ ){
        cvv[i] = 0.0;
        cvv_n[i] = 0.0;
        #ifdef VEL_DIR
        cvv_par[i] = 0.0;
        cvv_perp[i] = 0.0;
        #endif
    }
    
    // Calculate the numerator in the equation
    for( int i = t_offset; i < t_offset_2; i++ ){             // Run over all the timeseries
        
        for( int m1 = 0; m1 < M; m1++ ){                  // Run over cells to compute the distances
            for( int m2 = 0; m2 < M; m2++ ){
                
                // Calculate the distances between two cells with minimum image convention
                double dist_X = cells[m2].xpos[i] - cells[m1].xpos[i];
                dist_X = dist_X - Lx*dnearbyint( dist_X/Lx );
                
                double dist_Y = cells[m2].ypos[i] - cells[m1].ypos[i];
                dist_Y = dist_Y - Ly*dnearbyint( dist_Y/Ly );
                
                double dist_R = sqrt( dist_X*dist_X + dist_Y*dist_Y );
                
                // Get the bin corresponding to the distance
                int bin = inearbyint( dist_R );
                
                // Calculate the mean subtracted velocity fields
                double vx_2 = cells[m2].xvel[i] - vx_mean[i];   // x velocities
                double vx_1 = cells[m1].xvel[i] - vx_mean[i];
                double vy_2 = cells[m2].yvel[i] - vy_mean[i];   // y velocities
                double vy_1 = cells[m1].yvel[i] - vy_mean[i];
                
                #ifdef VEL_DIR
                // ++++++++
                // Directional velocity correlation
                
                double v_sqrt = sqrt(vx_1*vx_1 + vy_1*vy_1);
                double v_par = (vx_1*vx_2 + vy_1*vy_2)/v_sqrt;
                double v_perp = (-vy_1*vx_2 + vx_1*vy_2)/v_sqrt;
      
                #endif
                
                // Calculate the velocity correlation
                double v_temp = vx_2*vx_1 + vy_2*vy_1;
                
                cvv[bin] += v_temp;
                cvv_n[bin] += 1;
                
                #ifdef VEL_DIR
                cvv_par[bin] += v_par;
                cvv_perp[bin] += v_perp;
                #endif
                
                // Keep track of maximum distance between the cells
                if( dist_R > max_dist ){
                    max_dist = dist_R;
                }
                
            } // for second cell m2
        } // for first cell m1
        
    } // for running over time
    
    // Get the last bin
    int total_bins = (int)max_dist;    // Total number of bins

    ofstream sp_out("sp_vacf.txt");
    // Write the spatial velocity correlation to a file
    for( int i = 0; i < total_bins; i++ ){
        sp_out << i/(2*R_avg) << "\t" << cvv[i]/(cvv_n[i]*denom_avg) << "\n";
    }
    
    #ifdef VEL_DIR
    ofstream sp_dir_out("sp_dir_vacf.txt");
    for( int i = 0; i < total_bins; i++ ){
        sp_dir_out << i/(2*R_avg) << "\t" << cvv_par[i]/(cvv_n[i]) << "\t" << cvv_perp[i]/(cvv_n[i]) << "\n";
    }
    #endif
    
    
    #endif // end of SP_VACF
    
    
	#ifdef VEL_SEG
    
    t_offset = 0;
    t_offset_2 = T;
    t_offset_total = t_offset_2 - t_offset;

    int max_corr_length_avg = 0;
    
    // Calculate the speeds
    for( int ii = t_offset; ii < t_offset_2; ii++ ){
     
        Graph graph;
        graph.setNumNodes( (int)M/5 );
        int ust_limit = graph.getNumNodes();
        graph.setThreshold( Rcut, 2 );
        double th = graph.getThreshold();
        double threshold2 = th*th;
        
        vector<double> v;
        for(int m = 0; m < M; m++){
            v.push_back( sqrt( cells[m].xvel[ii]*cells[m].xvel[ii] + cells[m].yvel[ii]*cells[m].yvel[ii] )  );
        }
        
        // Get the fastest 20% of cells and store the cell indexes
//        ofstream debug_txt("debug.txt");
        int break_cnt = 0;
        set<int> idx;
        for( auto i: sort_indexes(v) ){
//            debug_txt << i << "\t" << v[i] << endl;
            idx.insert( i );
            break_cnt++;
            if( break_cnt == ust_limit ){
                break;
            }
        }
        
//        ofstream debug2_txt("debug2.txt");
//        for( int m = 0; m < M; m++ ){
//            debug2_txt << m << "\t" << v[m] << endl;
//        }

        // Build the graph between the fastest 20% of cells
        // ( Note that because of the "set" data structure 'idx',
        // .. list of cells is ordered. )
//        ofstream fastest_out("fastest_cells.txt");
        for( auto i: idx ) {
            
//            fastest_out << i << endl;
            set<int> destination;
            
            // Calculate the distances between cells with brute force
            for( auto j: idx ){
                
                if( i != j ){  // so that the graph is undirected ..
                    
                    double dx = cells[j].xpos[ii] - cells[i].xpos[ii];
                    dx -= Lx*dnearbyint( dx/Lx );
                    double dy = cells[j].ypos[ii] - cells[i].ypos[ii];
                    dy -= Ly*dnearbyint( dy/Ly );
                    double d2 = dx*dx + dy*dy;
                    if( d2 < threshold2 ) { destination.insert( j ); }
                    
                } // if
                
            } // loop over j
            
            graph.addEdge( i, destination );
        
        } // loop over i
        
        // Do a depth-first search to get the connected components
        vector<int> max_corr_length;
        vector<bool> visited(M);
        for( auto i: idx ) { visited[i] = false; }
        stack<int> nodes;
        
        // For all the nodes that are the fastest ..
        for( auto i: idx ){
            
            // If the node isn't visited ..
            if( !visited[i] ){
                max_corr_length.push_back( dfsFromNode( graph, i, nodes, visited ) );
//                cout << endl;
            }
        }
        
        max_corr_length_avg += *max_element( max_corr_length.begin(), max_corr_length.end() );
    }
    
    max_corr_length_avg /= t_offset_total;
    
    ofstream maxcorr_out;
    maxcorr_out.open( "max_corr_length.txt", ios::out | ios::app );
    maxcorr_out << dens << "\t" << Fmc << "\t" << eps << "\t" << max_corr_length_avg << endl;
    maxcorr_out.close();

	#endif // end of VEL_SEG

    #endif // end of VEL
    
} // end of the program
