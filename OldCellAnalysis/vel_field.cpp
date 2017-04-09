
//+++++++++++++++
//++ Data analysis with postprocessing for velocity field calculation
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
int dfsFromNode( Graph& graph, int i, stack<int>& S, vector<bool>& visited, ofstream& connected_comp_out ){
    
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
            connected_comp_out << i << "\t";
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
    int sample_data = 5;
    const int nsamp = samplePostData*sample_data;       // Sampling interval
    const double dtSamp = nsamp*dt;                     // Time units between two data points
    
    // Instantiate the data structure
    vector<Cell> cells( M );
    
    // Set filenames
    char * x_input_file = argv[1];        // Filename for the x data
    char * y_input_file = argv[2];        // Filename for the y data
    
    // Read the data to the cells
    int sample_cnt = -1;
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
    //++ Velocity field calculation
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
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
    
    int t_offset = 0;
    int t_offset_2 = T;
    int t_offset_total = t_offset_2 - t_offset;
    
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
        string fastest_cells_out_path = "img/dat/fastest_cells_" + to_string(ii) + ".txt";
        ofstream fastest_cells_out( fastest_cells_out_path );
        int break_cnt = 0;
        set<int> idx;
        for( auto i: sort_indexes(v) ){
            fastest_cells_out << ii*dtSamp << "\t" << i << "\t" << v[i] << "\n";
            idx.insert( i );
            break_cnt++;
            if( break_cnt == ust_limit ){
                break;
            }
        }
        
        string all_cells_out_path = "img/dat/all_cells_" + to_string(ii) + ".txt";
        ofstream all_cells_out( all_cells_out_path );
        for( int m = 0; m < M; m++ ){
            all_cells_out << ii*dtSamp << "\t" << m << "\t" << v[m] << "\t" << (cells[m].xpos[ii]-Lx*ifloor(cells[m].xpos[ii]/Lx)) << "\t" << (cells[m].ypos[ii]-Ly*ifloor(cells[m].ypos[ii]/Ly)) << "\t" << cells[m].xvel[ii] << "\t" << cells[m].yvel[ii] << "\t" << cells[m].xpos[ii] << "\t" << cells[m].ypos[ii] << "\n";
        }

        // Build the graph between the fastest 20% of cells
        // ( Note that because of the "set" data structure 'idx',
        // .. list of cells is ordered. )
        for( auto i: idx ) {
            
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
        vector<bool> visited(M);
        for( auto i: idx ) { visited[i] = false; }
        stack<int> nodes;
        
        vector<int> max_corr_length;
        int max_corr_length_hold;
        string max_corr_length_out_path = "img/dat/max_corr_lengths_" + to_string(ii) + ".txt";
        ofstream max_corr_length_out( max_corr_length_out_path );
        string connected_comp_out_path = "img/dat/connected_components_" + to_string(ii) + ".txt";
        ofstream connected_comp_out( connected_comp_out_path );
        
        // For all the nodes that are the fastest ..
        for( auto i: idx ){
            
            // If the node isn't visited ..
            if( !visited[i] ){
                max_corr_length.push_back( dfsFromNode( graph, i, nodes, visited, connected_comp_out ) );
                connected_comp_out << "\n";
            }
            max_corr_length_hold = *max_element( max_corr_length.begin(), max_corr_length.end() );
        }
        
        max_corr_length_out << ii*dtSamp << "\t" << max_corr_length_hold << "\n";
        
    } // End of run over time
    
} // End of the program
