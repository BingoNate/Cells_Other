

//+++++++++++++++++
//++ Data analysis with postprocessing for finding neighbors
//+++++++++++++++++


#include <set>
#include <map>
#include <stack>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <map>
#include <algorithm>
#include <boost/multi_array.hpp>
#include "variables.h"
#include "basic.h"
#include "kdtree.h"
#include "pointbox.h"

using namespace std;

#define pi M_PI


// +++++++
// DATA TYPES
// +++++++

// +++++
// Cell
// +++++

class SimBox;       // Forward declaration

class Cell
{
    
public:
    vector<double> xpos;
    vector<double> ypos;
    double xi;
    double yi;
    vector<double> xvel;
    vector<double> yvel;
    int numCell;
    void getImag(SimBox a, int i);
    
    Cell() {}
    ~Cell() {}
    
};

// +++++


// +++++
// The periodic box
// +++++

class SimBox
{
    
public:
    friend class Cell;
    SimBox() {}
    ~SimBox() {}
    void setWidth( double wdth );
    void setHeight( double hght );
    double getWidth() { return width; }
    double getHeight() { return height; }
    double getArea() { return width*height; }
    
private:
    double width;
    double height;
    
};

void SimBox::setWidth( double wdth ){
    width = wdth;
}

void SimBox::setHeight( double hght ){
    height = hght;
}

void Cell::getImag(SimBox a, int i){
    xi = xpos[i] - a.width*ifloor(xpos[i]/a.width);
    yi = ypos[i] - a.height*ifloor(ypos[i]/a.width);
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
    threshold = 2*R + carpan*cutoff;
}


// +++++


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
    const int nsamp = samplePostData;                   // Sampling interval
    const double dtSamp = nsamp*dt;                     // Time units between two data points
    
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
    SimBox box;
    box.setWidth( boxSize_x );
    box.setHeight( boxSize_y );
    
    const double Lx = box.getWidth();
    const double Ly = box.getHeight();
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++ Do the analysis
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    const int MM = 9*M;
    vector<Point<2> > cell_vec(MM);
    
    for(int m = 0; m < M; m++){
        
        cells[m].getImag(box, 0);
        
        cell_vec[m].x[0] = cells[m].xi;
        cell_vec[m].x[1] = cells[m].yi;
        
        cell_vec[M+m].x[0] = cells[m].xi - Lx;
        cell_vec[M+m].x[1] = cells[m].yi + Ly;
        
        cell_vec[2*M+m].x[0] = cells[m].xi;
        cell_vec[2*M+m].x[1] = cells[m].yi + Ly;
        
        cell_vec[3*M+m].x[0] = cells[m].xi + Lx;
        cell_vec[3*M+m].x[1] = cells[m].yi + Ly;

        cell_vec[4*M+m].x[0] = cells[m].xi + Lx;
        cell_vec[4*M+m].x[1] = cells[m].yi;
        
        cell_vec[5*M+m].x[0] = cells[m].xi + Lx;
        cell_vec[5*M+m].x[1] = cells[m].yi - Ly;

        cell_vec[6*M+m].x[0] = cells[m].xi;
        cell_vec[6*M+m].x[1] = cells[m].yi - Ly;

        cell_vec[7*M+m].x[0] = cells[m].xi - Lx;
        cell_vec[7*M+m].x[1] = cells[m].yi - Ly;
        
        cell_vec[8*M+m].x[0] = cells[m].xi - Lx;
        cell_vec[8*M+m].x[1] = cells[m].yi;
        
    }
    
    KDtree<2> kdtree(cell_vec);
    
//    int nmax = 100;
//    int *list = new int[nmax];
//    int a = kdtree.locatenear( cell_vec[0], 20, list, nmax );
//    cout << "The point is at " << cell_vec[0].x[0] << "\t" << cell_vec[0].x[1] << endl;
//    for(int i = 0; i < a; i++){
//        cout << list[i] << "\t" << cell_vec[list[i]].x[0] << "\t" << cell_vec[list[i]].x[1] << "\t";
//        cout << kdtree.disti( 0, list[i] ) << endl;
//    }
   
    int jpt = 0;
    int n = 6;
    int *nn = new int[n];
    double *dn = new double[n];
    kdtree.nnearest(jpt,nn,dn,n);
    cout << "The point is at " << cell_vec[jpt].x[0] << "\t" << cell_vec[jpt].x[1] << endl;
    for( int i = 0; i < n; i++ ){
        cout << nn[i] << "\t" << dn[i] << endl;
    }
    
//    Point<2> a;
//    a.x[0] = 33;
//    a.x[1] = 23;
//    int c = kdtree.nearest(a);
//    cout << c << endl;
//    cout << kdtree.locate(c) << endl;
//    cout << cell_vec[c].x[0] << "\t" << cell_vec[c].x[1] << endl;
    
}

