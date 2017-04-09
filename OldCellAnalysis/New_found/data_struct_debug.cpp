
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <set>
#include <stack>
#include <map>
#include <algorithm>
#include "variables.h"
#include "basic.h"
#include "pointbox.h"
#include "kdtree.h"
#include "delaunay.h"

using namespace std;

#define pi M_PI

// +++++++
// DATA TYPES
// +++++++

// Forward declarations
class SimBox;
class Cell;

// +++++
// Particle
// +++++

class Particle{
    
public:
    friend class Cell;
    vector<double> xpos;
    vector<double> ypos;
    double xold;
    double yold;
    double xi;
    double yi;
    int cell_idx;
    int totalSamp;
    
    void getImag(SimBox a, int i);
    void setOldPos(int i);
    
    Particle(int t_size, int c_idx) : xpos(t_size, 0), ypos(t_size, 0), cell_idx(c_idx), totalSamp(t_size) {}
    ~Particle() { xpos.clear(); ypos.clear(); }
};

void Particle::setOldPos(int i){
    xold = xpos[i];
    yold = ypos[i];
}

// +++++


// +++++
// Cell
// +++++

class Cell
{
    
public:
    int c_idx;
    int numCell;
    int numPart;
    int totalSamp;
    vector<double> xpos;
    vector<double> ypos;
    vector<double> xvel;
    vector<double> yvel;
    double xi;
    double yi;
    double radii;
    vector<int> p_idx;
    set<int> neighs;
    
    void calcNumPart();
    void getImag(SimBox a, int i);
    void setImag(SimBox a, vector<Point<2> > & cell_vec);
    
    Cell(int c_idx, int t_size, int n_size, int m_size) : totalSamp(t_size), c_idx(c_idx), xpos(t_size, 0), ypos(t_size, 0), xvel(t_size-1, 0), yvel(t_size-1, 0), numCell(m_size) { }
    ~Cell() { xpos.clear(); ypos.clear(); xvel.clear(); yvel.clear(); }
    
};

void Cell::calcNumPart(){
    numPart = (int)ceil(2*pi*radii/B);    // NB = 2piR (careful if B = 1 or not !!)
}

// +++++

// +++++
// Cell box with link lists
// +++++

class CellBox{
    
public:
    double separator;
    double Rx;
    double Ry;
    int Bx;
    int By;
    int BN;
    vector<int> head_box;
    vector<int> link_box;
    CellBox(double rsep, int npar) : separator(rsep), link_box(npar, -1) { }
    ~CellBox() {}
    
    void setBoxProp(SimBox a);
};

// +++++


// +++++
// The periodic box
// +++++

class SimBox
{
    
public:
    friend class CellBox;
    friend class Cell;
    friend class Particle;
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

void Particle::getImag(SimBox a, int i){
    xi = xpos[i] - a.width*ifloor(xpos[i]/a.width);
    yi = ypos[i] - a.height*ifloor(ypos[i]/a.height);
}

void Cell::getImag( SimBox a, int i ){
    xi = xpos[i] - a.width*ifloor(xpos[i]/a.width);
    yi = ypos[i] - a.height*ifloor(ypos[i]/a.height);
}

void Cell::setImag(SimBox a, vector<Point<2> > & cell_vec){
    
    cell_vec[c_idx].x[0] = xi;
    cell_vec[c_idx].x[1] = yi;
    
    cell_vec[numCell+c_idx].x[0] = xi - a.width;
    cell_vec[numCell+c_idx].x[1] = yi + a.height;
    
    cell_vec[2*numCell+c_idx].x[0] = xi;
    cell_vec[2*numCell+c_idx].x[1] = yi + a.height;
    
    cell_vec[3*numCell+c_idx].x[0] = xi + a.width;
    cell_vec[3*numCell+c_idx].x[1] = yi + a.height;
    
    cell_vec[4*numCell+c_idx].x[0] = xi + a.width;
    cell_vec[4*numCell+c_idx].x[1] = yi;
    
    cell_vec[5*numCell+c_idx].x[0] = xi + a.width;
    cell_vec[5*numCell+c_idx].x[1] = yi - a.height;
    
    cell_vec[6*numCell+c_idx].x[0] = xi;
    cell_vec[6*numCell+c_idx].x[1] = yi - a.height;
    
    cell_vec[7*numCell+c_idx].x[0] = xi - a.width;
    cell_vec[7*numCell+c_idx].x[1] = yi - a.height;
    
    cell_vec[8*numCell+c_idx].x[0] = xi - a.width;
    cell_vec[8*numCell+c_idx].x[1] = yi;

}


void CellBox::setBoxProp( SimBox a ){
    Rx = (double)ceil( a.width/ifloor(a.width/separator) );
    Ry = (double)ceil( a.height/ifloor(a.height/separator) );
    Bx = ceil( a.width/Rx );
    By = ceil( a.height/Ry );
    BN = Bx*By;
    for(int i = 0; i < BN; i++){
        head_box.push_back( -1 );
    }
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
    int getNumStep() { return numStep; }
    int getNumDelay() { return numDelay; }
    int getNumTotalDelay() { return numTotalDelay; }
    
private:
    int numStep;
    int numDelay;
    int numTotalDelay;
    
};

void Data::setNumStep( double nst ){
    numStep = (int)nst;
}

void Data::setNumDelay( double nd ){
    numDelay = (int)nd;
}

void Data::setNumTotalDelay(){
    numTotalDelay = (int)(numStep - numDelay);
}

// +++++


// +++++
// Graph
// +++++

// ( Note that edges are lines connecting nodes or vertices )

class Graph
{
    
public:
    Graph(int nds, double ths);
    ~Graph() { edges.clear(); }
    void addEdge( int nd_idx, set<int> dest );
    map<int,set<int> > edges; // [node, destination]
    int numNodes;
    double threshold;
    
};

Graph::Graph(int nds, double ths) : numNodes(nds), threshold(ths) {}


void Graph::addEdge( int nd_idx, set<int> dest ){
    edges.insert( pair<int,set<int> >( nd_idx, dest ) );
}


// +++++

// Function to sort a vector by value while keeping the indexes
bool sortRule(int i1, int i2, vector<double> &v){ return v[i1] > v[i2];  };

struct IdxCompare{
    const vector<double> & target;
    IdxCompare(const vector<double> & target): target(target) {}
    bool operator ()(int i1, int i2) const { return target[i1] > target[i2]; }
};

template <typename T>
vector<int> sort_indexes( const vector<T> &v ) {
    
    // Initialize original index locations
    vector<int> idx( v.size() );
    for( int i = 0; i != idx.size(); ++i ) idx[i] = i;
    
    // Sort indexes by comparing values in v
    //sort( idx.begin(), idx.end(), [&v](int i1, int i2) { return v[i1] > v[i2]; } );
    sort( idx.begin(), idx.end(), IdxCompare(v) );
    
    return idx;
}



// Function for depth-first search from a starting node
int dfsFromNode( Graph& graph, int i, stack<int>& S, vector<bool>& visited, ofstream & connected_comp_list_out ){
    
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
            connected_comp_list_out << i << "\t";
            connected_comp++;           // Add a connected component when you visit a new node
            set<int> neighbors;
            neighbors = graph.edges[i];
            for(set<int>::iterator it = neighbors.begin(); it != neighbors.end(); it++){
                i = *it;
                S.push( i );
            }
        } // if the node is visited before, get out
        
    } // while loop to check if the stack is empty or not
    
    return connected_comp;
    
}

// Function for reading the radii information
void readRadii(){
    
    ifstream r_in;
    r_in.open( "radii.dat", ios::binary | ios::in );
    r_in.read( (char *)R, M*sizeof(double) );
    r_in.close();
    L = 0;
    for(int m = 0; m < M; m++){
        A[m] = pi*R[m]*R[m];
        A0[m] = A0_*A[m];
        int Li = (int)ceil(2*pi*R[m]);
        Npart[m] = Li;
        L += Li;
    }
    
}

// Function for reading the positions
int readPosData( vector<Particle> & particles, char * file_name, int & sample_cnt, int sample_data ){
    
    // Inputs
    int part_cnt = -1;
    int t_cnt = 0;
    ifstream inputData;
    inputData.open( file_name );
    
    // Holder for each line
    string line_holder;
    
    // Yield error in case of a problem
    if( inputData.fail() ){
        cerr << "The file cannot be opened!\n";
        exit(1);
    }
    
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
                char char_holder;
                double double_holder1, double_holder2, double_holder3;
                
                // Now read the values in a line, finally ...
                while ( stream_holder >> char_holder >> double_holder1 >> double_holder2 >> double_holder3 )
                {
                    part_cnt++;
                    if( part_cnt == L ){
                        t_cnt++;
                        part_cnt = 0;
                    }
                    
                    // Each element in the line consists of particle coordinates at that given time instant
                    particles[part_cnt].xpos[t_cnt] = double_holder1;
                    particles[part_cnt].ypos[t_cnt] = double_holder2;
                    
                } // Read a particular line
                
            } // read line by line
            
        } // data sampling
    } // while data is good
    
    // Good, you're done!!
    inputData.close();
    
    return t_cnt;
    
}

// +++++++

// Calculate virial stress
double calcVirStress( double dx, double dy, double eps48, double sig2 ){
    
    double d2 = dx*dx + dy*dy;
    double Fx = 0.;
    double Fy = 0.;
    if( d2 < Rcut2 ){
        
        double p3 = sig2/d2;
        p3 = p3*p3*p3;
        double tempFx = eps48*( p3*p3 - p3/2. )*dx/d2;
        double tempFy = eps48*( p3*p3 - p3/2. )*dy/d2;
        Fx = tempFx;
        Fy = tempFy;
    }
    
    return dx*Fx + dy*Fy;
    
}


// +++++++

// BOX LINK LIST BASED VERLET LIST

void verletList( vector<Particle> & particles, CellBox & cellbox, const double Lx, const double Ly, vector<int> & numVerList, vector<int> & verList, vector<int> & verListOffset, int verListSize, SimBox box, int ti );

// Check if an update to the Verlet list is necessary
void updateList( vector<Particle> & particles, CellBox & cellbox, const double Lx, const double Ly, vector<int> & numVerList, vector<int> & verList, vector<int> & verListOffset, int verListSize, SimBox box, int ti ){
    
    double oldMax = 0;
    double newMax = 0;
    
    // Check if updates is required
    for(int i = 0; i < L; i++){
        double dx = particles[i].xpos[ti] - particles[i].xold;
        double dy = particles[i].ypos[ti] - particles[i].yold;
        double dr = sqrt( dx*dx + dy*dy );
        
        if( dr > newMax ){
            oldMax = newMax;
            newMax = dr;
        }
    }
    
    // Update the Verlet list if necessary
    if( (newMax + oldMax) > skin ){
        verletList( particles, cellbox, Lx, Ly, numVerList, verList, verListOffset, verListSize, box, ti );
    }
    
}

// Find the box number with periodic boundary conditions
int boxNum_periodic( double pos, const double siz, const int num ){
    int ara = ifloor( pos/siz );
    if( ara < 0 )
        ara += num;
    else if( ara == num )
        ara = 0;
    
    return ara;
}

// Initialize the box link lists
void initList( vector<Particle> & particles, CellBox & cellbox, SimBox box, int ti ){
    
    for(int i = 0; i < cellbox.BN; i++){
        cellbox.head_box[i] = -1;
    }
    
    for(int i = 0; i < L; i++){
        particles[i].getImag(box,ti);
        int kx = boxNum_periodic( particles[i].xi, cellbox.Rx, cellbox.Bx );
        int ky = boxNum_periodic( particles[i].yi, cellbox.Ry, cellbox.By );
        int kbox = kx*cellbox.By + ky;
        cellbox.link_box[i] = cellbox.head_box[kbox];
        cellbox.head_box[kbox] = i;
    }
    
}

// Create the box link list based Verlet list
void verletList( vector<Particle> & particles, CellBox & cellbox, const double Lx, const double Ly, vector<int> & numVerList, vector<int> & verList, vector<int> & verListOffset, int verListSize, SimBox box, int ti ){
    
    for(int i = 0; i < L; i++){
        numVerList[i] = 0;
        particles[i].setOldPos( ti );
    }
    
    int vlisti = 0;
    initList( particles, cellbox, box, ti );
    
    for(int i = 0; i < L; i++){
        verListOffset[i] = vlisti;
        particles[i].getImag(box, ti);
        int kx = boxNum_periodic( particles[i].xi, cellbox.Rx, cellbox.Bx );
        int ky = boxNum_periodic( particles[i].yi, cellbox.Ry, cellbox.By );
        
        for(int ncelxi = kx-1; ncelxi < kx+2; ncelxi++){
            for(int ncelyi = ky-1; ncelyi < ky+2; ncelyi++){
                
                int nbox_x = ncelxi - cellbox.Bx*ifloor( ncelxi/cellbox.Bx );
                if( nbox_x < 0 ){
                    nbox_x += cellbox.Bx;
                }
        
                int nbox_y = ncelyi - cellbox.By*ifloor( ncelyi/cellbox.By );
                if( nbox_y < 0 ){
                    nbox_y += cellbox.By;
                }
                
                int kbox = nbox_x*cellbox.By + nbox_y;
                int jp = cellbox.head_box[kbox];
                
                while( jp != -1 ){
                    
                    if( i != jp && particles[i].cell_idx != particles[jp].cell_idx ){
                        double dx = particles[jp].xpos[ti] - particles[i].xpos[ti];
                        dx -= Lx*dnearbyint( dx/Lx );
                        double dy = particles[jp].ypos[ti] - particles[i].ypos[ti];
                        dy -= Ly*dnearbyint( dy/Ly );
                        double d2 = dx*dx + dy*dy;
                        
                        if( d2 < Rver2 ){
                            numVerList[i]++;
                            verList[vlisti++] = jp;
                        }
                    }
                    
                    jp = cellbox.link_box[jp];
                    
                }
            }
        }
        
    }
    
    
}

// ++++++
