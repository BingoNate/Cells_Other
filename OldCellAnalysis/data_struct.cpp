
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
#include "variables.h"
#include "basic.h"
#include "pointbox.h"
#include "kdtree.h"

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
    int totalNumPart;
    void getImag(SimBox a, int i);
    void setOldPos(int i);
    
    Particle(int n_size, int t_size) : totalNumPart(n_size), xpos(t_size, 0.0), ypos(t_size, 0.0) {}
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
    int numCell;
    int numPart;
    vector<double> xpos;
    vector<double> ypos;
    vector<double> xi;
    vector<double> yi;
    double xcom;
    double ycom;
    double xcomi;
    double ycomi;
    double radii;
    void getPart(vector<Particle> & p, int pn, int tn);
    void getImag(SimBox a);
    void calcNumPart();
    void calcCOM(SimBox a);
    
    Cell(int m_size, int n_size) : numCell(m_size), xpos(n_size, 0.0), ypos(n_size, 0.0), xi(n_size, 0.0), yi(n_size, 0.0) { }
    ~Cell() { xpos.clear(); ypos.clear(); xi.clear(); yi.clear(); }
    
};

void Cell::calcNumPart(){
    numPart = (int)ceil(2*pi*radii);    // NB = 2piR (careful if B = 1 !!)
}

void Cell::getPart(vector<Particle> & p, int pn, int tn){
    for(int n = 0; n < numPart; n++){
        int i = pn+n;
        xpos[n] = p[i].xpos[tn];
        ypos[n] = p[i].ypos[tn];
    }
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

void Cell::getImag( SimBox a ){
    for(int n = 0; n < numPart; n++){
        xi[n] = xpos[n] - a.width*ifloor(xpos[n]/a.width);
        yi[n] = ypos[n] - a.height*ifloor(ypos[n]/a.height);
    }
}

void Cell::calcCOM( SimBox a ){
    xcom = 0; ycom = 0;
    for(int n = 0; n < numPart; n++){
        xcom += xpos[n];
        ycom += ypos[n];
    }
    xcom /= numPart;
    ycom /= numPart;
    xcomi = xcom - a.width*ifloor(xcom/a.width);
    ycomi = ycom - a.height*ifloor(ycom/a.height);
}

void CellBox::setBoxProp( SimBox a ){
    Rx = (double)ceil( a.width/ifloor(a.width/separator) );
    Ry = (double)ceil( a.height/ifloor(a.height/separator) );
    Bx = ifloor( a.width/Rx );
    By = ifloor( a.height/Ry );
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


// Function for reading the positions
void readPosData( vector<Particle> & particles, char * file_name, int & sample_cnt, int sample_data ){
    
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
                    
                    if( i != jp && i/N != jp/N ){
                        double dx = particles[i].xpos[ti] - particles[jp].xpos[ti];
                        dx -= Lx*dnearbyint( dx/Lx );
                        double dy = particles[i].ypos[ti] - particles[jp].ypos[ti];
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
