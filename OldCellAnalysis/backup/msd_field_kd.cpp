
//+++++++++++++++
//++ Data analysis with postprocessing for field calculation
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
#include "pointbox.h"
#include "kdtree.h"

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

// +++++++



int main( int argc, char *argv[] ){
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++ Get the data
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    // Choose the data to be used
    const int nsamp = samplePostData;       // Sampling interval
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
    SimBox box;
    box.setWidth( boxSize_x );
    box.setHeight( boxSize_y );
    
    const double Lx = box.getWidth();
    const double Ly = box.getHeight();
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++ Do the analysis
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    int t_offset = 0;
    int t_offset_2 = T;
    int t_offset_total = t_offset_2 - t_offset;
    
    double cell_diameter_2 = 4*R_avg*R_avg;
    
    int n = 4;      // Number of neighbors
    
    for( int ii = t_offset; ii < t_offset_2; ii++ ){
        
        // Create the points that will make up the tree
        const int MM = 9*M;
        vector<Point<2> > cell_vec(MM);
        
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //++ MSD calculation
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        string msd_out_path = "img/dat/msd_field_" + to_string(ii) + ".txt";
        ofstream msd_out( msd_out_path );
        for( int m = 0; m < M; m++ ){
            
            // MSD part
            double delta_x = cells[m].xpos[ii] - cells[m].xpos[0];
            double delta_y = cells[m].ypos[ii] - cells[m].ypos[0];
            msd_out << m << "\t" << (delta_x*delta_x + delta_y*delta_y)/cell_diameter_2 << "\n";
            
            // Image points
            cells[m].getImag(box, ii);
            
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
        
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //++ Local bond orientation
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        // Build the tree from the points
        KDtree<2> kdtree(cell_vec);
        
        vector<double> theta_real(M, 0.0);
        vector<double> theta_img(M, 0.0);
        for(int m = 0; m < M; m++){                 // Run over cells
            
            int *nn = new int[n];
            double *dn = new double[n];
            kdtree.nnearest(m,nn,dn,n);
            
            for(int j = 0; j < n; j++){             // Run over neighbors
                
                double dy = dist_xy(cell_vec[nn[j]], cell_vec[m], 1);
                double dx = dist_xy(cell_vec[nn[j]], cell_vec[m], 0);
                double theta = 4*atan2( dy, dx );
                theta_real[m] += cos( theta );
                theta_img[m] += sin( theta );
                
            }       // Run over neighbors
            
        }           // Run over cells

        string bond_out_path = "img/dat/bond_field_" + to_string(ii) + ".txt";
        ofstream bond_out( bond_out_path );
        for(int m = 0; m < M; m++){     // Run over cells
            
            bond_out << m << "\t" << sqrt( theta_real[m]*theta_real[m] + theta_img[m]*theta_img[m] )/n << endl;
            
        }   // End of run over cell
        
    } // End of run over time
    
} // End of the program
