
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
    double xmean;
    double ymean;
    int numCell;
    void meanSubtractedPosition( int i, double x_mean_i, double y_mean_i );
    
    Cell() {}
    ~Cell() {}
    
};

void Cell::meanSubtractedPosition( int i, double x_mean_i, double y_mean_i ){
    
    xmean = xpos[i] - x_mean_i;
    ymean = ypos[i] - y_mean_i;
    
}

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
    Box box;
    box.setWidth( boxSize_x );
    box.setHeight( boxSize_y );
    
    const double Lx = box.getWidth();
    const double Ly = box.getHeight();
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++ Do the analysis
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    
    // +++++++++++++++
    // ++ MSD calculation
    // +++++++++++++++
    
    // Declarations
    double M2D_time = 0;
    double M2D_ens = 0;
    double M4D_time = 0;
    double M4D_ens = 0;
    double M2D = 0;
    double M4D = 0;
    
    // Mean positions
    vector<double> x_mean( T, 0.0 );
    vector<double> y_mean( T, 0.0 );
    for(int i = 0; i < T; i++){
        for(int m = 0; m < M; m++){
            x_mean[i] += cells[m].xpos[i];
            y_mean[i] += cells[m].ypos[i];
        }
        x_mean[i] /= M;
        y_mean[i] /= M;
    }
    
    ofstream msdp_out("msd_post.txt");
    ofstream msd_out("msd_field.txt");
    ofstream m4d_out("nongauss.txt");
    
    // Calculate MSD
    for(int d = 0; d < D; d++){                 // Run over delay times
        for(int m = 0; m < M; m++){             // Run over cells
            for(int i = 0; i < TD; i++){        // Run over different time origins
               
                double x_origin = cells[m].xpos[i] - x_mean[i];
                double x_delay = cells[m].xpos[i+d] - x_mean[i+d];
                double y_origin = cells[m].ypos[i] - y_mean[i];
                double y_delay = cells[m].ypos[i+d] - y_mean[i+d];
                
                // Difference
                double delta_x = x_delay-x_origin;
                double delta_y = y_delay-y_origin;
                
                // Square of difference
                double delta_x_2 = delta_x*delta_x;
                double delta_y_2 = delta_y*delta_y;
                
                // 4th power of difference
                double delta_x_4 = delta_x_2*delta_x_2;
                double delta_y_4 = delta_y_2*delta_y_2;
                
                M2D_time += delta_x_2 + delta_y_2;
                M4D_time += delta_x_4 + delta_y_4 + 2*delta_x_2*delta_y_2;
                
            }
            
            
            // Take a time average
            M2D_ens += M2D_time/TD;
            M4D_ens += M4D_time/TD;
            
            // Write individual MSDs
            msd_out << M2D_time/TD << "\t"; // Average over time origins
            
            // Reset holders
            M2D_time = 0;
            M4D_time = 0;

        }
        
        // Take an ensemble average
        M2D = M2D_ens/M;
        M4D = M4D_ens/M;
        
        // Reset holders
        M2D_ens = 0;
        M4D_ens = 0;
        
        // Write non-Gaussian parameter
        m4d_out << d*dtSamp << "\t" << ( M4D/(2*M2D*M2D) ) - 1 << "\n";
        msdp_out << d*dtSamp << "\t" << M2D << endl; // Average over cells
        M4D_ens = 0;
        
        // Write an empty line in individual MSDs
        msd_out << "\n";
        
        
    }

    
    // ++++++++++++++
    
}
