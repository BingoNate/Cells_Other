
//+++++++++++++++++
//++ Data analysis with postprocessing for density correlators
//+++++++++++++++++

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

#define PAIR_CORR
#define STATIC_STRUCT
#define INT_SCATTER
#define DYNAMIC_STRUCT

// +++++++
// DATA TYPES
// +++++++

// +++++
// Cell
// +++++

class SimBox;

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
    yi = ypos[i] - a.height*ifloor(ypos[i]/a.height);
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
    const int nsamp = samplePostData;                   // Sampling interval
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
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++ Density correlators
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    // +++++++++
    // Pair correlation function
    // +++++++++
    
    #ifdef PAIR_CORR
    
    int t_offset = 0;
    int t_offset_2 = T;
    int t_offset_total = t_offset_2 - t_offset;
    
    double max_dist = 0.0;
    int sug_max_dist = (int)fmin( boxSize_x, boxSize_y );  // Suggest a sane approximation for maximum distance
    vector<double> n_cnt(sug_max_dist, 0.0);
    double rho = M/(Lx*Ly);
    double n_ideal = 0;
    
    for( int i = t_offset; i < t_offset_2; i++ ){                               // Run over all the timeseries
        
        for( int m1 = 0; m1 < M-1; m1++ ){                  // Run over cells to compute the distances in the first place
            for( int m2 = m1+1; m2 < M; m2++ ){
                
                // Calculate the distances between two cells
                double dist_X = cells[m2].xpos[i] - cells[m1].xpos[i];
                dist_X = dist_X - Lx*dnearbyint( dist_X/Lx );
                
                double dist_Y = cells[m2].ypos[i] - cells[m1].ypos[i];
                dist_Y = dist_Y - Ly*dnearbyint( dist_Y/Ly );
                
                double dist_R = sqrt( dist_X*dist_X + dist_Y*dist_Y );
                
                // Get the bin corresponding to the distance
                int bin = inearbyint( dist_R );
                
                // Count how many times you are making the calculation of correlation function
                n_cnt[bin] += 2;
                
            } // for m2
        } // for m1
        
    } // for time
    
    // Get the last bin
    int total_bins = (int)sug_max_dist;    // Total number of bins
    
    ofstream pair_out("pair_corr.txt");
    // Calculate the pair correlation function
    for( int i = 0; i < total_bins; i++ ){
        n_cnt[i] /= ( t_offset_total*M );
        n_ideal  = 2*pi*i*rho;
        pair_out << i/(2*R_avg) << "\t" << n_cnt[i]/n_ideal << "\n";
    }
    
    #endif // end of pair correlation function calculation

    
    
    // +++++++++++++++++++++++++++
    
    
    // +++++++++
    // Static structure factor
    // +++++++++
    
    #ifdef STATIC_STRUCT
    
    // Variable declarations
    int L = (int)fmin(Lx,Ly)/2;              // Effective box size
    double k;                              // Wavevector
    double tut;
    double cterm_x = 0;
    double sterm_x = 0;
    double cterm_y = 0;
    double sterm_y = 0;
    double first_avg_x = 0;
    double cos_avg_x = 0;
    double sin_avg_x = 0;
    double first_avg_y = 0;
    double cos_avg_y = 0;
    double sin_avg_y = 0;
    double Sx;
    double Sy;
    
    ofstream struct_out("static_struct.txt");

    // Calculate the static structure factor
    for(int i = 1; i < L; i++){  // Run over wavevectors
        
        k = 2*pi*i/L;

        for(int ii = t_offset; ii < t_offset_2; ii++){   // Run over time series
            
            for(int m = 0; m < M; m++){   // Run over cells
                
                
                // x wavevector -- Wavevector * x[m][ii] where m is the cell index, ii is the time index
                double xp = cells[m].xpos[ii];
                xp -= Lx*ifloor(xp/Lx);
                tut = k * xp;
                
                cterm_x += cos( tut );
                sterm_x += sin( tut );
    
                
                // y wavevector -- Wavevector * y[m][ii] where m is the cell index, ii is the time index
                double yp = cells[m].ypos[ii];
                yp -= Ly*ifloor(yp/Ly);
                tut = k * yp;
                
                cterm_y += cos( tut );
                sterm_y += sin( tut );
            
                
            } // Run over cells
            
            // x
            first_avg_x += cterm_x*cterm_x + sterm_x*sterm_x;
            
            cterm_x = 0;
            sterm_x = 0;
            
            // y
            first_avg_y += cterm_y*cterm_y + sterm_y*sterm_y;
            
            cterm_y = 0;
            sterm_y = 0;
            
        } // Run over time series
        
        // Take the average
        Sx = first_avg_x/t_offset_total;
        Sy = first_avg_y/t_offset_total;
        
        struct_out << k*R_avg << "\t" << (Sx+Sy)/(2*M) << endl;
        
        first_avg_x = 0;
        first_avg_y = 0;
        
    } // Run over wavevectors
    
    #endif // Static structure factor end
    
    
    //++++++++++++
    // Self intermediate scattering function
    //+++++++++++
    
    #ifdef INT_SCATTER
    
    // Calculate the normalization factor
    vector<double> Fr(D, 0.0), Fi(D, 0.0);

    ofstream inter_out;
    inter_out.open("inter_scattering.txt");
    inter_out << 0 << "\t" << 1.0 << "\n";
    inter_out.close();
    
    k = 0.5;
    
    cterm_x = 0;
    cterm_y = 0;
    sterm_x = 0;
    sterm_y = 0;
    
    double complex_norm_sq_x = 0;
    double complex_norm_sq_y = 0;
    double real_sq_x = 0;
    double real_sq_y = 0;
    double imag_sq_x = 0;
    double imag_sq_y = 0;

    // Calculate the intermediate scattering function
    for(int d = 1; d < D; d++){ // Run over delay time
        
        for(int ii = 0; ii < TD; ii++){   // Run over different time origins
            
            for(int m = 0; m < M; m++){   // Run over cells
                
                double dx = cells[m].xpos[ii+d]-cells[m].xpos[ii];
                cterm_x += cos( k*dx );
                sterm_x += sin( k*dx );
                
                double dy = cells[m].ypos[ii+d]-cells[m].ypos[ii];
                cterm_y += cos( k*dy );
                sterm_y += sin( k*dy );

                
            } // Run over cells
            
            real_sq_x += cterm_x;
            real_sq_y += cterm_y;
            imag_sq_x += sterm_x;
            imag_sq_y += sterm_y;
            complex_norm_sq_x += cterm_x*cterm_x + sterm_x*sterm_x;
            complex_norm_sq_y += cterm_y*cterm_y + sterm_y*sterm_y;

            cterm_x = 0;
            cterm_y = 0;
            sterm_x = 0;
            sterm_y = 0;
            
        } // Run over different time origins
        
        complex_norm_sq_x /= TD;
        complex_norm_sq_y /= TD;
        real_sq_x /= TD;
        real_sq_y /= TD;
        imag_sq_x /= TD;
        imag_sq_y /= TD;
        
        //Fr[d] = (complex_norm_sq_x+complex_norm_sq_y)/(2*M*M);
        Fr[d] = (real_sq_x+real_sq_y)/(2*M);
        Fi[d] = (imag_sq_x+imag_sq_y)/(2*M);
        
        inter_out.open("inter_scattering.txt", ios::out | ios::app);
        inter_out << d*dtSamp << "\t" << Fr[d] << "\n";
        inter_out.close();
        
        complex_norm_sq_x = 0;
        complex_norm_sq_y = 0;
        real_sq_x = 0;
        real_sq_y = 0;
        imag_sq_x = 0;
        imag_sq_y = 0;
        
    }   // Run over delay time
    
    #endif // Self intermediate scattering function end
    
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    
    //++++++++++++
    // Dynamic structure factor
    //+++++++++++
    
    #ifdef DYNAMIC_STRUCT
    
    k = 0.5;
    
    // Calculate the dynamic structure factor by taking the discrete Fourier transform in time domain
    
    ofstream dyns_out;
    dyns_out.open("dynamic_struct.txt");
    dyns_out.close();
    int maxK = (int)D/2;
    for(int j = 0; j < D; j++){ // Time series
        
        double sumreal = 0;
        double sumimag = 0;
        
        j -= maxK;
        
        double aci = 2*pi*j/D;

        for(int i = 0; i < D; i++){ // Frequency series
            
            sumreal += Fr[i]*cos(aci*i) + Fi[i]*sin(aci*i);
            sumimag += Fi[i]*cos(aci*i) - Fr[i]*sin(aci*i);
            
        } // Time series
        
        // Normalise
        sumreal /= D;
        sumimag /= D;
        
        dyns_out.open("dynamic_struct.txt", ios::out | ios::app);
        dyns_out << j*2*pi/D << "\t" << (sumreal*sumreal + sumimag*sumimag)/(2*pi) << "\n";
        dyns_out.close();
        
        j += maxK;
        
    } // Frequency series
    
    #endif // Dynamic structure factor end
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
} // end of the program
