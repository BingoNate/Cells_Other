
//+++++++++++++++
//++ Data analysis with postprocessing for block averaging
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
    const int nsamp = samplePostData;        // Sampling interval
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
    // ++ Block averaging for average displacement
    // +++++++++++++++
    
    ofstream block_out("block_displacement.txt");
    
    // Determine the block lengths
    vector<int> block_length;
    //int max_exponent = 13;
    //for( int i = 1; i < max_exponent; i++ ){
    //    block_length.push_back( pow(2,i) );
    //}
    int max_exponent = 5000;
    int l = 3;
    while( l < max_exponent ) { block_length.push_back( l ); l += 5; }
    
    // For each block length calculate the displacement average and store it for that particular block length
    for( auto i: block_length ){
        
        double M2D_ens = 0;
        double M2D = 0;
        vector<double> displacement_avg_per_block;
        
        // For the block length i, access block indices ..
        int n = 0;
        int max_n = i;
        int M2D_cnt = 0;
        int d = 1;  // Time delay
        double total_avg = 0;
        double total_avg_cnt = 0;
        
        NEXT_BLOCK:while( n < max_n-1 ){
            
            // .. and calculate displacement average for time delay d
            for( int m = 0; m < M; m++ ){         // Run over cells
                
                double x_origin = cells[m].xpos[n];
                double x_delay = cells[m].xpos[n+d];
                double y_origin = cells[m].ypos[n];
                double y_delay = cells[m].ypos[n+d];
                
                // Difference
                double delta_x = x_delay - x_origin;
                double delta_y = y_delay - y_origin;
                
                // Square of difference
                double delta_x_2 = delta_x * delta_x;
                double delta_y_2 = delta_y * delta_y;
                
                M2D_ens += sqrt( delta_x_2 + delta_y_2 );
                
            } // for cell indice
            
            // Take the ensemble average
            M2D += M2D_ens/M;
            M2D_cnt++;
            
            // Reset holders
            M2D_ens = 0;
            
            n++;
            
        } // while n -- inside the block
        
        // Store the displacement average and reset the holders
        displacement_avg_per_block.push_back( M2D/M2D_cnt );
        total_avg += M2D/M2D_cnt;
        total_avg_cnt++;
        M2D = 0;
        M2D_cnt = 0;
        
        // Go to the next block for this block length or calculate the standard deviation for this the block length
        max_n += i;
        if( max_n < T ){
            
            goto NEXT_BLOCK;
            
        }
        else{
            
            total_avg /= total_avg_cnt;
            
            int Tn = (int)T/i;   // Total number of blocks for block length i
            
            // Calculate average displacement over all blocks
            double displacement_avg_all_blocks = 0;
            for(int j = 0; j < Tn; j++) { displacement_avg_all_blocks += displacement_avg_per_block[j]; }
            displacement_avg_all_blocks /= (double)Tn;
            
            // Calculate the standard deviation for the block length i
            double vrnc = 0;
            double std_dev = 0;
            for(int j = 0; j < Tn; j++){
                //vrnc += (displacement_avg_per_block[j] - total_avg)*(displacement_avg_per_block[j] - total_avg);
                std_dev += (displacement_avg_per_block[j] - total_avg)*(displacement_avg_per_block[j] - total_avg);
            }
            //block_out << i << "\t" << vrnc/Tn << endl;
            // Write the block standard error for this block length i
            block_out << i << "\t" << sqrt( std_dev/Tn )/sqrt(Tn) << endl;
            total_avg = 0;
            total_avg_cnt = 0;
            
        }
        
    } // for auto i -- different block lengths
    
    
    // ++++++++++++++
    
}
