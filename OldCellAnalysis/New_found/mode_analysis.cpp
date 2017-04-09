
//+++++++++++++++
//++ Data analysis with postprocessing for spatial modes
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
#include <Eigen/Dense>
#include "variables.h"
#include "basic.h"

using namespace Eigen;
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
    
    // Decide whether to read the x value or the y value
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

// ++++++++


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
    
    // Read the data to the cells vector
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
    //++ Density of modes
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    // Time interval
    int t_offset = 0;
    int t_offset_2 = T;
    int t_offset_total = t_offset_2 - t_offset;
    
    // Mean positions
    vector<double> x_mean( M, 0.0 );
    vector<double> y_mean( M, 0.0 );
    
    // Covariance matrix
    const int ML = 2*M;
    vector<double> var_pos(M);
    MatrixXd Dab(ML,ML);
    for( int m1 = 0; m1 < ML; m1++ ){
        for(int m2 = 0; m2 < ML; m2++ ){
            Dab(m1,m2) = 0.0;
        }
        var_pos.push_back( 0.0 );
    }
    
    // Calculate the mean positions in time
    for( int i = t_offset; i < t_offset_2; i++ ){
        
        for( int m = 0; m < M; m++ ){
            x_mean[m] += cells[m].xpos[i];
            y_mean[m] += cells[m].ypos[i];
        }
        
    }
    for( int m = 0; m < M; m++ ){
        x_mean[m] /= t_offset_total;
        y_mean[m] /= t_offset_total;
    }

    // Time loop
    for( int i = t_offset; i < t_offset_2; i++ ){
        
        
        // Calculate the mean subtracted positions in x and y in each time frame
        for( int m = 0; m < M; m++ ){
            cells[m].meanSubtractedPosition( i, x_mean[m], y_mean[m] );
            
            // Calculate the variance in fluctuations around the mean position for each cell
            var_pos[m] += cells[m].xmean*cells[m].xmean + cells[m].ymean*cells[m].ymean;
        }
        
        // Construct the displacement covariance matrix
        for( int j1 = 0; j1 < ML; j1 += 2 ){
            for( int j2 = 0; j2 < ML; j2 += 2 ){
                
                int m1 = j1/2;
                int m2 = j2/2;
                
                // <ux_a,ux_b> at (j1,j2) - (0,0)
                Dab(j1,j2) += cells[m1].xmean * cells[m2].xmean;
                
                // <ux_a,uy_b> at (j1,j2+1) - (0,1)
                Dab(j1,j2+1) += cells[m1].xmean * cells[m2].ymean;
                
                // <uy_a,ux_b> at (j1+1,j2) - (1,0)
                Dab(j1+1,j2) += cells[m1].ymean * cells[m2].xmean;
                
                // <uy_a,uy_b> at (j1+1,j2+1) - (1,1)
                Dab(j1+1,j2+1) += cells[m1].ymean * cells[m2].ymean;

                
            } // Second cell
        } // First cell
        
        
    } // End of time loop
    
    // Variance file
    ofstream var_out("var_pos.txt");
    for( int m = 0; m < M; m++ ){
        var_out << var_pos[m]/t_offset_total << "\n";
    }
    
    // Normalize the covariance matrix
    for( int m1 = 0; m1 < ML; m1++ ){
        for(int m2 = 0; m2 < ML; m2++ ){
            Dab(m1,m2) /= t_offset_total;
        }
    }
    
    // Calculate the eigenvalues and eigenvectors of the displacement covariance matrix
    EigenSolver<MatrixXd> Dab_eig;
    Dab_eig.compute( Dab, true );
    
    
    // Calculate the density of modes
    string path = argv[3];
    string modes_out_path = "modes" + path + ".txt";
    ofstream modes_out( modes_out_path );
    string eigen_out_path = "eigen" + path + ".txt";
    ofstream eigen_out( eigen_out_path );
    for( int i = 0; i < ML; i++ ){
        
        double a = Dab_eig.eigenvalues()(i).real();
        modes_out << sqrt(1/a) << "\n";
        
        if( i == 0 ){
            int m = -1;
            for( int j = 0; j < ML; j+=2 ){
                m++;
                double a1 = Dab_eig.eigenvectors().col(i)[j].real();
                double a2 = Dab_eig.eigenvectors().col(i)[j+1].real();
                eigen_out << a << "\t" << cells[m].xpos[0] << "\t" << cells[m].ypos[0] << "\t" << a1 << "\t" << a2 << "\n";
            }
        }
        
        //vector<complex<double> > eigvecs(ML);
        //eigvecs = Dab_eig.eigenvectors().col(i);
        
//        for( int j = 0; j < M; j += 2 ){
//            
//            //complex<double> a1 = eigvecs[j].real();      // x value of the eigenvector
//            //complex<double> a2 = eigvecs[j+1].real();    // y value of the eigenvector
//            double a1 = Dab_eig.eigenvectors().col(i)[j].real();
//            double a2 = Dab_eig.eigenvectors().col(i)[j+1].real();
//            
//            eigen_out << a1 << "\t" << a2 << "\n";
//            
//        }

    }
    
    
    //cout << "Eigenvalues\t" << Dab_eig.eigenvalues()[0] << "\n";
    //cout << "Eigenvectors\t" << Dab_eig.eigenvectors().col(0)(1,0) << "\n";
    //complex<double> a = Dab_eig.eigenvalues()(0);
    //double b = abs(a);
    //cout << b << "\n";
    //cout << "Eigenvectors\t" << Dab_eig.eigenvectors().col(0) << "\t" << Dab_eig.eigenvectors().col(1) << "\n";
    //cout << "Eigenvalues\n" << Dab_eig.eigenvalues() << "\n\n";
    //cout << "Eigenvectors\n" << Dab_eig.eigenvectors().col(380) << "\n";

    
    
} // End of the program
