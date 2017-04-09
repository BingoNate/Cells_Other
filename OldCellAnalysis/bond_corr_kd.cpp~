

//+++++++++++++++++
//++ Data analysis with postprocessing for bond orientation
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
#include "pointbox.h"
#include "kdtree.h"

using namespace std;

#define pi M_PI

#define BOND_ORIENT
#define BOND_CORR

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
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++ Bond orientation
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    
    #ifdef BOND_ORIENT
    
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
    
    ofstream local_bond_out("local_bond.txt");
    vector<double> theta_real(M);
    vector<double> theta_img(M);
    vector<double> psi6(M);
    for(int m = 0; m < M; m++){
        theta_real[m] = 0;
        theta_img[m] = 0;
        psi6[m] = 0;
    }
    double glob_psi6 = 0;
    double glob_psi6_avg = 0;
    double psi6_bar_temp = 0;
    vector<double> psi6_bar;
    
    int n = 4;              // Number of neighbors
    
    int t_offset = 0;
    int t_offset_2 = T;
    int t_offset_total = t_offset_2 - t_offset;
    
    for( int i = t_offset; i < t_offset_2; i++){        // Run over time
        
        // Create the points that will make up the tree
        const int MM = 9*M;
        vector<Point<2> > cell_vec(MM);
        
        // Consider all the image points as well
        for(int m = 0; m < M; m++){
            
            cells[m].getImag(box, i);
            
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
        
        // Build the tree from the points
        KDtree<2> kdtree(cell_vec);
        
        for(int m = 0; m < M; m++){                 // Run over cells
            
            int *nn = new int[n];
            double *dn = new double[n];
            kdtree.nnearest(m,nn,dn,n);
            
            for(int j = 0; j < n; j++){             // Run over neighbors
                double dy = dist_xy(cell_vec[nn[j]], cell_vec[m], 1);
                double dx = dist_xy(cell_vec[nn[j]], cell_vec[m], 0);
                double theta = 4*atan2( dy, dx );
                double cos_theta = cos( theta );
                double sin_theta = sin( theta );
                theta_real[m] += cos_theta;
                theta_img[m] += sin_theta;
            }
            
        }
        
        for(int m = 0; m < M; m++){
            double tmp_value = sqrt( theta_real[m]*theta_real[m] + theta_img[m]*theta_img[m] )/n;
            //double tmp_value = theta_real[m]/n;
            psi6[m] += tmp_value;
            glob_psi6 += tmp_value;
            psi6_bar_temp += tmp_value;
            theta_real[m] = 0;
            theta_img[m] = 0;
        }
        
        glob_psi6_avg += glob_psi6/M;
        glob_psi6 = 0;
        psi6_bar.push_back( psi6_bar_temp/M );
        psi6_bar_temp = 0;
        
    } // for i (time)
    
    for(int m = 0; m < M; m++){
        psi6[m] /= t_offset_total;
        local_bond_out << psi6[m] << "\n";
    }
    
    ofstream glob_bond_out;
    glob_bond_out.open( "glob_bond.txt", ios::out | ios::app );
    glob_bond_out << dens << "\t" << Fmc << "\t" << eps << "\t" << glob_psi6_avg/t_offset_total << "\n";
    glob_bond_out.close();
    
    double pre_sq = 0;
    double pre_avg = 0;
    for( int i = t_offset; i < t_offset_2; i++ ){
        pre_sq += psi6_bar[i]*psi6_bar[i];
        pre_avg += psi6_bar[i];
    }
    pre_sq /= t_offset_total;
    pre_avg /= t_offset_total;
    ofstream susc_out;
    susc_out.open( "susceptibility.txt", ios::out | ios::app );
    susc_out << dens << "\t" << Fmc << "\t" << eps << "\t" << M*(pre_sq - pre_avg*pre_avg) << "\n";
    susc_out.close();
    
    #endif // End of bond orientation
    
    #ifdef BOND_CORR
    
    double max_dist = 0.0;
    int sug_max_dist = (int)fmax( boxSize_x, boxSize_y );  // Suggest a sane approximation for maximum distance
    vector<double> cvv(sug_max_dist, 0.0);
    vector<double> cvv_n(sug_max_dist, 0.0);

    for( int i = t_offset; i < t_offset_2; i++ ){           // Run over all the timeseries
        
        for( int m1 = 0; m1 < M-1; m1++ ){                  // Run over cells to compute the distances
            for( int m2 = m1+1; m2 < M; m2++ ){
                
                // Calculate the distances between two cells with minimum image convention
                double dist_X = cells[m2].xpos[i] - cells[m1].xpos[i];
                dist_X = dist_X - Lx*dnearbyint( dist_X/Lx );
                
                double dist_Y = cells[m2].ypos[i] - cells[m1].ypos[i];
                dist_Y = dist_Y - Ly*dnearbyint( dist_Y/Ly );
                
                double dist_R = sqrt( dist_X*dist_X + dist_Y*dist_Y );
                
                // Get the bin corresponding to the distance
                int bin = inearbyint( dist_R );
                
                // Calculate the bond correlation function
                cvv[bin] += psi6[m2]*psi6[m1];
                cvv_n[bin] += 1;
                
                // Keep track of maximum distance between the cells
                if( dist_R > max_dist ){
                    max_dist = dist_R;
                }
                
            } // for m2
        } // for m1
        
    } // for running over time
    
    // Get the last bin
    int total_bins = (int)max_dist;    // Total number of bins
    
    ofstream bond_corr_out("bond_corr.txt");
    // Write the bond correlation function
    for( int i = 0; i < total_bins; i++ ){
        bond_corr_out << i/(2*R_avg) << "\t" << cvv[i]/cvv_n[i] << "\n";
    }
    
    #endif // End of bond correlation
    
} // End of the program
