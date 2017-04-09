#!/usr/bin/python

import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("nsamp", help="Sampling frequency", type=float)
parser.add_argument("M", help="Number of cells", type=int)
parser.add_argument("cell_file", help="xyz file")
args = parser.parse_args()

nsamp = args.nsamp
dt = 1e-3
dtSamp = dt*nsamp

M = args.M

downlim = -5
uplim = 500

fl = open(args.cell_file, "r")
fr = -1

all_cells_file = open("img/dat/cells_0.txt")
xcm1 = []
ycm1 = []
for l1 in all_cells_file:
    A1 = l1.split()
    xcm1.append( float(A1[3]) )
    ycm1.append( float(A1[4]) )

L_cnt = 0
cnt = 0
for line in fl:
    A = line.split()
    L = int(A[0])
    break 

X = []
Y = []
Yc = []
trajx = []
trajy = []
trajx_abs = []
trajy_abs = []

# Box information for virial stress 
box_info_file = open("img/dat/box_info.txt", "r")
for l1 in box_info_file:
    A1 = l1.split()
    Bx = float(A1[0])
    By = float(A1[1])
    rsepx = float(A1[2])
    rsepy = float(A1[3])
    
bx, by = np.mgrid[slice(0, Bx*rsepx, rsepx),slice(0, By*rsepy, rsepy)]

for line in fl:
    
    A = line.split()
    
    if len(A) > 1:
    
        if A[0] == "C":
            
            X.append( float(A[1]) )
            Y.append( float(A[2]) )
            
            if fr < 0:
                Yc.append( float(A[2]) )
                
            L_cnt = L_cnt+1
    
            if L_cnt == L:
            
                # Set up the plot
                fig = plt.figure()
                ax1 = plt.subplot2grid((2,3), (0,0))
                ax2 = plt.subplot2grid((2,3), (1,0), sharex=ax1)
                ax3 = plt.subplot2grid((2,3), (1,1), sharey=ax2)
                ax4 = plt.subplot2grid((2,3), (0,1), sharex=ax3, sharey=ax1)
                ax5 = plt.subplot2grid((2,3), (1,2), sharex=ax3)
                ax6 = plt.subplot2grid((2,3), (0,2), sharex=ax4, sharey=ax5)

                # Next frame    
                fr = fr+1 
        
                # Set up center of mass positions of cells and their velocities
                all_cells_file = open("img/dat/cells_" + str(fr) + ".txt")
                v = []
                vx = []
                vy = []
                xcm = []
                ycm = []
                cell_idx = []
                time_idx = []
                xcm_abs = []
                ycm_abs = []
                for l1 in all_cells_file:
                    A1 = l1.split()
                    time_idx.append( int(A1[0]) )
                    cell_idx.append( int(A1[1]) )
                    v.append( float(A1[2]) )
                    xcm.append( float(A1[3]) )
                    ycm.append( float(A1[4]) )
                    vx.append( float(A1[5]) )
                    vy.append( float(A1[6]) )
                    xcm_abs.append( float(A1[7]) )
                    ycm_abs.append( float(A1[8]) )
                
                trajx.append(xcm)
                trajy.append(ycm)
                trajx_abs.append(xcm_abs)
                trajy_abs.append(ycm_abs)
                
                # Set up fastest cells
                fastest_cells_file = open("img/dat/fastest_cells_" + str(fr) + ".txt")
                cell_idx_fastest = []
                for l1 in fastest_cells_file:
                    A1 = l1.split()
                    cell_idx_fastest.append( int(A1[1]) )
                  
                # Set up neighbors
                neigh_file = open("img/dat/neighs_" + str(fr) + ".txt")
                cell_idx_neigh = []
                num_neighs = []
                for l1 in neigh_file:
                    A1 = l1.split()
                    cell_idx_neigh.append( int(A1[0]) )
                    num_neighs.append( float(A1[1]) )
                    
                num_neighs = np.asarray(num_neighs, dtype=float)
              
                  
                # Set up pressure
                pres_file = open("img/dat/local_pressure_" + str(fr) + ".txt")
                pres = []
                pres_idx = []
                for l1 in pres_file:
                    A1 = l1.split()
                    pres_idx.append( int(A1[0]) )
                    pres.append( float(A1[1]) )
                    
                pres = np.asarray(pres, dtype=float)
                
                
                # PLOTS
                
                # Monomer positions (ax1)
                fig.subplots_adjust(hspace=0.1)
                ax1.scatter(X,Y, s=0.3, c=Yc, cmap=plt.cm.brg, edgecolors='None', alpha=1)
                ax1.set_ylabel('y/l')
                ax1.set_title('Cell positions')
                ax1.set_xlim((downlim,uplim))
                ax1.set_ylim((downlim,uplim))
                ax1.set_aspect('equal')
                plt.setp(ax1.get_xticklabels(),visible=False)
                ax1.yaxis.set_ticks(np.arange(0,uplim,50))
                ax1.tick_params(axis='both', which='major', labelsize=8)
                
                
                # Center of mass trajectories (ax2)
                ax2.scatter(xcm_abs[70:75], ycm_abs[70:75], c='b', s=0.8, alpha=0.8)
                ax2.scatter(xcm_abs[270:275], ycm_abs[270:275], c='b', s=0.8, alpha=0.8)
                xx = np.asarray(trajx_abs)
                yy = np.asarray(trajy_abs)
                ax2.plot(xx[:,70:75], yy[:,70:75], color='b', alpha=0.3)
                ax2.plot(xx[:,270:275], yy[:,270:275], color='b', alpha=0.3)
                ax2.set_ylabel('y/l')
                ax2.set_xlabel('x/l')
                ax2.set_title('C.M. trajectories')
                ax2.set_xlim((downlim,uplim))
                ax2.set_ylim((downlim,uplim))
                ax2.set_aspect('equal')
                ax2.yaxis.set_ticks(np.arange(0,uplim,50))
                ax2.xaxis.set_ticks(np.arange(0,uplim,50))
                ax2.tick_params(axis='both',which='major',labelsize=8)
    
    
                # Velocity field (ax3)
                ax3.quiver(xcm, ycm, vx, vy, alpha=1, color='r', headaxislength=5)
                ax3.plot(xcm, ycm, 'o', ms=1.3, alpha=0.3)
                for idx in cell_idx_fastest:
                    ax3.plot(xcm[idx],ycm[idx], 'o', ms=1.3, alpha=0.7, color='b')
                ax3.set_xlabel('x/l')
                ax3.set_title('Velocity field')
                ax3.set_xlim((downlim,uplim))
                ax3.set_ylim((downlim,uplim))
                ax3.set_aspect('equal')
                plt.setp(ax3.get_yticklabels(),visible=False)
                ax3.tick_params(axis='both', which='major', labelsize=8)
                ax3.xaxis.set_ticks(np.arange(0,uplim,50))
                
      
                # Center of mass positions (ax4)       
                ax4.scatter(xcm, ycm, s=2, c=ycm1, cmap=plt.cm.brg, edgecolors='None', alpha=1)
                ax4.set_title('C.M. positions')
                ax4.set_xlim((downlim,uplim))
                ax4.set_ylim((downlim,uplim))
                ax4.set_aspect('equal')
                ax4.yaxis.set_ticks(np.arange(0,uplim,50))
                ax4.xaxis.set_ticks(np.arange(0,uplim,50))
                plt.setp(ax4.get_xticklabels(),visible=False)
                plt.setp(ax4.get_yticklabels(),visible=False)
                plt.figtext(0.5, 1, 't = ' + str(fr*10))
                ax4.tick_params(axis='both',which='major',labelsize=8)
                
                
                # Local pressure map (ax5)
                area_pres = np.pi * (1.0e-10 * pres)**2 
                line5 = ax5.scatter(xcm, ycm, c=pres, s=area_pres, cmap=plt.cm.brg, alpha=1)
                ax5.set_xlabel('x/l')
                ax5.set_title('P')
                ax5.set_xlim((downlim,uplim))
                ax5.set_ylim((downlim,uplim))
                ax5.set_aspect('equal')
                plt.setp(ax5.get_yticklabels(),visible=False)
                ax5.tick_params(axis='both', which='major', labelsize=8)
                ax5.xaxis.set_ticks(np.arange(0,uplim,50))
                plt.subplots_adjust(bottom=0.1, top=0.9, left=0.125, right=0.9)
                cax = plt.axes([0.95, 0.1, 0.02, 0.35])
                plt.colorbar(line5,cax=cax)
                cax.tick_params(width=1.1,labelsize=6)
    
    
                # Number of neighbors per cell  (ax6)    
                area_neigh = np.pi * (num_neighs)**2 
                line6 = ax6.scatter(xcm, ycm, c=num_neighs, s=num_neighs, cmap=plt.cm.brg, alpha=1)
                ax6.set_title('Interparticle Virial')
                ax6.set_xlim((downlim,uplim))
                ax6.set_ylim((downlim,uplim))
                ax6.set_aspect('equal')
                ax6.yaxis.set_ticks(np.arange(0,uplim,50))
                ax6.xaxis.set_ticks(np.arange(0,uplim,50))
                plt.setp(ax6.get_xticklabels(),visible=False)
                plt.setp(ax6.get_yticklabels(),visible=False)
                ax6.tick_params(axis='both',which='major',labelsize=8)
                cax = plt.axes([0.95, 0.55, 0.02, 0.35])
                plt.colorbar(line6,cax=cax)
                cax.tick_params(width=1.1,labelsize=6)
            
            
                plt.subplots_adjust(bottom=0.1, top=0.9, left=0.125, right=0.9, hspace=0.3, wspace=0.2)
                cax = plt.axes([0.95, 0.55, 0.02, 0.35])
                plt.colorbar(line6,cax=cax)
                cax.tick_params(width=1.1,labelsize=6)
                    
                    
                plt.savefig('./img/'+'frame-'+'{0:05d}'.format(fr)+'.png',dpi=200,bbox_inches='tight',pad_inches=0.08)
                plt.clf()
                X = []
                Y = []
        
                L_cnt = 0
                
