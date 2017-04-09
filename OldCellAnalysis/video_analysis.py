#!/usr/bin/python

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import subprocess
import sh
import glob
import os
import re
from drawnow import drawnow

# Simulation information
sample_data = 1
nsamp = 10000*sample_data
dt = 1e-3*nsamp
M = 400
N = 40 
N = M*N

fl = open("../Data/phi0.95eps10fm3/Data/Data1/cell.xyz", "r")
fr = 0

all_cells_file = open("img/dat/all_cells_0.txt")
xcm1 = []
ycm1 = []
for l1 in all_cells_file:
    A1 = l1.split()
    xcm1.append( float(A1[3]) )
    ycm1.append( float(A1[4]) )
   

X = []
Y = []
Yc = []
trajx = []
trajy = []
trajx_abs = []
trajy_abs = []

for line in fl:
    
    A = line.split()
    
    if len(A) > 0:
    
        if A[0] != "C":
            
            fig = plt.figure()
            ax1 = plt.subplot2grid((2,3), (0,0))
            #ax1.margins(1)
            ax2 = plt.subplot2grid((2,3), (1,0), sharex=ax1)
            #ax2.margins(1)
            ax3 = plt.subplot2grid((2,3), (1,1), sharey=ax2)
            #ax3.margins(1)
            ax4 = plt.subplot2grid((2,3), (0,1), sharex=ax3, sharey=ax1)
            #ax4.margins(1)
            ax5 = plt.subplot2grid((2,3), (1,2), sharex=ax3)
            #ax5.margins(1)
            ax6 = plt.subplot2grid((2,3), (0,2), sharex=ax4, sharey=ax5)
            #ax6.margins(1)

            fr = fr+1 
    
            all_cells_file = open("img/dat/all_cells_" + str(fr-1) + ".txt")
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
            
            fastest_cells_file = open("img/dat/fastest_cells_" + str(fr-1) + ".txt")
            cell_idx_fastest = []
            for l1 in fastest_cells_file:
                A1 = l1.split()
                cell_idx_fastest.append( int(A1[1]) )
                
            #cell_idx_fastest = np.asarray(cell_idx_fastest,dtype=int)
            
            msd_file = open("img/dat/msd_field_" + str(fr-1) + ".txt")
            msd = []
            msd_idx = []
            for l1 in msd_file:
                A1 = l1.split()
                msd_idx.append( int(A1[0]) )
                msd.append( float(A1[1]) )
                
            msd = np.asarray(msd, dtype=float)
                
            bond_file = open("img/dat/bond_field_" + str(fr-1) + ".txt")
            bond = []
            bond_idx = []
            for l1 in bond_file:
                A1 = l1.split()
                bond_idx.append( int(A1[0]) )
                bond.append( float(A1[1]) )
                
            bond = np.asarray(bond, dtype=float)
            
            # Monomer positions (ax1)
            fig.subplots_adjust(hspace=0.1)
            ax1.scatter(X,Y, s=0.3, c=Yc, cmap=plt.cm.brg, edgecolors='None', alpha=1)
            ax1.set_ylabel('y/l')
            ax1.set_title('Cell positions')
            ax1.set_xlim((-5,250))
            ax1.set_ylim((-5,250))
            ax1.set_aspect('equal')
            plt.setp(ax1.get_xticklabels(),visible=False)
            ax1.yaxis.set_ticks(np.arange(0,300,75))
            ax1.tick_params(axis='both', which='major', labelsize=8)
            
            # Center of mass trajectories (sax2)
            ax2.scatter(xcm_abs[70:75], ycm_abs[70:75], c='b', s=0.8, alpha=0.8)
            ax2.scatter(xcm_abs[270:275], ycm_abs[270:275], c='b', s=0.8, alpha=0.8)
            xx = np.asarray(trajx_abs)
            yy = np.asarray(trajy_abs)
            ax2.plot(xx[:,70:75], yy[:,70:75], color='b', alpha=0.3)
            ax2.plot(xx[:,270:275], yy[:,270:275], color='b', alpha=0.3)
            ax2.set_ylabel('y/l')
            ax2.set_xlabel('x/l')
            ax2.set_title('C.M. trajectories')
            ax2.set_xlim((-10,240))
            ax2.set_ylim((-10,240))
            ax2.set_aspect('equal')
            ax2.yaxis.set_ticks(np.arange(0,230,50))
            ax2.xaxis.set_ticks(np.arange(0,230,50))
            ax2.tick_params(axis='both',which='major',labelsize=8)


            # Velocity field (ax3)
            ax3.quiver(xcm, ycm, vx, vy, alpha=1, color='r', headaxislength=5)
            ax3.plot(xcm, ycm, 'o', ms=1.3, alpha=0.3)
            for idx in cell_idx_fastest:
                ax3.plot(xcm[idx],ycm[idx], 'o', ms=1.3, alpha=0.7, color='b')
            ax3.set_xlabel('x/l')
            ax3.set_title('Velocity field')
            ax3.set_xlim((-5,250))
            ax3.set_ylim((-5,250))
            ax3.set_aspect('equal')
            plt.setp(ax3.get_yticklabels(),visible=False)
            ax3.tick_params(axis='both', which='major', labelsize=8)
            ax3.xaxis.set_ticks(np.arange(0,300,75))
            
  
            # Center of mass positions (ax4)       
            ax4.scatter(xcm, ycm, s=2, c=ycm1, cmap=plt.cm.brg, edgecolors='None', alpha=1)
            ax4.set_title('C.M. positions')
            ax4.set_xlim((-5,250))
            ax4.set_ylim((-5,250))
            ax4.set_aspect('equal')
            ax4.yaxis.set_ticks(np.arange(0,300,75))
            ax4.xaxis.set_ticks(np.arange(0,300,75))
            plt.setp(ax4.get_xticklabels(),visible=False)
            plt.setp(ax4.get_yticklabels(),visible=False)
            plt.figtext(0.5, 1, 't = ' + str(fr*10))
            ax4.tick_params(axis='both',which='major',labelsize=8)
            #ax4.annotate('Test', xy=(1, 0), xycoords='axes fraction', fontsize=16,
                #xytext=(0, 0), textcoords='offset points',
                #ha='right', va='top')
            
            # Local bond correlation (ax5)
            area_bond = np.pi * (3 * bond)**2 
            line5 = ax5.scatter(xcm, ycm, c = bond, s=area_bond, cmap=plt.cm.brg, alpha=1)
            ax5.set_xlabel('x/l')
            ax5.set_title('$\psi_{4}$')
            ax5.set_xlim((-5,250))
            ax5.set_ylim((-5,250))
            ax5.set_aspect('equal')
            plt.setp(ax5.get_yticklabels(),visible=False)
            ax5.tick_params(axis='both', which='major', labelsize=8)
            ax5.xaxis.set_ticks(np.arange(0,300,75))
            plt.subplots_adjust(bottom=0.1, top=0.9, left=0.125, right=0.9)
            cax = plt.axes([0.95, 0.1, 0.02, 0.35])
            plt.colorbar(line5,cax=cax)
            cax.tick_params(width=1.1,labelsize=6)

            # Mean square displacement (ax6)    
            area_msd = np.pi * (1 * msd)**2 
            line6 = ax6.scatter(xcm, ycm, c = msd, s=area_msd, cmap=plt.cm.brg, alpha=0.7)
            ax6.set_title('$MSD/4R^2$')
            ax6.set_xlim((-5,250))
            ax6.set_ylim((-5,250))
            ax6.set_aspect('equal')
            ax6.yaxis.set_ticks(np.arange(0,300,75))
            ax6.xaxis.set_ticks(np.arange(0,300,75))
            plt.setp(ax6.get_xticklabels(),visible=False)
            plt.setp(ax6.get_yticklabels(),visible=False)
            ax6.tick_params(axis='both',which='major',labelsize=8)
        
            plt.subplots_adjust(bottom=0.1, top=0.9, left=0.125, right=0.9, hspace=0.3, wspace=0.2)
            cax = plt.axes([0.95, 0.55, 0.02, 0.35])
            plt.colorbar(line6,cax=cax)
            cax.tick_params(width=1.1,labelsize=6)
                    
            plt.savefig('./img/'+'frame-'+'{0:05d}'.format(fr)+'.png',dpi=200,bbox_inches='tight',pad_inches=0.08)
            plt.clf()
            X = []
            Y = []
    
            cnt = 0
        
        elif A[0] == "C":
            
            cnt = cnt+1
            
            if cnt <= N:
                
                X.append( float(A[1]) )
                Y.append( float(A[2]) )
                
                if fr <= 1:
                    Yc.append( float(A[2]) )
                    
            cnt = 0
                
                