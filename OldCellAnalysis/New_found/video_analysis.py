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
uplim = 510

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

# Box information for virial stress 
box_info_file = open("img/dat/box_info.txt", "r")
for l1 in box_info_file:
    A1 = l1.split()
    Bx = float(A1[0])
    By = float(A1[1])
    rsepx = float(A1[2])
    rsepy = float(A1[3])
    
bx, by = np.mgrid[slice(0, Bx*rsepx, rsepx),slice(0, By*rsepy, rsepy)]

X = []
Y = []
Yc = []
trajx = []
trajy = []
trajx_abs = []
trajy_abs = []

tick_interval = 100

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
                    
                # Set up MSD
                msd_file = open("img/dat/msd_field_" + str(fr) + ".txt")
                msd = []
                msd_idx = []
                for l1 in msd_file:
                    A1 = l1.split()
                    msd_idx.append( int(A1[0]) )
                    msd.append( float(A1[1]) )
                    
                msd = np.asarray(msd, dtype=float)
                    
                    
                # Set up bond order parameter
                bond_file = open("img/dat/bond_field_" + str(fr) + ".txt")
                bond = []
                bond_idx = []
                for l1 in bond_file:
                    A1 = l1.split()
                    bond_idx.append( int(A1[0]) )
                    bond.append( float(A1[1]) )
                    
                bond = np.asarray(bond, dtype=float)
                
                
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
                
                
                # Set up virial stress
                vir_file = open("img/dat/local_vir_stress_" + str(fr) + ".txt")
                vir = []
                vir_idx_x = []
                vir_idx_y = []
                for l1 in vir_file:
                    A1 = l1.split()
                    vir_idx_x.append( float(A1[0])*rsepx )
                    vir_idx_y.append( float(A1[1])*rsepy )                
                    vir.append( float(A1[2]) )
                
                vir = np.asarray(vir, dtype=float)
                
                
                # PLOTS
                
                # Set up the plot
                fig = plt.figure()

                s_len = 0.5
                sx_start = 0.1
                sy_start = 0.1
                xint = 0.1
                yint = 0.1
                sxcnt = 0
                sycnt = 0
                
                ax2 = fig.add_axes([sx_start, sy_start, s_len, s_len])
                sxcnt = sxcnt+1
                sycnt = sycnt+1
                        
                ax1 = fig.add_axes( \
                [sx_start, sy_start+sycnt*s_len+sycnt*yint, s_len, s_len])
                
                ax4 = fig.add_axes( \
                [sx_start+sxcnt*s_len+sxcnt*xint, sy_start+sycnt*s_len+sycnt*yint, s_len, s_len])      
                
                cax4 = plt.axes( \
                [sx_start+(sxcnt+1)*s_len+sxcnt*xint+0.02, sy_start+sycnt*s_len+sycnt*yint, 0.02, s_len])

                ax3 = fig.add_axes( \
                [sx_start+sxcnt*s_len+sxcnt*xint, sy_start, s_len, s_len])
                sxcnt = sxcnt+1
                      
                ax6 = fig.add_axes( \
                [sx_start+sxcnt*s_len+sxcnt*xint, sy_start+sycnt*s_len+sycnt*yint, s_len, s_len])
                
                cax6 = plt.axes( \
                [sx_start+(sxcnt+1)*s_len+sxcnt*xint+0.02, sy_start+sycnt*s_len+sycnt*yint, 0.02, s_len])
                
                ax5 = fig.add_axes( \
                [sx_start+sxcnt*s_len+sxcnt*xint, sy_start, s_len, s_len])
                
                cax5 = plt.axes( \
                [sx_start+(sxcnt+1)*s_len+sxcnt*xint+0.02, sy_start, 0.02, s_len])
                sxcnt = sxcnt+1

                ax8 = fig.add_axes( \
                [sx_start+sxcnt*s_len+sxcnt*xint, sy_start+sycnt*s_len+sycnt*yint, s_len, s_len])
                cax8 = plt.axes( \
                [sx_start+(sxcnt+1)*s_len+sxcnt*xint+0.02, sy_start+sycnt*s_len+sycnt*yint, 0.02, s_len])
                
                ax7 = fig.add_axes( \
                [sx_start+sxcnt*s_len+sxcnt*xint, sy_start, s_len, s_len])
                
                cax7 = plt.axes( \
                [sx_start+(sxcnt+1)*s_len+sxcnt*xint+0.02, sy_start, 0.02, s_len])
                sxcnt = sxcnt+1                
            

                
                # Monomer and center of mass positions (ax1)
                ax1.scatter(xcm, ycm, s=2, c=ycm1, cmap=plt.cm.brg, edgecolors='None', alpha=1)
                ax1.scatter(X, Y, s=0.3, c=Yc, cmap=plt.cm.brg, edgecolors='None', alpha=0.7)
                ax1.set_ylabel('y/l')
                ax1.set_title('Cell positions')
                ax1.set_xlim((downlim,uplim))
                ax1.set_ylim((downlim,uplim))
                ax1.yaxis.set_ticks(np.arange(0,uplim,tick_interval))
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
                ax2.yaxis.set_ticks(np.arange(0,uplim,tick_interval))
                ax2.xaxis.set_ticks(np.arange(0,uplim,tick_interval))
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
                ax3.tick_params(axis='both', which='major', labelsize=8)
                ax3.xaxis.set_ticks(np.arange(0,uplim,tick_interval))
                
      
                # Mean square displacement (ax4)
                area_msd = np.pi * (0.5 * msd)**2 
                line4 = ax4.scatter(xcm, ycm, c = msd, s=area_msd, cmap=plt.cm.brg, alpha=0.7, vmin=0, vmax=50)
                ax4.set_title('$MSD/4R^2$')
                ax4.set_xlim((downlim,uplim))
                ax4.set_ylim((downlim,uplim))
                ax4.yaxis.set_ticks(np.arange(0,uplim,tick_interval))
                ax4.xaxis.set_ticks(np.arange(0,uplim,tick_interval))
                ax4.tick_params(axis='both',which='major',labelsize=8)
                
                plt.colorbar(line4,cax=cax4)
                cax4.tick_params(width=1.1,labelsize=6)
    

                
                # Text
                plt.figtext(1.15, 1.3, 't = ' + str(fr*dtSamp))
                
                
                # Local bond order (ax5)
                area_bond = np.pi * (3 * bond)**2 
                line5 = ax5.scatter(xcm, ycm, c = bond, s=area_bond, cmap=plt.cm.brg, alpha=1, vmin=0, vmax=1)
                ax5.set_xlabel('x/l')
                ax5.set_title('$\psi_{4}$')
                ax5.set_xlim((downlim,uplim))
                ax5.set_ylim((downlim,uplim))
                ax5.tick_params(axis='both', which='major', labelsize=8)
                ax5.xaxis.set_ticks(np.arange(0,uplim,tick_interval))
                
                plt.colorbar(line5,cax=cax5)
                cax5.tick_params(width=1.1,labelsize=6)
    
    
                # Number of neighbors per cell (ax6)    
                area_neigh = np.pi * (20*num_neighs)**2 
                line6 = ax6.scatter(xcm, ycm, c=num_neighs, s=num_neighs, cmap=plt.cm.brg, alpha=1, vmin=0, vmax=9)
                ax6.set_title('Num. of neigh.')
                ax6.set_xlim((downlim,uplim))
                ax6.set_ylim((downlim,uplim))
                ax6.yaxis.set_ticks(np.arange(0,uplim,tick_interval))
                ax6.xaxis.set_ticks(np.arange(0,uplim,tick_interval))
                ax6.tick_params(axis='both',which='major',labelsize=8)
                
                plt.colorbar(line6,cax=cax6)
                cax6.tick_params(width=1.1,labelsize=6)
                
                
                # Local pressure map (ax7)
                area_pres = np.pi * (1e-1*pres)**2 
                line7 = ax7.scatter(xcm, ycm, c=pres, s=area_pres, cmap=plt.cm.brg, alpha=1, vmin=0, vmax=20)
                ax7.set_xlabel('x/l')
                ax7.set_title('P')
                ax7.set_xlim((downlim,uplim))
                ax7.set_ylim((downlim,uplim))
                ax7.xaxis.set_ticks(np.arange(0,uplim,tick_interval))
                ax7.tick_params(axis='both', which='major', labelsize=8)
                
                plt.colorbar(line7,cax=cax7)
                cax7.tick_params(width=1.1,labelsize=6)
    
    
                # Local interparticle virial field (ax8)    
                area_vir = np.pi * (1e-3*vir)**2 
                line8 = ax8.scatter(vir_idx_x, vir_idx_y, c=vir, s=area_vir, cmap=plt.cm.brg, alpha=0.7, vmin=0, vmax=3)
                ax8.set_title('Inter. Virial')
                ax8.set_xlim((downlim,uplim))
                ax8.set_ylim((downlim,uplim))
                ax8.yaxis.set_ticks(np.arange(0,uplim,tick_interval))
                ax8.xaxis.set_ticks(np.arange(0,uplim,tick_interval))

                ax8.tick_params(axis='both',which='major',labelsize=8)
                
                plt.colorbar(line6,cax=cax8)
                cax8.tick_params(width=1.1,labelsize=6)
                

                plt.setp(ax3.get_yticklabels(),visible=False)
                plt.setp(ax4.get_xticklabels(),visible=False)
                plt.setp(ax4.get_yticklabels(),visible=False)
                plt.setp(ax5.get_yticklabels(),visible=False)
                plt.setp(ax6.get_xticklabels(),visible=False)
                plt.setp(ax6.get_yticklabels(),visible=False)
                plt.setp(ax7.get_yticklabels(),visible=False)
                plt.setp(ax8.get_xticklabels(),visible=False)
                plt.setp(ax8.get_yticklabels(),visible=False)
                

                        
                    
                plt.savefig('./img/'+'frame-'+'{0:05d}'.format(fr)+'.png',dpi=200,bbox_inches='tight',pad_inches=0.08)
                plt.clf()
                X = []
                Y = []
        
                L_cnt = 0
                
