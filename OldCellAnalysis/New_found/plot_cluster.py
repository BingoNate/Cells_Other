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
tick_interval = 100

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
                
                xcm = np.asarray(xcm,dtype=float)
                ycm = np.asarray(ycm,dtype=float)
                
                
                # Set up positions of cells making up the cluster 
                cluster_file = open("img/dat/cluster_list_" +str(fr) + ".txt")
                node_idx = []
                edge_idx = []
                for l1 in cluster_file:
                    A1 = l1.split()
                    node_idx.append( int(A1[0]) )
                    for nidx in np.arange( np.size(A1) ):
                        edge_idx.append( int(A1[nidx]) )
                    edge_idx.append(-1)

                clusters = dict()
                node_cnt = -1
                edges = []
                for edge in edge_idx:
                    if edge == -1:
                        node_cnt = node_cnt+1
                        clusters[node_idx[node_cnt]] = edges
                        edges = []
                    else:
                        edges.append(edge)

                    
                        
                
                # PLOTS
                
                # Set up the plot
                fig = plt.figure()
                
                ax1 = fig.add_axes([0.1,0.1,0.5,0.5])
                ax2 = fig.add_axes([0.65,0.1,0.5,0.5])
                                
                
                # Monomer positions (ax1)
                ax1.scatter(xcm, ycm, s=2, c=ycm1, cmap=plt.cm.brg, edgecolors='None', alpha=1)
                ax1.scatter(X, Y, s=0.3, c=Yc, cmap=plt.cm.brg, edgecolors='None', alpha=1)
                ax1.set_xlabel('x/l')
                ax1.set_ylabel('y/l')
                ax1.set_title('Cell positions')
                ax1.set_xlim((downlim,uplim))
                ax1.set_ylim((downlim,uplim))
                ax1.set_aspect('equal')
                ax1.xaxis.set_ticks(np.arange(0,uplim,tick_interval))
                ax1.yaxis.set_ticks(np.arange(0,uplim,tick_interval))
                ax1.tick_params(axis='both', which='major', labelsize=8)
                
                cmap = plt.cm.Set1
                plot_cnt = -1
                for key, value in clusters.iteritems():
                    plot_cnt = plot_cnt+1
                    #print key, value
                    value = np.asarray(value,dtype=int)
                    
#                    area = np.pi * (0.01 * np.size(value))**2 
#                    ax2.scatter(xcm[value], ycm[value], s=area, color=cmap(plot_cnt/float(node_cnt+1)))

                    if np.size(value) > 4:
                        ax2.scatter(xcm[value], ycm[value], color=cmap(plot_cnt/float(node_cnt+1)))
                    else:
                        ax2.scatter(xcm[value], ycm[value], color='r', s=2)
                    
                ax2.set_xlabel('x/l')
                ax2.set_title('Clusters')
                ax2.set_xlim((downlim,uplim))
                ax2.set_ylim((downlim,uplim))
                ax2.set_aspect('equal')
                plt.setp(ax2.get_yticklabels(),visible=False)
                ax2.tick_params(axis='both', which='major', labelsize=8)
                ax2.xaxis.set_ticks(np.arange(0,uplim,tick_interval))
                
                
                # Text
                plt.figtext(0.6, 0.65, 't = ' + str(fr*dtSamp))

            
            
                    
                plt.savefig('./cluster/'+'frame-'+'{0:05d}'.format(fr)+'.png',dpi=200,bbox_inches='tight',pad_inches=0.08)
                plt.clf()
                X = []
                Y = []
        
                L_cnt = 0
                
