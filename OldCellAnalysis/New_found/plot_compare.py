#!/usr/bin/python

import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import shutil as shu
import os

parser = argparse.ArgumentParser()
parser.add_argument("nsamp", help="Sampling frequency", type=float)
parser.add_argument("M", help="Number of cells", type=int)
parser.add_argument("cell_file", help="xyz file")
args = parser.parse_args()

# Simulation information
nsamp = args.nsamp
dt = 1e-3*nsamp
Dr = 0.01
tr = 1/Dr 
M = args.M


# Function for reading data from file
def readData( indata ):

	#d = indata.shape[1] # dimension information
	N = indata.shape[0] # number of data points

	x = np.zeros(N)
	y = np.zeros(N)
	for i in np.arange(N):
		x[i] = indata[i][0]
		y[i] = indata[i][1]
	return x, y, N

# Functions for plotting specific analysis modules
def linlog_plotter( x, y, ax, label_name, color_idx ):

	ax[0].plot(x,y,label=label_name,linewidth=1.0,color=color_idx)
	ax[1].loglog(x,y,label=label_name,linewidth=1.0,color=color_idx)

	return ax

def semilogx_plotter( x, y, ax, label_name, color_idx ):

	ax.semilogx(x,y,label=label_name,linewidth=1.0,color=color_idx)

	return ax

def lin_plotter( x, y, ax, label_name, color_idx ):

	ax.plot(x,y,label=label_name,linewidth=1.0,color=color_idx)

	return ax


#
# Read the data
#

# Folder names
#phi = [0.3, 0.4, 0.5, 0.6, 0.7]
#eps = [10, 12, 14]
#Fm = [1, 2]

phi = [0.3, 0.4, 0.5, 0.6, 0.7]
eps = [14]
Fm = [2]

#dir_name = 'eps' + str(eps[0]) + 'Fm' + str(Fm[0])
dir_name = args.cell_file + 'eps' + str(eps[0]) + 'fm' + str(Fm[0])
os.mkdir(dir_name)

folders = []
for P in phi:
    for E in eps:
        for F in Fm:
            folders.append('phi'+str(P)+'eps'+str(E)+'fm'+str(F))


# Analysis modules
analysis = ['msd_post','inter_scattering','dynamic_struct','static_struct', \
'pair_corr','vacf_cm','sp_vacf','bond_corr']


for a in analysis:
    
    x = []
    y = []
    cnt = 0
    
    # Load the data from the files
    for folder in folders:

        path = '../Graphs/' + folder + '/' + a + '.txt'
        indata = np.loadtxt(path,dtype=float)
        outdata = readData( indata )
        x.append( outdata[0] )
        y.append( outdata[1] )
        N = outdata[2]
        cnt += 1
        print path, cnt
        
    # Set general graph properties
    plt.figure(figsize=(20,17))
    colors = plt.cm.spectral(np.linspace(0,1,cnt))
    if a == 'msd_post':
        fig, ax = plt.subplots(1,2)
    else:
        fig, ax = plt.subplots(1,1)
    
    # Plot the data
    for i in np.arange(cnt):
        if a == 'msd_post':
            ax = linlog_plotter(x[i]/tr,y[i],ax,folders[i],colors[i])
        elif a == 'inter_scattering':
            ax = semilogx_plotter(x[i]/tr,y[i],ax,folders[i],colors[i])
        elif a == 'vacf_cm':
            ax = lin_plotter(x[i]/tr,y[i],ax,folders[i],colors[i])
        else:
            ax = lin_plotter(x[i],y[i],ax,folders[i],colors[i])
            
    # Customize the plots
    if a == 'msd_post':
        ax[1].legend(bbox_to_anchor=(1.05,1), loc=2, borderaxespad=0., prop={'size': 9})
        ax[0].tick_params(axis='both',which='major',labelsize=11)
        ax[1].tick_params(axis='both',which='major',labelsize=11)
    else:
        ax.legend(bbox_to_anchor=(1.05,1), loc=2, borderaxespad=0., prop={'size': 9})
        ax.tick_params(axis='both',which='major',labelsize=15)
        
    if a == 'static_struct':
        ax.set_title('Static structure factor')
        ax.set_xlabel('k*R')
        ax.set_ylabel('S(k)')
    elif a == 'inter_scattering':
        ax.set_title('Intermediate scattering function')
        ax.set_xlabel('$t/\\tau_{R}$')
        ax.set_ylabel('F(k,t)') 
    elif a == 'dynamic_struct':
        ax.set_title('Dynamic structure factor')
        ax.set_xlabel('w')
        ax.set_ylabel('DS(k,w)')
    elif a == 'msd_post':
        ax[0].set_xlabel('$t/\\tau_{R}$')
        ax[0].set_ylabel('$\Delta r^2/4R^2$')
        ax[1].set_xlabel('$t/\\tau_{R}$')
    elif a == 'pair_corr':
        ax.set_title('Pair correlation function')
        ax.set_xlabel('r/2R')
        ax.set_ylabel('g(r)') 
        ax.set_xlim((0,9)) 
    elif a == 'sp_vacf':
        ax.set_title('Spatial velocity correlation')
        ax.set_xlabel('r/2R')
        ax.set_ylabel('$C_{vv}(r)$')
    elif a == 'vacf_cm':
        ax.set_title('Velocity autocorrelation')
        ax.set_xlabel('$t/\\tau_{R}$')
        ax.set_ylabel('$C_{vv}(t)$') 
    elif a == 'bond_corr':
        ax.set_title('Spatial orientational order correlation')
        ax.set_xlabel('r/2R')
        ax.set_ylabel('$g_6(r)$') 
        
    plt.savefig(a+'.eps', bbox_inches='tight')
    
    shu.move(a+'.eps',dir_name)

