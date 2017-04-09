#!/usr/bin/python

import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("nsamp", help="Sampling frequency", type=float)
parser.add_argument("M", help="Number of cells", type=int)
parser.add_argument("folder_path", help="folder where the analysis file resides")
parser.add_argument("analysis_file", help="which analysis is being made")
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
 

path = args.folder_path + args.analysis_file + '.txt'
indata = np.loadtxt(path,dtype=float)
outdata = readData(indata)
x = outdata[0]
y = outdata[1]

plt.figure(figsize=(12,10))
fig = plt.subplot(111)

if args.analysis_file == 'num_neigh_dist':
    
    plt.plot(x,y,'-o')
    fig.tick_params(axis='both', which='major', labelsize=20)
    plt.xlabel('n',fontsize=30)
    plt.ylabel('P(n)',fontsize=30)
    plt.title('Avg. num. neigh. dist.',fontsize=20,fontweight='bold')
    
elif args.analysis_file == 'cluster_size_dist':
    
    plt.plot(x,y,'-o')
    fig.tick_params(axis='both', which='major', labelsize=20)
    plt.xlabel('c(s)',fontsize=30)
    plt.ylabel('P(c(s))',fontsize=30)
    plt.title('Avg. cluster size dist.',fontsize=20,fontweight='bold')
    
elif args.analysis_file == 'num_neigh_per_frame':
    
    plt.plot(x/tr,y)
    fig.tick_params(axis='both', which='major', labelsize=20)
    plt.xlabel('$t/\\tau_{R}$',fontsize=30)
    plt.ylabel('n',fontsize=30)
    plt.title('Avg. num. neigh. per frame',fontsize=20,fontweight='bold')

plt.savefig(args.analysis_file+'.png',bbox_inches='tight')


 
 