#!/usr/local/bin/python

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import subprocess
import sh
import glob
import os
import re

# Simulation information
sample_data = 1
nsamp = 10000*sample_data
dt = 1e-3*nsamp

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

indata = np.loadtxt("bond_corr.txt",dtype=float)
outdata = readData( indata )
x = []
y = []
x.append( outdata[0] )
y.append( outdata[1] )
N = outdata[2]

plt.figure(figsize=(16,13))
fig, ax = plt.subplots(1,1)
ax.plot(x,y,'*')
