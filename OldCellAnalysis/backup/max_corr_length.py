#!/usr/bin/python

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

# Simulation information
sample_data = 1
nsamp = 10000*sample_data
dt = 1e-3*nsamp

# Function to read the data from text file
def readData( indata ):
    
    d = indata.shape[1] # dimension information
    N = indata.shape[0] # number of data points
        
    if d == 2:
        x = np.zeros(N)
        y = np.zeros(N)
        for i in np.arange(N):
            x[i] = indata[i][0]
            y[i] = indata[i][1]
        return x, y, d, N
    elif d == 1:
        x = np.zeros(N)
        for i in np.arange(N):
            x[i] = indata[i][0]
        return x, d, N

# Read the data

# First one
indata = np.loadtxt('max_corr_length.txt',dtype=float)
outdata = readData( indata )

x = outdata[0]
y = outdata[1]
d = outdata[2]
N = outdata[3]


# Plot the data
plt.figure(figsize=(12,10))

fig = plt.subplot(111)

plt.plot(x,y,'o')
plt.plot(x,y,linewidth=2.0)
plt.xlabel('$\epsilon/k_{B}T$',fontsize=30)
plt.ylabel('$\eta_s$ (Cells)',fontsize=30)
plt.title('Maximum cluster size',fontsize=30,fontweight='bold')
#plt.xlim((25,41))
#plt.ylim((25,41))
fig.tick_params(axis='both', which='major', labelsize=30)

plt.savefig('max_cluster_size.eps', bbox_inches='tight')
