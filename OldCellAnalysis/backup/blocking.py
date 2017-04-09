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
indata = np.loadtxt('block_displacement.txt',dtype=float)
outdata = readData( indata )

if len(outdata) == 3:
    x = outdata[0]
    d = outdata[1]
    N = outdata[2]
elif len(outdata) == 4:
    x = outdata[0]
    y = outdata[1]
    d = outdata[2]
    N = outdata[3]


# Plot the data
plt.figure(figsize=(12,10))

fig = plt.subplot(111)

plt.plot(x[10:],y[10:],'o')
#plt.xlim((-1,5010))
#plt.ylim((-1,1.5))
plt.xlabel('Block length',fontsize=30)
plt.ylabel('BSE',fontsize=30)
plt.title('Standard error of average displacement with blocking',fontsize=20,fontweight='bold')
fig.tick_params(axis='both', which='major', labelsize=20)

plt.savefig('block_error.png', bbox_inches='tight')
