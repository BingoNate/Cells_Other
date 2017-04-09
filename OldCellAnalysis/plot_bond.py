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

	N = indata.shape[0]

	x = np.zeros(N)
	y = np.zeros(N)
	for i in np.arange(N):
		x[i] = indata[i][0]
		y[i] = indata[i][1]
    
	return x, y, N


# Read the data

# Bond orientational order
indata = np.loadtxt('glob_bond.txt')
outdata = readData( indata )

x = outdata[0]
y = outdata[1]
N = outdata[2]

# Plot the data
plt.figure(figsize=(12,10))

fig = plt.subplot(111)

plt.plot(x,y)
plt.plot(x,y,'o')
plt.xlabel('$\epsilon$',fontsize=30)
plt.ylabel('$\psi_{6}$',fontsize=30)
plt.title('Bond orientational order',fontsize=25,fontweight='bold')
fig.tick_params(axis='both', which='major', labelsize=20)

plt.savefig('glob_bond.png', bbox_inches='tight')

# Susceptibility
indata = np.loadtxt('susceptibility.txt')
outdata = readData( indata )

x = outdata[0]
y = outdata[1]
N = outdata[2]

# Plot the data
plt.figure(figsize=(12,10))

fig = plt.subplot(111)

plt.plot(x,y,linewidth=2.0)
plt.xlabel('$\epsilon$',fontsize=30)
plt.ylabel('$\chi$',fontsize=30)
plt.title('Susceptibility of bond orientation',fontsize=25,fontweight='bold')
fig.tick_params(axis='both', which='major', labelsize=20)

plt.savefig('suscp.png', bbox_inches='tight')

