#!/usr/bin/python

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

nsamp = 50000
dt = 1e-3*nsamp
Dr = 0.01
tr = 1/Dr
M = 600

sfac = 1
vmin_val = 0
vmax_val = 6

indata = np.loadtxt('avg_neigh_change.txt',dtype=float)
N = indata.shape[0]
phi = np.zeros(N)
eps = np.zeros(N)
t1 = np.zeros(N)

for i in np.arange(N):
	fm = indata[i][1]
	if fm == 1:
		phi[i] = indata[i][0]
		eps[i] = indata[i][2]
		t1[i] = indata[i][3]

plt.figure()
area = np.pi * (sfac*t1)**2
line = plt.scatter(phi,eps,c=t1,s=area,vmin=vmin_val,vmax=vmax_val)
plt.title('Avg. number of neigh. changes',fontsize=20,fontweight='bold')
plt.xlabel('$\phi$',fontsize=30)
plt.ylabel('$\epsilon$',fontsize=30)
plt.colorbar(line)
plt.savefig('neigh_change_plot.png',bbox_inches='tight')



