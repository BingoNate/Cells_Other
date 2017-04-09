#!/usr/bin/python

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

sample_data = 1
nsamp = 10000*sample_data
dt = 1e-3
dtSamp = dt*nsamp

# Mean square displacment
indata = np.loadtxt('msd_post.txt',dtype=float)

d = indata.shape[1] # dimensions
N = indata.shape[0] # number of data points

t = np.zeros(N)
for i in np.arange(N):
    t[i] = indata[i][0]

v = np.zeros(N)
for i in np.arange(N):
    v[i] = indata[i][1]/((2*20/np.pi)*(2*20/np.pi))

plt.figure(figsize=(12,10))

fig = plt.subplot(211)
plt.plot(t,v)
plt.ylabel('$<\\Delta r^{2}>/4R^2$',fontsize=30)
plt.title('Mean square displacement',fontsize=35,fontweight='bold')
fig.tick_params(axis='both', which='major', labelsize=20)

fig = plt.subplot(212)
plt.loglog(t,v)
plt.xlabel('$\\Delta t$',fontsize=30)
plt.ylabel('$<\\Delta r^{2}>/4R^2$',fontsize=30)
fig.tick_params(axis='both', which='major', labelsize=20)

plt.savefig('msd_post.png',bbox_inches='tight')


# Plot exponent coefficients
indata = np.loadtxt('msd_post.txt',dtype=float)

d = indata.shape[1] # dimensions
N = indata.shape[0] # number of data points

t = np.zeros(N)
for i in np.arange(N):
    t[i] = indata[i][0]

v = np.zeros(N)
for i in np.arange(N):
    v[i] = indata[i][1]

t_log = np.log(t[1:np.size(t)])
msd_log = np.log(v[1:np.size(v)])
egim_offset = 10
N = 90
egim = np.zeros(N)
t_egim = np.zeros(N)
i = 0
for i in np.arange(N):
    coef = np.polyfit(t_log[i:i+egim_offset],msd_log[i:i+egim_offset],1)
    egim[i] = coef[0]
    t_egim[i] = dtSamp*i

plt.figure(figsize=(12,10))
fig = plt.subplot(111)
plt.plot(t_egim,egim,'o')
plt.title('Exponent fits within a window of %s' %(egim_offset), fontsize=20,fontweight='bold')
plt.xlabel('t',fontsize=20)
plt.ylabel('$\\alpha$',fontsize=20)
fig.tick_params(axis='both', which='major', labelsize=20)

plt.savefig('exponents.png',bbox_inches='tight')

