#!/usr/bin/python

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

msd = np.loadtxt('msd.txt',dtype=float)
info = np.loadtxt('rawinfo.txt',dtype=float)
Dt = info[0]
Dr = info[1]
v0 = info[2]
Tr = 1/Dr

dt = 1e-3
sampleMSD = 10
totalMSD = msd.shape[0]
dtMSD = dt*sampleMSD
t = dtMSD*np.arange(totalMSD)

msd_analytic = np.zeros(totalMSD)
Dt = Dt*np.sqrt(dt)
for i in np.arange(totalMSD):
    msd_analytic[i] = 4*Dt*t[i] + 2*v0*v0*Tr*Tr*(t[i]/Tr + np.exp(-t[i]/Tr) - 1)

# Plot MSD
plt.figure(figsize=(12,10))

plt.subplot(211)
plt.plot(t,msd,label='result')
plt.plot(t,msd_analytic,label='analytic')
plt.legend(loc='upper left')
plt.ylabel('MSD',fontsize=20)
plt.title('Mean square displacement',fontsize=20,fontweight='bold')

plt.subplot(212)
plt.loglog(t,msd,label='result')
plt.loglog(t,msd_analytic,label='analytic')
plt.legend(loc='upper left')
plt.xlabel('Lag time',fontsize=20)
plt.ylabel('MSD',fontsize=20)

plt.savefig('msd_detailed.png',bbox_inches='tight')


# Plot exponent coefficients
t_log = np.log(t[1:np.size(t)])
msd_log = np.log(msd[1:np.size(msd)])
egim_offset = 10
N = 200
#N = totalMSD-10
egim = np.zeros(N)
t_egim = np.zeros(N)
i = 0
#while i < N:
for i in np.arange(N):
    coef = np.polyfit(t_log[i:i+egim_offset],msd_log[i:i+egim_offset],1)
    egim[i] = coef[0]
    t_egim[i] = dtMSD*i
#i = i+egim_offset

plt.figure(figsize=(12,10))
plt.plot(t_egim,egim,'o')
plt.title('Exponent fits within a window of %s' %(egim_offset), fontsize=20,fontweight='bold')
plt.xlabel('t',fontsize=20)
plt.ylabel('$\\alpha$',fontsize=20)

plt.savefig('diff.png',bbox_inches='tight')
