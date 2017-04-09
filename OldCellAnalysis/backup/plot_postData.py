#!/usr/bin/python

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

sample_data = 1
nsamp = 10000*sample_data
dt = 1e-3
dtSamp = dt*nsamp


# Velocity autocorrelation
indata = np.loadtxt('vacf_cm.txt',dtype=float)

d = indata.shape[1] # dimensions
N = indata.shape[0] # number of data points

t = np.zeros(N)
for i in np.arange(N):
    t[i] = indata[i][0]

v = np.zeros(N)
for i in np.arange(N):
    v[i] = indata[i][1]

plt.figure(figsize=(12,10))

fig = plt.subplot(211)
plt.plot(t,v,linewidth=2.0)
plt.ylabel('<v(t)v(0)>/<v(0)v(0)>',fontsize=20)
plt.title('Velocity autocorrelation',fontsize=30,fontweight='bold')
fig.tick_params(axis='both', which='major', labelsize=20)

fig = plt.subplot(212)
plt.loglog(t,v,linewidth=2.0)
plt.xlabel('t',fontsize=30)
plt.ylabel('<v(t)v(0)>/<v(0)v(0)>',fontsize=30)
fig.tick_params(axis='both', which='major', labelsize=20)

plt.savefig('vacf_cm.png',bbox_inches='tight')


# Spatial velocity correlation
indata = np.loadtxt('sp_vacf.txt',dtype=float)

d = indata.shape[1] # dimensions
N = indata.shape[0] # number of data points

rd = np.zeros(N)
for i in np.arange(N):
    rd[i] = indata[i][0]

cv = np.zeros(N)
for i in np.arange(N):
    cv[i] = indata[i][1]

plt.figure(figsize=(12,10))

fig = plt.subplot(111)
plt.plot(rd,cv,linewidth=2.0)
plt.title('Spatial velocity correlation',fontsize=30,fontweight='bold')
plt.xlabel('$\\frac{r}{2R}$',fontsize=30)
plt.ylabel('$C_{vv}$',fontsize=30)
fig.tick_params(axis='both', which='major', labelsize=20)

plt.savefig('vacf_sp.png',bbox_inches='tight')


# Pair distribution function
indata = np.loadtxt('pair_corr.txt',dtype=float)

d = indata.shape[1] # dimensions
N = indata.shape[0] # number of data points

k = np.zeros(N)
for i in np.arange(N):
    k[i] = indata[i][0]

s = np.zeros(N)
for i in np.arange(N):
    s[i] = indata[i][1]

plt.figure(figsize=(12,10))

fig = plt.subplot(111)
plt.plot(k,s,linewidth=2.0)
plt.xlabel('$\\frac{r}{2R}$',fontsize=30)
plt.ylabel('g(r)',fontsize=30)
plt.title('Pair distribution function',fontsize=30,fontweight='bold')
fig.tick_params(axis='both', which='major', labelsize=20)

plt.savefig('pair_corr.png',bbox_inches='tight')



# Static structure factor
indata = np.loadtxt('static_struct.txt',dtype=float)

d = indata.shape[1] # dimensions
N = indata.shape[0] # number of data points

k = np.zeros(N)
for i in np.arange(N):
    k[i] = indata[i][0]

s = np.zeros(N)
for i in np.arange(N):
    s[i] = indata[i][1]

plt.figure(figsize=(12,10))

fig = plt.subplot(111)
plt.plot(k,s,linewidth=2.0)
plt.xlabel('kR',fontsize=30)
plt.ylabel('S(k)',fontsize=30)
fig.tick_params(axis='both', which='major', labelsize=20)
plt.title('Static structure factor',fontsize=30,fontweight='bold')

plt.savefig('static_struct.png',bbox_inches='tight')



# Intermediate scattering function
indata = np.loadtxt('inter_scattering.txt',dtype=float)

d = indata.shape[1] # dimensions
N = indata.shape[0] # number of data points

delt = np.zeros(N)
for i in np.arange(N):
    delt[i] = indata[i][0]

f = np.zeros(N)
for i in np.arange(N):
    f[i] = indata[i][1]

plt.figure(figsize=(12,10))

fig = plt.subplot(111)
plt.plot(delt,f,linewidth=2.0)
plt.xscale("log")
plt.xlabel('t',fontsize=30)
plt.ylabel('F(k,t)',fontsize=30)
plt.title('Intermediate scattering',fontsize=30,fontweight='bold')
fig.tick_params(axis='both', which='major', labelsize=20)


plt.savefig('inter_scattering.png',bbox_inches='tight')



# Dynamic structure factor
indata = np.loadtxt('dynamic_struct.txt',dtype=float)

d = indata.shape[1] # dimensions
N = indata.shape[0] # number of data points

w = np.zeros(N)
for i in np.arange(N):
    w[i] = indata[i][0]

DS = np.zeros(N)
for i in np.arange(N):
    DS[i] = indata[i][1]

plt.figure(figsize=(12,10))

fig = plt.subplot(111)
plt.plot(w,DS,linewidth=2.0)
#plt.xlim((-1,1))
plt.xlabel('w',fontsize=30)
plt.ylabel('S(k,w)/S(k)',fontsize=30)
plt.title('Dynamic structure factor',fontsize=30,fontweight='bold')
fig.tick_params(axis='both', which='major', labelsize=20)

plt.savefig('dynamic_struct.png',bbox_inches='tight')



# Bond correlation
indata = np.loadtxt('bond_corr.txt',dtype=float)

d = indata.shape[1] # dimensions
N = indata.shape[0] # number of data points

rd = np.zeros(N)
for i in np.arange(N):
    rd[i] = indata[i][0]

cv = np.zeros(N)
for i in np.arange(N):
    cv[i] = indata[i][1]

plt.figure(figsize=(12,10))

fig = plt.subplot(111)
plt.plot(rd,cv,linewidth=2.0)
plt.title('Spatial bond correlation',fontsize=30,fontweight='bold')
plt.xlabel('$\\frac{r}{2R}$',fontsize=30)
plt.ylabel('$C_{bb}$',fontsize=30)
fig.tick_params(axis='both', which='major', labelsize=20)

plt.savefig('bond_corr.png',bbox_inches='tight')


# Directional spatial velocity correlation
indata = np.loadtxt('sp_dir_vacf.txt',dtype=float)

d = indata.shape[1] # dimensions
N = indata.shape[0] # number of data points

r = np.zeros(N)
for i in np.arange(N):
    r[i] = indata[i][0]

par = np.zeros(N)
for i in np.arange(N):
    par[i] = indata[i][1]

perp = np.zeros(N)
for i in np.arange(N):
    perp[i] = indata[i][2]

#data_son = np.zeros(N)
#data_son = np.append(par,perp,cv)

plt.figure(figsize=(12,10))

fig = plt.subplot(111)
plt.plot(r,par,label='Parallel',linewidth=2.0)
plt.plot(r,perp,label='Perpendicular',linewidth=2.0)
plt.xlabel('$\\frac{r}{2R}$',fontsize=30)
plt.ylabel('$C_{vv}$',fontsize=30)
plt.legend(loc='upper left')
#plt.ylim((-0.1,1))
fig.tick_params(axis='both', which='major', labelsize=20)

plt.savefig('vacf_sp_dir.png',bbox_inches='tight')


