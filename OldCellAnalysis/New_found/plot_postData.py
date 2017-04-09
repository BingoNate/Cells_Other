#!/usr/bin/python

import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("nsamp", help="Sampling frequency", type=float)
parser.add_argument("M", help="Number of cells", type=int)
parser.add_argument("cell_file", help="xyz file")
args = parser.parse_args()

nsamp = args.nsamp
dt = 1e-3
dtSamp = dt*nsamp
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

	#ax[1].semilogy(x,y,label=label_name,linewidth=1.0,color=color_idx)
	ax[1].loglog(x,y,label=label_name,linewidth=1.0,color=color_idx)


	return ax

def semilogx_plotter( x, y, ax, label_name, color_idx ):

	ax.semilogx(x,y,label=label_name,linewidth=1.0,color=color_idx)

	return ax

def lin_plotter( x, y, ax, label_name, color_idx ):

	ax.plot(x,y,label=label_name,linewidth=1.0,color=color_idx)

	return ax


# Velocity autocorrelation
indata = np.loadtxt('vacf_cm.txt',dtype=float)
outdata = readData( indata )
t = outdata[0]
v = outdata[1]
N = outdata[2]

plt.figure(figsize=(12,10))

fig = plt.subplot(111)
plt.plot(t/tr,v,linewidth=2.0)
plt.ylabel('<v(t)v(0)>',fontsize=20)
plt.title('Velocity autocorrelation',fontsize=30,fontweight='bold')
fig.tick_params(axis='both', which='major', labelsize=20)

plt.savefig('vacf_cm.png',bbox_inches='tight')


# Spatial velocity correlation
indata = np.loadtxt('sp_vacf.txt',dtype=float)
outdata = readData( indata )
rd = outdata[0]
cv = outdata[1]
N = outdata[2]

plt.figure(figsize=(12,10))

fig = plt.subplot(111)
plt.plot(rd,cv,linewidth=2.0)
plt.title('Spatial velocity correlation',fontsize=30,fontweight='bold')
plt.xlabel('$\\frac{r}{2R}$',fontsize=30)
plt.ylabel('$C_{vv}$',fontsize=30)
fig.tick_params(axis='both', which='major', labelsize=20)

plt.savefig('vacf_sp.png',bbox_inches='tight')


# Pair correlation function
indata = np.loadtxt('pair_corr.txt',dtype=float)
outdata = readData( indata )
k = outdata[0]
s = outdata[1]
N = outdata[2]

plt.figure(figsize=(12,10))

fig = plt.subplot(111)
plt.plot(k,s,linewidth=2.0)
plt.xlabel('$\\frac{r}{2R}$',fontsize=30)
plt.ylabel('g(r)',fontsize=30)
plt.title('Pair correlation function',fontsize=30,fontweight='bold')
fig.tick_params(axis='both', which='major', labelsize=20)

plt.savefig('pair_corr.png',bbox_inches='tight')


# Static structure factor
indata = np.loadtxt('static_struct.txt',dtype=float)
outdata = readData( indata )
k = outdata[0]
s = outdata[1]
N = outdata[2]

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
outdata = readData( indata )
delt = outdata[0]
f = outdata[1]
N = outdata[2]

plt.figure(figsize=(12,10))

fig = plt.subplot(111)
plt.plot(delt/tr,f,linewidth=2.0)
plt.xscale("log")
plt.xlabel('$\\Delta t$',fontsize=30)
plt.ylabel('F(k,t)',fontsize=30)
plt.title('Intermediate scattering',fontsize=30,fontweight='bold')
fig.tick_params(axis='both', which='major', labelsize=20)

plt.savefig('inter_scattering.png',bbox_inches='tight')



# Dynamic structure factor
indata = np.loadtxt('dynamic_struct.txt',dtype=float)
outdata = readData( indata )
w = outdata[0]
DS = outdata[1]
N = outdata[2]

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
outdata = readData( indata )
rd = outdata[0]
cv = outdata[1]
N = outdata[2]

plt.figure(figsize=(12,10))

fig = plt.subplot(111)
plt.plot(rd,cv,linewidth=2.0)
plt.title('Spatial bond correlation',fontsize=30,fontweight='bold')
plt.xlabel('$\\frac{r}{2R}$',fontsize=30)
plt.ylabel('$C_{bb}$',fontsize=30)
fig.tick_params(axis='both', which='major', labelsize=20)

plt.savefig('bond_corr.png',bbox_inches='tight')



# Mean square displacment
indata = np.loadtxt('msd_post.txt',dtype=float)
outdata = readData( indata )
tmsd = outdata[0]
vmsd = outdata[1]
N = outdata[2]

plt.figure(figsize=(12,10))

fig = plt.subplot(211)
plt.plot(tmsd/tr,vmsd)
plt.ylabel('$<\\Delta r^{2}>/4R^2$',fontsize=30)
plt.title('Mean square displacement',fontsize=35,fontweight='bold')
fig.tick_params(axis='both', which='major', labelsize=20)

fig = plt.subplot(212)
plt.loglog(tmsd/tr,vmsd)
plt.xlabel('$\\Delta t$',fontsize=30)
plt.ylabel('$<\\Delta r^{2}>/4R^2$',fontsize=30)
fig.tick_params(axis='both', which='major', labelsize=20)

plt.savefig('msd_post.png',bbox_inches='tight')


## Plot exponent coefficients
#t_log = np.log(t[1:np.size(tmsd)])
#msd_log = np.log(v[1:np.size(vmsd)])
#egim_offset = 10
#N = 90
#egim = np.zeros(N)
#t_egim = np.zeros(N)
#i = 0
#for i in np.arange(N):
#    coef = np.polyfit(t_log[i:i+egim_offset],msd_log[i:i+egim_offset],1)
#    egim[i] = coef[0]
#    t_egim[i] = dtSamp*i
#
#plt.figure(figsize=(12,10))
#fig = plt.subplot(111)
#plt.plot(t_egim,egim,'o')
#plt.title('Exponent fits within a window of %s' %(egim_offset), fontsize=20,fontweight='bold')
#plt.xlabel('t',fontsize=20)
#plt.ylabel('$\\alpha$',fontsize=20)
#fig.tick_params(axis='both', which='major', labelsize=20)
#
#plt.savefig('exponents.png',bbox_inches='tight')

