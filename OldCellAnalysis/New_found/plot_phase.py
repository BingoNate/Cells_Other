#!/usr/bin/python

import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()
parser.add_argument("nsamp", help="Sampling frequency", type=float)
parser.add_argument("M", help="Number of cells", type=int)
parser.add_argument("cell_file", help="Analysis file")
parser.add_argument("analysis_module", help="Analysis module to use")
args = parser.parse_args()

# Simulation information
nsamp = args.nsamp
dt = 1e-3*nsamp
Dr = 0.01
tr = 1/Dr 
M = args.M


# Function to read the data from text file
def readData( indata ):
    
    d = indata.shape[1]     # Dimension information
    N = indata.shape[0]     # Number of data points
    
    o1 = np.zeros(N)
    o2 = np.zeros(N)
    o3 = np.zeros(N)
    o4 = np.zeros(N)
    
    for i in np.arange(N):
        o1[i] = indata[i][0]
        o2[i] = indata[i][1]
        o3[i] = indata[i][2]
        o4[i] = indata[i][3]
        
    return o1, o2, o3, o4


# Function to put the data in correct format
def formData( zi, X, Y, Z, I ):
    
    x = []
    y = []
    z = []
    cnt = -1
    for i in I:
        cnt += 1
        if i == zi:
            x.append( X[cnt] )
            y.append( Y[cnt] )
            z.append( Z[cnt] )
            
    return x, y, z
    
    
# Function to append the data into lists
def appData( x, y, z, temp ):
    
    x.append( temp[0] )
    y.append( temp[1] )
    z.append( temp[2] )
    
    return x, y, z



# Folder names
P = [0.3, 0.4, 0.5, 0.6, 0.7]
E = [10, 12, 14, -1, -1]
F = [1, 2, -1, -1, -1]

c = args.analysis_module

if c == 'stress':
    sfac = 1
    vmin_val = -4
    vmax_val = -2
    path = 'phase/phase_stress'
    
elif c == 'max':
    sfac = 0.5
    vmin_val = 5
    vmax_val = 40
    path = 'phase/phase_max'
    
elif c == 'bond':
    sfac = 100
    vmin_val = 0.01
    vmax_val = 0.05
    path = 'phase/phase_bond'
    
elif c == 'vel':
    sfac = 100
    vmin_val = 0.01
    vmax_val = 0.1
    path = 'phase/phase_speed'
    
elif c == 'neigh_chn':
    sfac = 0.1
    vmin_val = 0.1
    vmax_val = 105
    path = 'phase/phase_neigh_chn'
    
elif c == 'neigh_num':
    sfac = 1
    vmin_val = 3
    vmax_val = 6
    path = 'phase/phase_num_neigh'
    
elif c == 't1':
    sfac = 1
    vmin_val = 0.001
    vmax_val = 26
    path = 'phase/phase_t1'      
    
elif c == 'pressure':
    sfac = 1
    vmin_val = 2
    vmax_val = 2.4
    path = 'phase/phase_pressure'   

# Read the data
indata = np.loadtxt(args.cell_file,dtype=float)
outdata = readData( indata )
phi = outdata[0]
Fm = outdata[1]
eps = outdata[2]
ob = outdata[3]

# Plot properties
fig = plt.figure()
fi, axi = plt.subplots(nrows=3,ncols=5)
fi.subplots_adjust(bottom=0.1, top=0.9, left=0.125, right=0.9, hspace=0.3, wspace=0.2)

# Set up the tick labels
plt.setp([a.get_xticklabels() for a in axi[0, 1:]], visible=False)
plt.setp([a.get_xticklabels() for a in axi[1, 1:]], visible=False)
plt.setp([a.get_xticklabels() for a in axi[2, 1:]], visible=False)
plt.setp([a.get_yticklabels() for a in axi[:, 1]], visible=False)
plt.setp([a.get_yticklabels() for a in axi[:, 2]], visible=False)
plt.setp([a.get_yticklabels() for a in axi[:, 3]], visible=False)
plt.setp([a.get_yticklabels() for a in axi[:, 4]], visible=False)

# Set up the tick sizes
plt.setp([a.get_yticklabels() for a in axi[0, :]], fontsize=5)
plt.setp([a.get_yticklabels() for a in axi[1, :]], fontsize=5)
plt.setp([a.get_yticklabels() for a in axi[2, :]], fontsize=5)
plt.setp([a.get_xticklabels() for a in axi[1, :]], fontsize=5)
plt.setp(axi[0,0].get_xticklabels(),fontsize=5)
plt.setp(axi[1,0].get_xticklabels(),fontsize=5)
plt.setp(axi[2,0].get_xticklabels(),fontsize=5)

axi[0,0].xaxis.set_ticks([0.3,0.5,0.7])
axi[1,0].xaxis.set_ticks([0.3,0.5,0.7])
axi[2,0].xaxis.set_ticks([10,12,14])


# Plot
row_cnt = 0
col_cnt = -1

# As a function of eps (phi vs. Fm)
for e in E:
    
    col_cnt += 1
    x = []
    y = []
    z = []
    
    temp = formData(e,phi,Fm,ob,eps)
    x.append( temp[0] )
    y.append( temp[1] )
    z.append( temp[2] )
    z = np.asarray(z, dtype=float)
    
    area = np.pi * (sfac*z)**2
    line1 = axi[row_cnt, col_cnt].scatter(x,y,c=z,s=area,vmin=vmin_val,vmax=vmax_val)
    
    if col_cnt == 3:
        cax = plt.axes([0.93, 0.1, 0.02, 0.8])
        plt.colorbar(line1,cax=cax)
        cax.tick_params(width=1.1,labelsize=5)
#        bbox_props = dict(boxstyle="rarrow,pad=0.01", fc="cyan", ec="b", lw=0.5)
#        t = axi[row_cnt, col_cnt].text(-12, 4.5, "$+ \epsilon$", ha="center", va="center", rotation=0, \
#        bbox=bbox_props, fontsize=8)
        
    if col_cnt == 2:
        axi[row_cnt, col_cnt].set_title("x = $\phi$, y = $F_m$",fontsize=8)
    
      
# As a function of Fm (phi vs. eps)
row_cnt = 1
col_cnt = -1
for f in F:
    
    col_cnt += 1
    x = []
    y = []
    z = []
    
    temp = formData(f,phi,eps,ob,Fm)
    x.append( temp[0] )
    y.append( temp[1] )
    z.append( temp[2] )
    z = np.asarray(z, dtype=float)

    area = np.pi * (sfac*z)**2
    axi[row_cnt, col_cnt].scatter(x,y,c=z,s=area,vmin=vmin_val,vmax=vmax_val) 
    
#    if col_cnt == 3:
#        bbox_props = dict(boxstyle="rarrow,pad=0.01", fc="cyan", ec="b", lw=0.5)
#        t = axi[row_cnt, col_cnt].text(-12, 4.5, "$+ F_m$", ha="center", va="center", rotation=0, \
#        bbox=bbox_props, fontsize=6)
    
    if col_cnt == 2:
        axi[row_cnt, col_cnt].set_title("x = $\phi$, y = $\epsilon$",fontsize=8)
        
        
# As a function of phi (eps vs. Fm)
row_cnt = 2
col_cnt = -1 
for p in P:
    
    col_cnt += 1
    x = []
    y = []
    z = []
    
    temp = formData(p,eps,Fm,ob,phi)
    x.append( temp[0] )
    y.append( temp[1] )
    z.append( temp[2] )
    z = np.asarray(z, dtype=float)

    area = np.pi * (sfac*z)**2
    axi[row_cnt, col_cnt].scatter(x,y,c=z,s=area,vmin=vmin_val,vmax=vmax_val)
    
#    if col_cnt == 3:
#        bbox_props = dict(boxstyle="rarrow,pad=0.01", fc="cyan", ec="b", lw=0.5)
#        t = axi[row_cnt, col_cnt].text(-12, 4.5, "$+ \phi$", ha="center", va="center", rotation=0, \
#        bbox=bbox_props, fontsize=6)
        
    if col_cnt == 2:
        axi[row_cnt, col_cnt].set_title("x = $\epsilon$, y = $F_m$",fontsize=8)
        
   
plt.savefig(path + '.png',bbox_inches='tight')

                
