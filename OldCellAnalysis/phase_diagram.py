#!/usr/bin/python

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


# Simulation information
sample_data = 1
nsamp = 10000*sample_data
dt = 1e-3*nsamp
R = 20/np.pi


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
P = [0.4, 0.5, 0.6, 0.7, 0.8, 0.95]
E = [0.05, 0.5, 2, 4, 6, 8, 10]
F = [0, 0.5, 1, 1.5, 2, 2.5, 3]

# Read the data
indata = np.loadtxt('stress_int.txt',dtype=float)
outdata = readData( indata )
phi = outdata[0]
Fm = outdata[1]
eps = outdata[2]
ob = outdata[3]

# Plot properties
fig = plt.figure()
fi, axi = plt.subplots(nrows=3,ncols=7)
fi.subplots_adjust(bottom=0.1, top=0.9, left=0.125, right=0.9, hspace=0.3, wspace=0.2)
#axi[0,0].set_title('$\phi$ - $F_{m}$')
#axi[1,0].set_title('$\phi$ - $\epsilon$')
#axi[2,0].set_title('$\epsilon$ - $F_{m}$')


# Set up the tick labels
plt.setp([a.get_xticklabels() for a in axi[0, 1:]], visible=False)
plt.setp([a.get_xticklabels() for a in axi[1, 1:]], visible=False)
plt.setp([a.get_xticklabels() for a in axi[2, 1:]], visible=False)
plt.setp([a.get_yticklabels() for a in axi[:, 1]], visible=False)
plt.setp([a.get_yticklabels() for a in axi[:, 2]], visible=False)
plt.setp([a.get_yticklabels() for a in axi[:, 3]], visible=False)
plt.setp([a.get_yticklabels() for a in axi[:, 4]], visible=False)
plt.setp([a.get_yticklabels() for a in axi[:, 5]], visible=False)
plt.setp([a.get_yticklabels() for a in axi[:, 6]], visible=False)

# Set up the tick sizes
plt.setp([a.get_yticklabels() for a in axi[0, :]], fontsize=5)
plt.setp([a.get_yticklabels() for a in axi[1, :]], fontsize=5)
plt.setp([a.get_yticklabels() for a in axi[2, :]], fontsize=5)
plt.setp([a.get_xticklabels() for a in axi[1, :]], fontsize=5)
plt.setp(axi[0,0].get_xticklabels(),fontsize=5)
plt.setp(axi[1,0].get_xticklabels(),fontsize=5)
plt.setp(axi[2,0].get_xticklabels(),fontsize=5)

axi[0,0].xaxis.set_ticks([0.4,0.9])
axi[1,0].xaxis.set_ticks([0.4,0.9])
axi[2,0].xaxis.set_ticks([0.05,5,10])

c = 'stress'

if c == 'max':
    
    row_cnt = 0
    col_cnt = -1
    sfac = 1000
    vmin_val = 5
    vmax_val = 55
    
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
        
        line1 = axi[row_cnt, col_cnt].scatter(x,y,c=z,s=sfac*z,vmin=vmin_val,vmax=vmax_val)
        
        if col_cnt == 6:
            cax = plt.axes([0.93, 0.1, 0.02, 0.8])
            plt.colorbar(line1,cax=cax)
            cax.tick_params(width=1.1,labelsize=5)
            bbox_props = dict(boxstyle="rarrow,pad=0.01", fc="cyan", ec="b", lw=0.5)
            t = axi[row_cnt, col_cnt].text(-2, 4.3, "$+ \epsilon$", ha="center", va="center", rotation=0, \
            bbox=bbox_props, fontsize=8)
            
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
        
        axi[row_cnt, col_cnt].scatter(x,y,c=z,s=sfac*z,vmin=vmin_val,vmax=vmax_val) 
        
        if col_cnt == 6:
            bbox_props = dict(boxstyle="rarrow,pad=0.01", fc="cyan", ec="b", lw=0.5)
            t = axi[row_cnt, col_cnt].text(-2, 14, "$+ F_m$", ha="center", va="center", rotation=0, \
            bbox=bbox_props, fontsize=6)
        
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
        
        axi[row_cnt, col_cnt].scatter(x,y,c=z,s=sfac*z,vmin=vmin_val,vmax=vmax_val)
        
        if col_cnt == 5:
            bbox_props = dict(boxstyle="rarrow,pad=0.01", fc="cyan", ec="b", lw=0.5)
            t = axi[row_cnt, col_cnt].text(-26, 4, "$+ \phi$", ha="center", va="center", rotation=0, \
            bbox=bbox_props, fontsize=6)
            
        if col_cnt == 2:
            axi[row_cnt, col_cnt].set_title("x = $\epsilon$, y = $F_m$",fontsize=8)
       
    
        
    plt.savefig('phase_plot_max.png',bbox_inches='tight')
    
elif c == 'bond':
    
    row_cnt = 0
    col_cnt = -1
    sfac = 100000
    vmin_val = 0.01
    vmax_val = 0.2
    
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
        
        line1 = axi[row_cnt, col_cnt].scatter(x,y,c=z,vmin=vmin_val,vmax=vmax_val)
        
        if col_cnt == 6:
            cax = plt.axes([0.93, 0.1, 0.02, 0.8])
            plt.colorbar(line1,cax=cax)
            cax.tick_params(width=1.1,labelsize=5)
            bbox_props = dict(boxstyle="rarrow,pad=0.01", fc="cyan", ec="b", lw=0.5)
            t = axi[row_cnt, col_cnt].text(-2, 4.3, "$+ \epsilon$", ha="center", va="center", rotation=0, \
            bbox=bbox_props, fontsize=8)
            
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
        
        axi[row_cnt, col_cnt].scatter(x,y,c=z,vmin=vmin_val,vmax=vmax_val) 
        
        if col_cnt == 6:
            bbox_props = dict(boxstyle="rarrow,pad=0.01", fc="cyan", ec="b", lw=0.5)
            t = axi[row_cnt, col_cnt].text(-2, 14, "$+ F_m$", ha="center", va="center", rotation=0, \
            bbox=bbox_props, fontsize=6)
        
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
        
        axi[row_cnt, col_cnt].scatter(x,y,c=z,vmin=vmin_val,vmax=vmax_val)
        
        if col_cnt == 5:
            bbox_props = dict(boxstyle="rarrow,pad=0.01", fc="cyan", ec="b", lw=0.5)
            t = axi[row_cnt, col_cnt].text(-26, 4, "$+ \phi$", ha="center", va="center", rotation=0, \
            bbox=bbox_props, fontsize=6)
            
        if col_cnt == 2:
            axi[row_cnt, col_cnt].set_title("x = $\epsilon$, y = $F_m$",fontsize=8)
       
    
        
    plt.savefig('phase_plot_bond_2.eps',bbox_inches='tight')
    
    
elif c == 'vel':

    row_cnt = 0
    col_cnt = -1
    sfac = 1000000
    vmin_val = 0.001
    vmax_val = 0.23
 
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
        
        line1 = axi[row_cnt, col_cnt].scatter(x,y,c=z,vmin=vmin_val,vmax=vmax_val)
        
        if col_cnt == 6:
            cax = plt.axes([0.93, 0.1, 0.02, 0.8])
            plt.colorbar(line1,cax=cax)
            cax.tick_params(width=1.1,labelsize=5)
            bbox_props = dict(boxstyle="rarrow,pad=0.01", fc="cyan", ec="b", lw=0.5)
            t = axi[row_cnt, col_cnt].text(-2, 4.3, "$+ \epsilon$", ha="center", va="center", rotation=0, \
            bbox=bbox_props, fontsize=8)
            
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
        
        axi[row_cnt, col_cnt].scatter(x,y,c=z,vmin=vmin_val,vmax=vmax_val) 
        
        if col_cnt == 6:
            bbox_props = dict(boxstyle="rarrow,pad=0.01", fc="cyan", ec="b", lw=0.5)
            t = axi[row_cnt, col_cnt].text(-2, 14, "$+ F_m$", ha="center", va="center", rotation=0, \
            bbox=bbox_props, fontsize=6)
        
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
        
        axi[row_cnt, col_cnt].scatter(x,y,c=z,vmin=vmin_val,vmax=vmax_val)
        
        if col_cnt == 5:
            bbox_props = dict(boxstyle="rarrow,pad=0.01", fc="cyan", ec="b", lw=0.5)
            t = axi[row_cnt, col_cnt].text(-26, 4, "$+ \phi$", ha="center", va="center", rotation=0, \
            bbox=bbox_props, fontsize=6)
            
        if col_cnt == 2:
            axi[row_cnt, col_cnt].set_title("x = $\epsilon$, y = $F_m$",fontsize=8)
       
    
        
    plt.savefig('phase_plot_vel.eps',bbox_inches='tight')

elif c == 'stress':
    
    row_cnt = 0
    col_cnt = -1
    sfac = 100000
    #vmin_val = -5e-5
    #vmax_val = 1e-7
    vmin_val = -0.5
    vmax_val = 0.7
    
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
        
        line1 = axi[row_cnt, col_cnt].scatter(x,y,c=z,vmin=vmin_val,vmax=vmax_val)
        
        if col_cnt == 6:
            cax = plt.axes([0.93, 0.1, 0.02, 0.8])
            plt.colorbar(line1,cax=cax)
            cax.tick_params(width=1.1,labelsize=5)
            bbox_props = dict(boxstyle="rarrow,pad=0.01", fc="cyan", ec="b", lw=0.5)
            t = axi[row_cnt, col_cnt].text(-2, 4.3, "$+ \epsilon$", ha="center", va="center", rotation=0, \
            bbox=bbox_props, fontsize=8)
            
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
        
        axi[row_cnt, col_cnt].scatter(x,y,c=z,vmin=vmin_val,vmax=vmax_val) 
        
        if col_cnt == 6:
            bbox_props = dict(boxstyle="rarrow,pad=0.01", fc="cyan", ec="b", lw=0.5)
            t = axi[row_cnt, col_cnt].text(-2, 14, "$+ F_m$", ha="center", va="center", rotation=0, \
            bbox=bbox_props, fontsize=6)
        
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
        
        axi[row_cnt, col_cnt].scatter(x,y,c=z,vmin=vmin_val,vmax=vmax_val)
        
        if col_cnt == 5:
            bbox_props = dict(boxstyle="rarrow,pad=0.01", fc="cyan", ec="b", lw=0.5)
            t = axi[row_cnt, col_cnt].text(-26, 4, "$+ \phi$", ha="center", va="center", rotation=0, \
            bbox=bbox_props, fontsize=6)
            
        if col_cnt == 2:
            axi[row_cnt, col_cnt].set_title("x = $\epsilon$, y = $F_m$",fontsize=8)
       
    
        
    plt.savefig('phase_plot_stress_int.eps',bbox_inches='tight')
