
""" Build a picture based phase diagram"""

### example command line arguments: 
###    -fl=/local/duman/SIMULATIONS/Cells_in_LAMMPS/.../ 
###         -sb=/usr/users/iff_th2/duman/Cells_in_LAMMPS/MOVIES/
###             -ti=10 -tf=2000

##############################################################################

import sys
sys.path.append('../Utility')

import argparse
import numpy as np
import os
import matplotlib as mpl
mpl.use('Agg', warn=False)
import matplotlib.pyplot as plt
import read_write
import misc_tools
import glob
import pandas as pd
from string import atof
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from matplotlib._png import read_png
import colorsys 
import seaborn as sns
sns.set()
sns.set(style="white",context='notebook',
        font_scale=1.4,font="Open Sans",
        rc={'mathtext.default': 'regular','font.size': 20, 
            'font.family': 'sans',"figure.dpi":120,"axes.formatter.limits":(-2,2)})

##############################################################################

def plot_slices(Data, sliceby, namex, namey, save=False, scale=.15, savesuffix='', binning=False):
    
    labels={'areak':r"$\kappa_{A}$",
           'fp':r"$f_{m}$",
           'eps':'$\epsilon$'}
    
    if isinstance(binning, int) and binning:
        bins = np.linspace(Data[sliceby].min(), Data[sliceby].max(), binning)
        Data=Data.copy()
        Data['digi']=Data[sliceby].apply(lambda x:np.digitize(x,bins))
        group = Data.groupby('digi')
    else:
        group=Data.groupby(sliceby)

    for i, grp in group:
        
        grp=grp.copy()
        fwidth=10.
        fig,ax=plt.subplots(1,1,figsize=(fwidth,fwidth*3/4.))
        ax.set_title(r'{1:s}: {0:.2f}'.format(grp[sliceby].mean(),labels[sliceby]))

        namex2x=dict([(j,i) for i,j in enumerate(sorted(grp[namex].unique()))])
        namey2y=dict([(j,i) for i,j in enumerate(sorted(grp[namey].unique()))])
        grp['x']=grp[namex].apply(lambda x:namex2x[x])
        grp['y']=grp[namey].apply(lambda x:namey2y[x])


        ax.scatter(grp.x,grp.y)
        ax.set_xlabel(labels[namex])
        ax.set_ylabel(labels[namey])

        show_screenshot=True
        if show_screenshot:
            for fname,x,y in zip(grp.folder,grp.x,grp.y):

                if os.path.exists(fname):
                    arr_hand = read_png(fname)
                
                    zoom=scale
                    if arr_hand.shape[0]==1024:
                        zoom/=2.
                    imagebox = OffsetImage(arr_hand, zoom=zoom)
    
                    xy = [x,y]               # coordinates to position this image
    
                    ab = AnnotationBbox(imagebox, xy,
                        xybox=(0., -0.),
                        xycoords='data',
                        boxcoords="offset points",frameon=1,pad=.1)                                  
                    ax.add_artist(ab)
    
            ax.xaxis.set_ticks(sorted(grp.x.unique()))
            ax.yaxis.set_ticks(sorted(grp.y.unique()))
    
            xticks=sorted(grp[namex].unique().tolist())
            yticks=sorted(grp[namey].unique().tolist())
            xlabels=["{0:.1f}".format(j) for j in xticks]
            ylabels=["{0:.2f}".format(j) for j in yticks]
            ax.xaxis.set_ticklabels(xlabels)
            ax.yaxis.set_ticklabels(ylabels)
    
            ax.grid(1,color='#cccccc',linestyle='--')
            if save:
                fig.savefig("/usr/users/iff_th2/duman/Cells_in_LAMMPS/POVRAY//set_{1:s}{0:.2f}.pdf".format(grp[sliceby].mean(),sliceby+savesuffix),bbox_inches='tight',dpi=300)    

##############################################################################

def main():

    ### index the data

    eps = [0.05, 0.5, 1.0, 5.0, 10.0, 20.0]
    fp = [0.0, 0.5, 1.0, 3.0, 5.0, 10.0]
    areak = [1.0, 10.0, 100.0]
    
    base = '/usr/users/iff_th2/duman/Cells_in_LAMMPS/POVRAY/'
    folders = []
    for e in eps:
        for f in fp:
            for a in areak:
                folders.append( base + 'eps_' + str(e) + '_fp_' + str(f) + '_areak_' + str(a) )
                
            
    database = '/local/duman/SIMULATIONS/Cells_in_LAMMPS/density_0.8/'
    
    ### load the pictures
    
    files = []
    tmp = {}
    for folder in folders:
        
        ### get the file path address which encapsulates all the parameters
        
        filekeys = folder.split('/')[-1]
        print filekeys.split('_')
        
        ### add the parameters and the files to the dictionary
        
        tmp = {}
        ename = filekeys.split('_')[0]
        evalue = atof(filekeys.split('_')[1])
        edict = {}
        edict[ename] = evalue
        
        fname = filekeys.split('_')[2]
        fvalue = atof(filekeys.split('_')[3])
        fdict = {}
        fdict[fname] = fvalue
        
        aname = filekeys.split('_')[4]
        avalue = atof(filekeys.split('_')[5])
        adict = {}
        adict[aname] = avalue
        
        tmp.update(edict)
        tmp.update(fdict)
        tmp.update(adict)
        
        tmp.update({'folder':folder+'/pov-frame-00999.png', 'eps':tmp['eps'], 'fp':tmp['fp'], 'areak':tmp['areak']})
        files.append(tmp)
    
    Data = pd.DataFrame(files, columns=files[0].keys())
    sliceby = 'eps'
    bins = np.linspace(Data[sliceby].min(), Data[sliceby].max(), 5)
    
    plot1={'sliceby':'fp',
    'namex':'eps',
    'namey':'areak'}
    
    plot_slices(Data, save=True, **plot1)
    
##############################################################################
    
if __name__ == '__main__':
    main()
    
##############################################################################
    
