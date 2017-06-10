
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
from matplotlib.lines import Line2D
import seaborn as sns
sns.set(style="white",context='paper',
        font_scale=1.2,font="Open Sans",
        rc={'mathtext.default': 'regular','font.size': 30, 
            'font.family': 'sans',"figure.dpi":300,
            "xtick.major.size": 8, "ytick.major.size": 8,
            'grid.linestyle': '--'})   
#sns.set()
#sns.set(style="white",context='notebook',
#        font_scale=1.4,font="Open Sans",
#        rc={'mathtext.default': 'regular','font.size': 20, 
#            'font.family': 'sans',"figure.dpi":300,"axes.formatter.limits":(-2,2)})

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

    return

##############################################################################
    
def plot_data(data, param_choice, args):
    """ plot the frames per parameter"""
     
    ### set general plot properties
    
    savebase = '/usr/users/iff_th2/duman/Cells_in_LAMMPS/POVRAY/'
    #downlim = -1
    #uplim = sim.lx/4.
    num_ticks = 5
    ax_len = 1.0                          # Length of one subplot square box
    ax_b = 0.0                            # Beginning/offset of the subplot in the box
    ax_sep = 0.0                          # Separation length between two subplots
    total_subplots_in_x = 1               # Total number of subplots    
    fig = plt.figure()
    subp = misc_tools.Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x) 
    ax0 = subp.addSubplot()
 
    name = ''
    pname = ''
    if param_choice == 'areak': 
        name = 'AREAK'
        pname = name + '_eps_' + str(args.eps) + '_fp_' + str(args.fp) + \
            '_kappa_' + str(args.kappa)
        xlab = '$\kappa_A$'
        tit = '$\epsilon=$' + str(args.eps) + ',$f_m=$' + str(args.fp) + \
            ',$\kappa=$' + str(args.kappa)
    elif param_choice == 'eps':
        name = 'EPS'
        pname = name + '_fp_' + str(args.fp) + '_areak_' + str(args.areak) + \
            '_kappa_' + str(args.kappa)
        xlab = '$\epsilon$'
        tit = '$f_m=$' + str(args.fp) + ',$\kappa_A=$' + str(args.areak) + \
            ',$\kappa=$' + str(args.kappa)        
    elif param_choice == 'fp':
        name = 'FP'
        pname = name + '_eps_' + str(args.eps) + '_areak_' + str(args.areak) + \
            '_kappa_' + str(args.kappa)
        xlab = '$f_{m}$'
        tit = '$\epsilon=$' + str(args.eps) + ',$\kappa_A=$' + str(args.areak) + \
            ',$\kappa=$' + str(args.kappa)          
    elif param_choice == 'kappa':
        name = 'KAPPA'
        pname = name + '_eps_' + str(args.eps) + '_fp_' + str(args.fp) + \
            '_areak_' + str(args.areak)
        xlab = '$\kappa$'
        tit = '$\epsilon=$' + str(args.eps) + ',$f_m=$' + str(args.fp) + \
            ',$\kappa_A=$' + str(args.areak) 
    base = savebase + name + '/'
    os.system("mkdir -p " + base)  
                         
    ### plot 

    subp = misc_tools.Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x) 
    ax0 = subp.addSubplot()
    
    x = data.keys()
    y = [1 for j in range(len(data.keys()))]
    print x
    ax0.scatter(x, y)
    ax0.set_xscale('log')
    ax0.set_yscale('log')
    
    for j, p in enumerate(data.keys()):
        
        fname = data[p]         
        
        if os.path.exists(fname):
            arr_hand = read_png(fname)
        
            zoom=0.099
            imagebox = OffsetImage(arr_hand, zoom=zoom)

            xy = [x[j], y[j]]               # coordinates to position this image

            ab = AnnotationBbox(imagebox, xy,
                xybox=(0., -0.),
                xycoords='data',
                boxcoords="offset points",frameon=1,pad=.1)  
                                
            ax0.add_artist(ab)
            
    ### title
    
    ax0.set_title(tit, fontsize=30)
    
    ### labels
        
    ax0.set_xlabel(xlab, fontsize=30)
    #ax0.set_ylabel("$F_{s}(q,\\Delta t)$", fontsize=40)

    ### limits

    #ax0.set_xlim((-1, 15))
    ax0.set_ylim((0.9999, 1.0001))
    
    ax0.grid(1, color='#cccccc', linestyle='--')
    ax0.set_frame_on(False)
    ax0.get_xaxis().tick_bottom()
    ax0.axes.get_yaxis().set_visible(False)
    xmin, xmax = ax0.get_xaxis().get_view_interval()
    ymin, ymax = ax0.get_yaxis().get_view_interval()
    ax0.add_artist(Line2D((xmin, xmax), (ymin, ymin), color='black', linewidth=2))       
    ### ticks
    
    #ax0.xaxis.set_ticks(np.linspace(0, 15, num_ticks, endpoint=True))
    #ax0.yaxis.set_ticks(np.linspace(0, uplim, num_ticks, endpoint=True))
    plt.setp(ax0.get_yticklabels(),visible=False)   
    ax0.tick_params(axis='both', which='major', labelsize=30)
    
    ### legend

#    ax0.legend(bbox_to_anchor=(1.005, 0.,0.65, 1.), loc=2, borderaxespad=0., \
#        prop={'size': 20}, mode="expand", frameon=False)
    
    ### save  
    
    savepath = base + "images_per_" + pname + ".pdf"
    print savepath
    plt.savefig(savepath, dpi=300, bbox_inches='tight', pad_inches=0.08)        
    fig.clf()                             
    
    return
    
##############################################################################
    
def gen_folder(base, e, f, a, k):
    """ index the folder addresses based on the parameters and the paths"""
                    
    return base + 'eps_' + str(e) + '_fp_' + str(f) + \
                   '_areak_' + str(a) + '_kappa_' + str(k)
                
##############################################################################

def main():

    ### set the premilinary info
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--eps", type=float, nargs="?", const=-1, \
                        help="Strength of LJ potential")
    parser.add_argument("-f", "--fp", type=float, nargs="?", const=-1, \
                        help="Propulsion force")
    parser.add_argument("-a", "--areak", type=float, nargs="?", const=-1, \
                        help="Area constraint potential strength")
    parser.add_argument("-k", "--kappa", type=float, nargs="?", const=-1, \
                        help="Bending rigidity")     
    args = parser.parse_args()

    ### detect the parameter choice 
    
    param = []
    param_choice = ''
    if args.eps == -1:
        param_choice = 'eps'
        param = [0.05, 1.0, 5.0, 20.0]
    if args.fp == -1:
        param_choice = 'fp'
        param = [0.5, 1.0, 5.0]
    if args.areak == -1:
        param_choice = 'areak'
        param = [1.0, 10.0, 100.0]
    if args.kappa == -1:
        param_choice = 'kappa'
        param = [1.0, 10.0, 100.0, 1000.0]

    ### get the data per detected parameter

    data = {}       # carries the data per parameter set
    base = '/usr/users/iff_th2/duman/Cells_in_LAMMPS/POVRAY/'
    
    for p in param:
        
        if param_choice == 'areak':
            datafolder = gen_folder(base, args.eps, args.fp, p, args.kappa)
        elif param_choice == 'eps':
            datafolder = gen_folder(base, p, args.fp, args.areak, args.kappa)
        elif param_choice == 'fp':            
            datafolder = gen_folder(base, args.eps, p, args.areak, args.kappa)
        elif param_choice == 'kappa':            
            datafolder = gen_folder(base, args.eps, args.fp, args.areak, p)
            
        data[p] = datafolder + '/pov-frame-00800.png'
                
    plot_data(data, param_choice, args)            
    
#    ### load the pictures
#    
#    files = []
#    tmp = {}
#    for folder in folders:
#        
#        ### get the file path address which encapsulates all the parameters
#        
#        filekeys = folder.split('/')[-1]
#        print filekeys.split('_')
#        
#        ### add the parameters and the files to the dictionary
#        
#        tmp = {}
#        ename = filekeys.split('_')[0]
#        evalue = atof(filekeys.split('_')[1])
#        edict = {}
#        edict[ename] = evalue
#        
#        fname = filekeys.split('_')[2]
#        fvalue = atof(filekeys.split('_')[3])
#        fdict = {}
#        fdict[fname] = fvalue
#        
#        aname = filekeys.split('_')[4]
#        avalue = atof(filekeys.split('_')[5])
#        adict = {}
#        adict[aname] = avalue
#        
#        tmp.update(edict)
#        tmp.update(fdict)
#        tmp.update(adict)
#        
#        tmp.update({'folder':folder+'/pov-frame-00999.png', 'eps':tmp['eps'], 'fp':tmp['fp'], 'areak':tmp['areak']})
#        files.append(tmp)
#    
#    Data = pd.DataFrame(files, columns=files[0].keys())
#    sliceby = 'eps'
#    bins = np.linspace(Data[sliceby].min(), Data[sliceby].max(), 5)
#    
#    plot1={'sliceby':'fp',
#    'namex':'eps',
#    'namey':'areak'}
#    
#    plot_slices(Data, save=True, **plot1)
    
##############################################################################
    
if __name__ == '__main__':
    main()
    
##############################################################################
    
