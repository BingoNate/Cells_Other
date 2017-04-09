
""" Generate a movie of the beads of cells in povray format"""

### example command line arguments: 
###    -fl=/local/duman/SIMULATIONS/Cells_in_LAMMPS/.../ 
###         -sb=/usr/users/iff_th2/duman/Cells_in_LAMMPS/POVRAY/
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
import vapory
import matplotlib.cm as cm
import matplotlib.colors as mplcolors 
 
##############################################################################
    
def gen_img_settings_quality(l):
    """ generate the general povray settings"""
    
    lhalf = 0.5*l
    
    ### sphere radius
    
    sphere_radius = 0.7
    #sphere_rgbcolor = [0.25,0.65,0.65]
    
    ### RESOLUTION
    
    img_widthpx = 1024
    img_heightpx = 1024

    ### includes and defaults

    povray_includes = ["colors.inc", "textures.inc", "shapes.inc"]
    povray_defaults = [vapory.Finish( 'ambient', 0.1,
	     			  'diffuse', 0.65,
		    		  'specular', 0.5,
			    	  'shininess', 0.53,
				  'opacity', 1.0)]


    ### light sources

    sun1 = vapory.LightSource([lhalf, lhalf, -1.01*lhalf], 'color', 'White')
    sun2 = vapory.LightSource([lhalf, lhalf, -1.01*lhalf], 'color', [0.7, 0.7, 0.7])

    ### background

    background = vapory.Background('color', [1,1,1])

    ### camera

    #povray_cam = vapory.Camera('angle', 75, 'location',  [-15 , 15.0+0.5,15.0-0.25],'look_at', [0.25 , 15.0+0.5, 15.0-0.25])
    povray_cam = vapory.Camera('location', [lhalf, lhalf, -1.01*lhalf], 'look_at', [lhalf,lhalf,0], 'angle', 90)

    ### text
    # If desired include this in the povray_objects - array declared in the loop
    #text1 = vapory.Text( 'ttf', '"timrom.ttf"' ,'"Division:"', 0.01, 0.0, 'scale', [0.5,0.5,0.5],'rotate', [0,90,0], 'translate' , [0.0 , 15.0+2.75-1 , 15.0+1.5], vapory.Pigment('Black') ) 

    ### render quality

    quality = 10
    
    return sphere_radius, img_widthpx, img_heightpx, povray_includes, povray_defaults, sun1, sun2, background, povray_cam, quality

##############################################################################

def gen_colors(nbeads):
    """ generate an array with colors for the spheres"""

    # use colormap (blue, red, green)
    #  blue 0.0, 0.0, 1.0
    #  red: 1.0, 0.0, 0.0
    #  green: 0.0, 1.0, 1.0

    my_cmap = cm.get_cmap('gnuplot')
    minval = 0
    maxval = nbeads - 1
    norm = mplcolors.Normalize(minval, maxval)
    
    sphere_rgbcolor = []
    for i in range(nbeads):
        si = my_cmap(norm(i))
        sphere_rgbcolor.append([si[0], si[1], si[2]])

    return sphere_rgbcolor
       
##############################################################################

def plot_frames(beads, sim, ti, tf, savebase):
    """ plot frames within the specified time window"""
    
    ### define the color for the spheres

    print 'defining colors'
    sphere_rgbcolor = gen_colors(sim.nbeads)

    ### create povray settings

    print 'creating povray settings'
    sphere_radius, img_widthpx, img_heightpx, povray_includes, \
        povray_defaults, sun1, sun2, background, povray_cam, quality \
            = gen_img_settings_quality(sim.lx)
    
    zi = np.zeros((sim.nbeads))
        
    ### set general plot properties

    savebase += 'eps_' + str(sim.eps) + '_fp_' + str(sim.fp) + '_areak_' + str(sim.areak) + '/'
    os.system("mkdir -p " + savebase)
    
    ### plot the frames
    
    for step in range(ti, tf):
        
        time = step*sim.dt
        print 'Step / Total : ', step, tf
        
        ### create povray items
        
        print 'generating povray item'
        particles = vapory.Object( \
            vapory.Union( \
                *[ vapory.Sphere([beads.xi[step, 0, j], beads.xi[step, 1, j],zi[j]], \
                    sphere_radius, vapory.Texture( \
                        vapory.Pigment('color', sphere_rgbcolor[j]), \
                            vapory.Finish('phong',1)) ) for j in range(0, sim.nbeads ) ] ) )

        ### generate povray objects

        print 'generating povray objects'
        povray_objects = [sun1, sun2, background, particles]
        ### create the scene
        scene = vapory.Scene( camera = povray_cam,
                       objects = povray_objects, 
                       included = povray_includes, 
                       defaults = povray_defaults )
                       
        ### render image
                           
        print 'rendering scene'
        savename = "pov-frame-" + "{0:05d}".format(int(step)) + ".png"
        scene.render(outfile=savename, width=img_widthpx, height=img_heightpx, \
            antialiasing=0.001, quality=quality, remove_temp=True)
            
        ### move the image to the correct destination
            
        os.system('mv ' + savename + ' ' + savebase)
        
    return
        
##############################################################################

def main():

    ### get the data folder
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-fl", "--folder", \
                        help="Folder containing data, as in /local/duman/SIMULATIONS/Cells_in_LAMMPS/.../")
    parser.add_argument("-sb", "--savebase", nargs="?", \
                        const = "/usr/users/iff_th2/duman/Cells_in_LAMMPS/MOVIES/", \
                        help="Folder to save the data, as in /usr/users/iff_th2/duman/Cells_in_LAMMPS/POVRAY/")     
    parser.add_argument("-ti","--init_time", nargs="?", const=999, type=int, \
                        help="First frame of the video (in terms of frame number), you can also leave it empty")
    parser.add_argument("-tf","--fin_time", nargs="?", const=1000, type=int, \
                        help="Last frame of the video (in terms of frame number), you can also leave it empty")
    args = parser.parse_args()
    
    ### read the data and general information from the folder
    
    sim, cells, beads = read_write.read_h5_file(args.folder)
    beads.get_img_pos(sim.lx)
    print "folder = ", args.folder
        
    ### plot the data in the given time window
    
    plot_frames(beads, sim, args.init_time, args.fin_time, args.savebase)
    
    return
    
##############################################################################

if __name__ == '__main__':
    main()    
    
##############################################################################
