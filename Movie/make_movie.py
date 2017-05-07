
""" Generate a movie of the beads of cells in povray format"""

### example command line arguments: 
###    -fl=/usr/users/iff_th2/duman/Cells_in_LAMMPS/MOVIES/eps_0.05_fp_1.0_areak_1.0/

##############################################################################

import subprocess
import os
import argparse

##############################################################################

def order_files(folder):
    """ order the filenames for ffmpeg to work"""

    ### change/order the filenames for ffmpeg to work
    
    final_vid = 1
    enum = []
    
    if 1 == final_vid:
    	
    	for ff in os.listdir(folder):	
    	
    		if ".png" == ff[-4:] and "frame" == ff[0:5]:
    			if "." == ff[11]:
    				num = int(ff[6:11])
    			elif "." == ff[12]:
    				num = int(ff[6:12])
    			elif "." == ff[13]:
    				num = int(ff[6:13])
    			elif "." == ff[14]:
    				num = int(ff[6:14])
    			enum.append([num,ff])
    
    		enum = sorted(enum,key=lambda x:x[0])
    
    	for i in range(0,len(enum)):
    		print enum[i][1]
    		os.rename(enum[i][1],"frame-"+ "%05d" % i +".png")
      
    return
      
##############################################################################

def gen_video(folder, savefolder):
    """ generate the video with ffmpeg""" 	
    
    ### generate the video 
    
    os.system('mkdir -p ' + savefolder)
    v_name = savefolder + 'detailed.mp4'

    # -r : framerate(fps), -s is resoluation
    # -i is input %04 is to pad with zeros until 4th string pic0001, pic0002, pic0020, ...
    # -crf is quality, the lower the better
    # pix_fmt is the pixel format
    # -y is to force overwrite

#    subprocess.call(['ffmpeg','-r','40','-f','image2','-s','720:720','-i',folder+'frame-%05d.png','-y',
#                     '-vcodec','libx264','-crf','25','-pix_fmt','yuv420p','-vf','scale=720:-1',v_name])
    subprocess.call(['ffmpeg','-f','image2','-s','720:720','-i',folder+'frame-%05d.png','-y',
                     '-vcodec','libx264','-crf','25','-pix_fmt','yuv420p','-vf','scale=720:-1',v_name])  
#    subprocess.call(['ffmpeg','-f','image2','-s','720:720','-i',folder+'frame-%05d.png','-y',
#                     '-vcodec','libx264','-crf','25','-pix_fmt','yuv420p','-vf','scale=720:trunc(ow/a/2)*2',v_name])  

    return

##############################################################################

def main():

    ### get the data folder
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-fl", "--folder", \
                        help="Folder containing images, as in /local/duman/Cells_in_LAMMPS/MOVIES/eps_0.05_fp_1.0_areak_1.0/")
    parser.add_argument("-sfl", "--savefolder", \
                        help="Folder containing images, as in /usr/users/iff_th2/duman/Cells_in_LAMMPS/MOVIES/eps_0.05_fp_1.0_areak_1.0/")                        
    args = parser.parse_args()
    
    ### generate the video in the given folder

    order_files(args.folder)
    gen_video(args.folder, args.savefolder)
    
    return
    
##############################################################################

if __name__ == '__main__':
    main()

##############################################################################
    
