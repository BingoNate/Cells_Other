#!/bin/python
 
#Copyleft Arvind Ravichandran
#Tue Feb 10 11:21:45 CET 2015
#st2_video.py
#Description:
 
import subprocess
import os

final_vid = 1;
enum = [];

if 1 == final_vid:
	
	for ff in os.listdir("."):	
	
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

# print 'rendering video'
v_name = 'video-3'

subprocess.call(['ffmpeg','-i','frame-%05d.png',
'-vcodec','mpeg2video','-qscale','0','-filter:v','setpts=2.0*PTS','-vf','scale=640:480','./'+v_name+'.mpeg'])


