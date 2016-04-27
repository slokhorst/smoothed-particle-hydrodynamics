import os
import glob
import subprocess

num_frames = 100
res = (1280, 720)
dlen = len(str(num_frames))
subprocess.run(['povray', 'Output_File_Name=data/img.png', '+W{}'.format(res[0]), '+H{}'.format(res[1]), '+KFF{}'.format(num_frames), 'sea.pov'])
subprocess.run(['ffmpeg', '-y', '-framerate', '10', '-i', 'data/img%0{}d.png'.format(dlen), '-c:v', 'libx264', 'video.mp4'])