#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 23:05:33 2017

@author: peterjanvanleeuwen
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import imageio
from matplotlib import cm
from scipy.io import FortranFile

plt.close('all')

dx = 100000.
nxx = 257
nyy = 129
nl = 2
H1 = 5000.
H2 = 5000.

filenames = []
images = []

#Read file at time step t=1
ft = FortranFile('stream1','r')
x = ft.read_reals()
ft.close()


sbt = np.zeros([nxx,nyy])
sbc = np.zeros([nxx,nyy])
s1 = np.reshape(x[0:nxx*nyy],(nyy,nxx))
s2 = np.reshape(x[nxx*nyy:2*nxx*nyy],(nyy,nxx))


#Plot top layer (s1) or bottom layer (s2)
####plt.contourf(s1,20,cmap=cm.gist_ncar)
plt.contourf(s2,20,cmap=cm.gist_ncar)
plt.colorbar()

plt.title('Stream Function Field at t=1', fontsize=18)
plt.xlabel('Zonal', fontsize=16)
plt.ylabel('Meridional', fontsize=16)

plt.savefig('1.png')

images.append('1.png')


for i in range(10,2000,10):
  if i < 99.:
    ft = FortranFile('stream000'+str(i),'r')
  elif i < 999.:
    ft = FortranFile('stream00'+str(i),'r')
  else:
    ft = FortranFile('stream0'+str(i),'r')  
  x = ft.read_reals()
  ft.close()

  sbt = np.zeros([nxx,nyy])
  sbc = np.zeros([nxx,nyy])
  s1 = np.reshape(x[0:nxx*nyy],(nyy,nxx))
  s2 = np.reshape(x[nxx*nyy:2*nxx*nyy],(nyy,nxx))
  
  #Plot top layer (s1) or bottom layer (s2)
  ####plt.contourf(s1,20,cmap=cm.gist_ncar)
  plt.contourf(s2,20,cmap=cm.gist_ncar)
  plt.title('Stream Function Field at t='+str(i), fontsize=18)
  plt.xlabel('Zonal', fontsize=16)
  plt.ylabel('Meridional', fontsize=16)

  plt.savefig(str(i)+'.png')

  images.append(str(i)+'.png') 
  
#images=images.append(['1.png'])
#for j in range(100,102):
#  images = images.append([str(j)+'.png']) 
  
  
# Build animated GIF
with imageio.get_writer('stream_s2_animated.gif', mode='I') as writer:
    #for filename in ['1.png', '100.png', '101.png']:
    #for images in [str(j)+'.png']:
    for filename in images:
        image = imageio.imread(filename)
        writer.append_data(image)
        
#plt.show()
