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

ft = FortranFile('stream1','r')
x = ft.read_reals()
ft.close()
#print (x[0],x[1],x[100], np.size(x))
nxx = 257
nyy = 129
nl = 2
H1 = 5000.
H2 = 5000.
sbt = np.zeros([nxx,nyy])
sbc = np.zeros([nxx,nyy])
s1 = np.reshape(x[0:nxx*nyy],(nyy,nxx))

plt.contourf(s1,20,cmap=cm.gist_ncar)
plt.colorbar()
plt.savefig('1.png')


ft = FortranFile('stream00100','r')
x = ft.read_reals()
ft.close()
#print (x[0],x[1],x[100], np.size(x))
nxx = 257
nyy = 129
nl = 2
H1 = 5000.
H2 = 5000.
sbt = np.zeros([nxx,nyy])
sbc = np.zeros([nxx,nyy])
s1 = np.reshape(x[0:nxx*nyy],(nyy,nxx))

plt.contourf(s1,20,cmap=cm.gist_ncar)
plt.savefig('2.png')

# Build GIF
with imageio.get_writer('mygif.gif', mode='I') as writer:
    for filename in ['1.png', '2.png']:
        image = imageio.imread(filename)
        writer.append_data(image)
        
#plt.show()
