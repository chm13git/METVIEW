#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 23:05:33 2017

@author: peterjanvanleeuwen
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.io import FortranFile

plt.close('all')

dx = 100000.

ft = FortranFile('stream1','r')
x = ft.read_reals()
ft.close()
print (x[0],x[1],x[100], np.size(x))
nxx = 257
nyy = 129
nl = 2
H1 = 5000.
H2 = 5000.
sbt = np.zeros([nxx,nyy])
sbc = np.zeros([nxx,nyy])
s1 = np.reshape(x[0:nxx*nyy],(nyy,nxx))
s2 = np.reshape(x[nxx*nyy:2*nxx*nyy],(nyy,nxx))
sbt = (H1*s1 + H2*s2)/(H1+H2)
sbc = s1 - s2
plt.figure(2)
plt.contourf(s1,20,cmap=cm.gist_ncar)
plt.colorbar()
plt.figure(3)
plt.contourf(s2,20,cmap=cm.gist_ncar)
plt.colorbar()
plt.figure(6)
plt.plot(s1[:,0],'g')
print(np.sum(s1[nyy-1,:]-s1[nyy-2,:]),np.sum(s2[nyy-1,:]-s2[nyy-2,:]))
print(np.sum(sbt[nyy-1,:]-sbt[nyy-2,:]),np.sum(sbc[nyy-1,:]-sbc[nyy-2,:]))

ft = FortranFile('random02000','r')
x = ft.read_reals()
print (x[0],x[1],x[100], np.size(x))
s1 = np.reshape(x[0:nxx*nyy],(nyy,nxx))
s2 = np.reshape(x[nxx*nyy:2*nxx*nyy],(nyy,nxx))
sbt = (H1*s1 + H2*s2)/(H1+H2)
sbc = s1 - s2
plt.plot(s1[:,0],'or')
plt.figure(4)
plt.contourf(s1,20,cmap=cm.gist_ncar)
plt.colorbar()
plt.figure(5)
plt.contourf(s2,20,cmap=cm.gist_ncar)
plt.colorbar()
print(np.sum(s1[nyy-1,:]-s1[nyy-2,:]),np.sum(s2[nyy-1,:]-s2[nyy-2,:]))
print(np.sum(sbt[nyy-1,:]-sbt[nyy-2,:]),np.sum(sbc[nyy-1,:]-sbc[nyy-2,:]))

plt.show()
