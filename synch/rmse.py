#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
plt.ion()

#print 'Please, do not forget to set nvar, nens, run and variable!'

#nvar = 1000
#nens = 20
#run = 1000


##psi = np.loadtxt('psi.txt')
fort = np.loadtxt('fort.10')
##print fort.shape
##print fort

ts = fort[0:3999,0]
##print ts
rmse = fort[0:3999,1]
##print rmse
spread = fort[0:3999,2]
##print spread

plt.plot(ts,rmse)
plt.plot(ts,spread)
plt.legend(['RMSE','Spread'], loc='best')
###plt.plot(fort)
##plt.show()
##psi = psi.reshape(-1, nvar, nens, order='C')
##truth = np.loadtxt('truth.1.0')
##print psi.shape

##def rmse_spread(i):
    #fig, ax = plt.subplots(figsize=[6,5])
##    fig, ax = plt.subplots(figsize=[20,4])

    #ax.plot(np.sqrt( ((psi[:, i, :].T-truth[1:, i].T)**2).mean(axis=0)), label='RMSE')   #for every 2 obs (1,3,5...)
##    ax.plot(np.sqrt( ((psi[0:run, i, :].T-truth[0:run, i].T)**2).mean(axis=0)), label='RMSE') #for every 3, 5, 10 obs, change 0:___ in truth
##    ax.plot(psi[0:run, i, :].std(axis=1), label='ensemble spread')

##    ax.legend(loc='best')

##    for tick in ax.xaxis.get_ticklabels():
        #tick.set_fontsize('large')
##        tick.set_fontsize(16)
        #tick.set_fontname('Times New Roman')
        #tick.set_color('blue')
        #tick.set_weight('bold')
##    for tick in ax.yaxis.get_ticklabels():
        #tick.set_fontsize('large')
##        tick.set_fontsize(16)

##    return fig, ax

### put in number of variable here
##rmse_spread(100)


#plt.savefig('testpsifig.png', dpi=300)
plt.title('RMSE and Ensemble spread', fontsize=20)
plt.xlabel('Time steps', fontsize=16)        
##plt.ylabel('Value', fontsize=16)
#plt.ylim((0,1.21))
#plt.xlabel('Timestep')
#plt.ylabel('Value')
plt.savefig('rmse.png')
plt.show()

