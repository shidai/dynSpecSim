#!/usr/bin/env python

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
from matplotlib import rc

rc('text', usetex=True)

# read data
#data = np.loadtxt('test.dat')
#data = np.loadtxt('acf.dat')
data = np.loadtxt('acffftshift.dat')

plt.figure(figsize=(8,8), dpi=80)
ax = plt.subplot(111)

(n,m) = data.shape
ax.set_xlim(0,m)
ax.set_ylim(0,n)
ax.pcolor(data)
#ax.pcolor(data, cmap='YlOrBr')
#plt.minorticks_on()
#plt.savefig("dynSpec.png",dpi=80)
#plt.savefig("acf.png",dpi=80)
#plt.savefig("acffftshift.png",dpi=80)

#plt.savefig("dynSpec.ps",dpi=80)
#plt.savefig("acf.ps",dpi=80)
plt.savefig("acffftshift.ps",dpi=80)
plt.show()
