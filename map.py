#!/usr/bin/env python

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
from matplotlib import rc

rc('text', usetex=True)

# read data
data = np.loadtxt('test.dat')

plt.figure(figsize=(8,8), dpi=80)
ax = plt.subplot(111)

(n,m) = data.shape
n1 = n/1
m1 = m/1
#print data[300:(300+10),200:(200+10)]
#ax.set_xlim(300,500)
#ax.set_ylim(200,400)
ax.set_xlim(0,m1)
ax.set_ylim(0,n1)
#ax.pcolor(data[20:(20+n1),20:(20+m1)])
#ax.pcolor(data)
ax.pcolor(data, cmap='YlOrBr')
#plt.minorticks_on()
plt.savefig("dynSpec.png",dpi=80)

plt.show()
