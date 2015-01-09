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
plt.minorticks_on()
ax = plt.subplot(111)

(n,m) = data.shape
n1 = n/1
m1 = m/1
#print n1, m1
#print data[300:(300+10),200:(200+10)]
#ax.set_xlim(300,500)
#ax.set_ylim(200,400)
ax.set_title("Freq: 1369 MHz BW: -256 MHz Length: 3840 s")
x0 = 1369.0+128.0
ax.set_xlim(0,m1)
ax.set_ylim(0,n1)
ticks = np.arange((x0-1450.0)*1024/256,(x0-1200)*1024/256,50.0*1024/256)
labels = [1450,1400,1350,1300,1250]
ax.set_yticks(ticks)
ax.set_yticklabels(labels)
ax.set_xticks(np.arange(0,64,10))
ax.set_xticklabels(np.arange(0,64,10))
ax.set_xlabel("Subintegration",fontsize=14)
ax.set_ylabel("Frequency (MHz)",fontsize=14)
#ax.pcolor(data[20:(20+n1),20:(20+m1)])
ax.pcolor(data, cmap='YlOrBr')
#ax.pcolor(data[100:200,100:300], cmap='YlOrBr')
#plt.minorticks_on()
plt.savefig("dynSpec.png",dpi=80)

plt.show()
