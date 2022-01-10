import numpy as np
import struct
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.gridspec import GridSpec
import glob
import sys


# read parallel information ----------------
file_prl = open('./output/parallel_info.dat','rb')

skip = struct.unpack("f",file_prl.read(4))

npe = int((struct.unpack("f",file_prl.read(4)))[0])
nvar = int((struct.unpack("f",file_prl.read(4)))[0])

skip = struct.unpack("f",file_prl.read(4))

file_prl.close()
#-----------------------------------------



data = np.loadtxt('./output/rms.dat')

time_rms = data[:,0]

average_fields = data[:,1:(nvar+1)]
rms = data[:,(nvar+1):(2*nvar+1)]
Ek = data[:,(2*nvar+1):(2*nvar+4)] / 2.0

Ebtot = (rms[:,4] + rms[:,5] + rms[:,6])/2.0



Ektot = Ek[:,0] + Ek[:,1] + Ek[:,2]


fig = plt.figure()
sub = fig.add_subplot(111)
sub.plot(time_rms, Ektot, label=r'$E_k$')
sub.plot(time_rms, Ebtot, label=r'$E_b$')
sub.legend()
sub.set_xlabel(r'$t$')
plt.show()