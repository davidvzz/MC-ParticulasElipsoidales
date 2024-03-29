from matplotlib import pylab as plt
import numpy as np
import os

print(os.getcwd())

n=np.loadtxt("EvsN.dat",usecols=0)
energ=np.loadtxt("EvsN.dat", usecols=1)

plt.plot(n,energ)
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.ylabel("U / kT")
plt.xlabel("Pasos MC")
plt.savefig("EvsN")

