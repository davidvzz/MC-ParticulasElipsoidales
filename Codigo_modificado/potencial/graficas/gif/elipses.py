import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

N=100
for i in range(20):
   x,y,a,b,ang=np.loadtxt("gif.dat",skiprows=i*N,max_rows=N,unpack=True)
