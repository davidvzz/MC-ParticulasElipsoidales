import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import imageio

N=500
nframes=60

def elipse(x,y,a,b,ang,i,N):
   fig = plt.figure(figsize=(5,5))
   ax = fig.add_subplot(111, aspect='equal')
   ax.set_xlim(-0.05, 1.05)
   ax.set_ylim(-0.05, 1.05)
   ax.axis("off")
   for j in range(N):
      elipse=Ellipse(xy=(x[j],y[j]), width=a[j], height=b[j], angle= ang[j])
      ax.add_artist(elipse)
   plt.savefig("frame{}.png".format(i),bbox_inches='tight',dpi=200)
   plt.clf()
   plt.close(fig)

frames=[]

for i in range(nframes):
   x,y,a,b,ang=np.loadtxt("gif.dat",skiprows=i*N,max_rows=N,unpack=True)
   elipse(x,y,a,b,ang,i,N)
   image=imageio.v2.imread('frame{}.png'.format(i))
   frames.append(image)

imageio.mimsave("ejemplo.gif",frames, duration=nframes/7000)
