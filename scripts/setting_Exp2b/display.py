import numpy as nu 
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import make_axes_locatable
import csv 

# latex 
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
# import data for post-processing 
epII = nu.genfromtxt('epII.txt', delimiter=',')
xp   = nu.genfromtxt('xp.txt'  , delimiter=',')
up   = nu.genfromtxt('up.txt'  , delimiter=',')
sig  = nu.genfromtxt('sig.txt' , delimiter=',')
p    = -(sig[:,0]+sig[:,1]+sig[:,2])/3
# 
fig, ax = plt.subplots(figsize=(5,3), dpi=80)
im = ax.scatter(xp[:,0], xp[:,2], s=1, c=epII, alpha=1.0, cmap='cividis')
plt.gca().set_aspect('equal', adjustable='box')
divider = make_axes_locatable(ax)
cax = divider.new_vertical(size = '25%', pad = 0.5)
fig.add_axes(cax)
cb=fig.colorbar(im, cax = cax, orientation = 'horizontal',extend='max',pad=0.2,label=r'$\epsilon_{II}$ [-]')
plt.savefig('epII.png', dpi=300)
plt.show()