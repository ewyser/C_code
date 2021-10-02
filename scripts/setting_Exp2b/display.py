import math
import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import csv 

# latex 
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
plt.xkcd()
# colorbar
cBtype = 'cividis'
# import data for post-processing 
epII = np.genfromtxt('epII.txt', delimiter=',')
xp   = np.genfromtxt('xp_883.txt'  , delimiter=',')
up   = np.genfromtxt('up.txt'  , delimiter=',')
sig  = np.genfromtxt('sig.txt' , delimiter=',')
p    = -(sig[:,0]+sig[:,1]+sig[:,2])/3
du   = up[:,0]*up[:,0]+up[:,1]*up[:,1]+up[:,2]*up[:,2]

print('')
print('o---------------------------------------------o')
print('|          ** Display CPU results **          |')
print('o---------------------------------------------o')
print('Field variable(s): export to fig.png')

du   = np.sqrt(du)
# 
fig, ax = plt.subplots(figsize=(5,3), dpi=80)
im = ax.scatter(xp[:,0], xp[:,2], s=1, c=du, alpha=1.0, cmap=cBtype)
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel(r'$x$ [m]')
plt.ylabel(r'$z$ [m]')
divider = make_axes_locatable(ax)
cax = divider.new_vertical(size = '25%', pad = 0.5)
fig.add_axes(cax)
cb=fig.colorbar(im, cax = cax, orientation = 'horizontal',extend='max',pad=0.2,label=r'$\Delta u$ [m]')
plt.savefig('du.png', dpi=300, bbox_inches='tight')
plt.show()
print('                 : ' u'\u2713' ' displacement (du) ')

fig, ax = plt.subplots(figsize=(5,3), dpi=80)
im = ax.scatter(xp[:,0], xp[:,2], s=1, c=epII, alpha=1.0, cmap=cBtype)
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel(r'$x$ [m]')
plt.ylabel(r'$z$ [m]')
divider = make_axes_locatable(ax)
cax = divider.new_vertical(size = '25%', pad = 0.5)
fig.add_axes(cax)
cb=fig.colorbar(im, cax = cax, orientation = 'horizontal',extend='max',pad=0.2,label=r'$\epsilon_{II}$ [-]')
plt.savefig('epII.png', dpi=300, bbox_inches='tight')
plt.show()
print('                 : ' u'\u2713' ' plastic strain (epII)')

fig, ax = plt.subplots(figsize=(5,3), dpi=80)
im = ax.scatter(xp[:,0], xp[:,2], s=1, c=p/1e3, alpha=1.0, cmap=cBtype)
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel(r'$x$ [m]')
plt.ylabel(r'$z$ [m]')
divider = make_axes_locatable(ax)
cax = divider.new_vertical(size = '25%', pad = 0.5)
fig.add_axes(cax)
cb=fig.colorbar(im, cax = cax, orientation = 'horizontal',extend='max',pad=0.2,label=r'$p$ [kPa]')
plt.savefig('p.png', dpi=300, bbox_inches='tight')
plt.show()
print('                 : ' u'\u2713' ' pressure (p)')