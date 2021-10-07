import math
import numpy as np
import matplotlib.pyplot as plt 

def topol(nnx,nny,nnz,nelx,nely,nelz,nn):
	nno = nnx*nny*nnz
	nel = nelx*nely*nelz
	gnumbers = (np.arange(1,nno+1,1,dtype=int)).reshape((nnz,nnx,nny),order='F')	
	for k in range(nny):
		gnumbers[:,:,k]=np.flipud(gnumbers[:,:,k])
	g_num  = np.zeros((nel,nn),dtype=int)
	iel    = 1
	for k in range(nely): 
		for i in range(nelx):
			for j in range(nelz):
				if(i>0 and i<(nelx-1) and j>0 and j<(nelz-1) and k>0 and k<(nely-1)):
					g_num[iel,0 ]= gnumbers[j-1,i-1,k-1]
					g_num[iel,1 ]= gnumbers[j-0,i-1,k-1]
					g_num[iel,2 ]= gnumbers[j+1,i-1,k-1]
					g_num[iel,3 ]= gnumbers[j+2,i-1,k-1]
					g_num[iel,4 ]= gnumbers[j-1,i  ,k-1]
					g_num[iel,5 ]= gnumbers[j-0,i  ,k-1]
					g_num[iel,6 ]= gnumbers[j+1,i  ,k-1]
					g_num[iel,7 ]= gnumbers[j+2,i  ,k-1]
					g_num[iel,8 ]= gnumbers[j-1,i+1,k-1]
					g_num[iel,8 ]= gnumbers[j-0,i+1,k-1]
					g_num[iel,10]= gnumbers[j+1,i+1,k-1]
					g_num[iel,11]= gnumbers[j+2,i+1,k-1]
					g_num[iel,12]= gnumbers[j-1,i+2,k-1]
					g_num[iel,13]= gnumbers[j-0,i+2,k-1]
					g_num[iel,14]= gnumbers[j+1,i+2,k-1]
					g_num[iel,15]= gnumbers[j+2,i+2,k-1]
					g_num[iel,16]= gnumbers[j-1,i-1,k  ]
					g_num[iel,17]= gnumbers[j-0,i-1,k  ]
					g_num[iel,18]= gnumbers[j+1,i-1,k  ]
					g_num[iel,19]= gnumbers[j+2,i-1,k  ]
					g_num[iel,20]= gnumbers[j-1,i  ,k  ]
					g_num[iel,21]= gnumbers[j-0,i  ,k  ]
					g_num[iel,22]= gnumbers[j+1,i  ,k  ]
					g_num[iel,23]= gnumbers[j+2,i  ,k  ]
					g_num[iel,24]= gnumbers[j-1,i+1,k  ]
					g_num[iel,25]= gnumbers[j-0,i+1,k  ]
					g_num[iel,26]= gnumbers[j+1,i+1,k  ]
					g_num[iel,27]= gnumbers[j+2,i+1,k  ]
					g_num[iel,28]= gnumbers[j-1,i+2,k  ]
					g_num[iel,29]= gnumbers[j-0,i+2,k  ]
					g_num[iel,30]= gnumbers[j+1,i+2,k  ]
					g_num[iel,31]= gnumbers[j+2,i+2,k  ]
					g_num[iel,32]= gnumbers[j-1,i-1,k+1]
					g_num[iel,33]= gnumbers[j-0,i-1,k+1]
					g_num[iel,34]= gnumbers[j+1,i-1,k+1]
					g_num[iel,35]= gnumbers[j+2,i-1,k+1]
					g_num[iel,36]= gnumbers[j-1,i  ,k+1]
					g_num[iel,37]= gnumbers[j-0,i  ,k+1]
					g_num[iel,38]= gnumbers[j+1,i  ,k+1]
					g_num[iel,39]= gnumbers[j+2,i  ,k+1]
					g_num[iel,40]= gnumbers[j-1,i+1,k+1]
					g_num[iel,41]= gnumbers[j-0,i+1,k+1]
					g_num[iel,42]= gnumbers[j+1,i+1,k+1]
					g_num[iel,43]= gnumbers[j+2,i+1,k+1]
					g_num[iel,44]= gnumbers[j-1,i+2,k+1]
					g_num[iel,45]= gnumbers[j-0,i+2,k+1]
					g_num[iel,46]= gnumbers[j+1,i+2,k+1]
					g_num[iel,47]= gnumbers[j+2,i+2,k+1]
					g_num[iel,48]= gnumbers[j-1,i-1,k+2]
					g_num[iel,49]= gnumbers[j-0,i-1,k+2]
					g_num[iel,50]= gnumbers[j+1,i-1,k+2]
					g_num[iel,51]= gnumbers[j+2,i-1,k+2]
					g_num[iel,52]= gnumbers[j-1,i  ,k+2]
					g_num[iel,53]= gnumbers[j-0,i  ,k+2]
					g_num[iel,54]= gnumbers[j+1,i  ,k+2]
					g_num[iel,55]= gnumbers[j+2,i  ,k+2]
					g_num[iel,56]= gnumbers[j-1,i+1,k+2]
					g_num[iel,57]= gnumbers[j-0,i+1,k+2]
					g_num[iel,58]= gnumbers[j+1,i+1,k+2]
					g_num[iel,59]= gnumbers[j+2,i+1,k+2]
					g_num[iel,60]= gnumbers[j-1,i+2,k+2]
					g_num[iel,61]= gnumbers[j-0,i+2,k+2]
					g_num[iel,62]= gnumbers[j+1,i+2,k+2]
					g_num[iel,63]= gnumbers[j+2,i+2,k+2]
				iel+=1
	return(g_num)

def fun(nelx,lx,ly,lz):
	nely = nelx
	nelz = nely
	lz   = math.ceil(lz)
	L    = np.array([lx,ly,lz])
	h    = np.array([L[0]/nelx,L[0]/nelx,L[0]/nelx])

	x    = np.arange(0.0-2.0*h[0], lx+2.0*h[0], h[0])
	y    = np.arange(0.0-2.0*h[0], ly+2.0*h[0], h[0])
	z    = np.arange(0.0-2.0*h[0], lz+2.0*h[0], h[0])

	zv, xv, yv = np.meshgrid(z, x, y, sparse=False, indexing='ij')

	nnz  = xv.shape[0]
	nnx  = xv.shape[1]
	nny  = xv.shape[2]
	nno  = nnz*nnx*nny
	nelx = nnx-1
	nely = nny-1
	nelz = nnz-1
	nel  = nelx*nely*nelz
	nn   = 64

	x    = xv.flatten('F')
	y    = yv.flatten('F')
	z    = (np.flip(zv,0)).flatten('F')

	e2n  = topol(nnx,nny,nnz,nelx,nely,nelz,nn)

	xB   = [min(x)+2*h[0],max(x)-2*h[0],0.0,np.inf,min(y)+2*h[0],max(y)-2*h[0]]
	bcx  = np.hstack([np.where(x<=xB[0]),np.where(x>=xB[1])]).flatten('F')
	bcy  = np.hstack([np.where(y<=xB[4]),np.where(y>=xB[5])]).flatten('F')
	bcz  = np.array(np.where(z<=xB[2])).flatten('F')					 
	

	fig, ax = plt.subplots(figsize=(5,3), dpi=80)
	plt.gca().set_aspect('equal', adjustable='box')
	ax = plt.axes(projection='3d')
	#ax.scatter3D(x, y, z, c=z, cmap='Blues');
	ax.scatter3D(x[bcx], y[bcx], z[bcx], c=z[bcx], cmap='Reds');
	ax.scatter3D(x[bcy], y[bcy], z[bcy], c=z[bcy], cmap='Reds');
	ax.scatter3D(x[bcz], y[bcz], z[bcz], c=z[bcz], cmap='Blues');
	plt.show()

	np.savetxt("test.txt", e2n.flatten('F'), fmt="%d", delimiter="\n")

	return(x,y)
# -----------------------------------------------------------------#
## main program for init
# non-dimensional constant
nu   = 0.3
# dimensional constant
g    = 9.81                                                        # gravitationnal acceleration [m/s^2]
E    = 1.0e6                                                       # young's modulus             [Pa]
Gc   = E/(2*(1+nu))                                                # shear modulus               [Pa]
Kc   = E/(3*(1-2*nu))                                              # bulk modulus                [Pa]
rho0 = 2700                                                        # density                     [kg/m^3]
n0   = 0.0                                                         # initial porosity            [-]
yd   = math.sqrt((Kc+4/3*Gc)/rho0)                                 # elastic wave speed          [m/s]
coh0 = 20.0e3                                                      # cohesion                    [Pa]
phi0 = 20.0*math.pi/180                                            # friction angle              [Rad]
psi0 =  0.0*math.pi/180                                            #
H    = -60e3                                                       # softening modulus           [Pa]
cohr =  4.0e3                                                      # residual cohesion           [Pa]
phir = 7.0*math.pi/180                                             # residual friction angle     [Rad]
t    = 15.0                                                        # simulation time             [s]
te   = 8.0                                                         # elastic loading             [s]
tg   = te/1.5   
# mesh and point initialization
nel  = 40   
lx   = 64.1584                                                     
ly   = lx/4                                                        
lz   = 12.80 
x,y=fun(nel,lx,ly,lz)

