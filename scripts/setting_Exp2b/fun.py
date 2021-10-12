import math
import numpy as np
import matplotlib.pyplot as plt 
# compile python3 -B init.py, -B flag hides __pycache__
def topol(nnx,nny,nnz,nelx,nely,nelz,nn):
	nno = nnx*nny*nnz
	nel = nelx*nely*nelz
	gnumbers = (np.arange(0,nno,1,dtype=int)).reshape((nnz,nnx,nny),order='F')	
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

def mesh(nelx,lx,ly,lz0,ni,rho0,coh0,cohr,phi0,phir):
	nely = nelx
	nelz = nely
	lz   = math.ceil(lz0)
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

	xn   = xv.flatten('F')
	yn   = yv.flatten('F')
	zn   = (np.flip(zv,0)).flatten('F')

	e2n  = topol(nnx,nny,nnz,nelx,nely,nelz,nn)

	xB   = [min(xn)+2*h[0],max(xn)-2*h[0],0.0,np.inf,min(yn)+2*h[0],max(yn)-2*h[0]]
	bcx  = np.hstack([np.where(xn<=xB[0]),np.where(xn>=xB[1])]).flatten('F')
	bcy  = np.hstack([np.where(yn<=xB[4]),np.where(yn>=xB[5])]).flatten('F')
	bcz  = np.array(np.where(zn<=xB[2])).flatten('F')					 
	BC   = np.ones((nno*3,1),dtype=int)
	bcx1 = bcx+0*nno
	bcy2 = bcy+1*nno
	bcz3 = bcz+2*nno
	bc   = np.hstack([bcx1,bcy2,bcz3])
	BC[bc] = 0

	fig, ax = plt.subplots(figsize=(5,3), dpi=80)
	plt.gca().set_aspect('equal', adjustable='box')
	ax = plt.axes(projection='3d')
	#ax.scatter3D(x, y, z, c=z, cmap='Blues');
	ax.scatter3D(xn[bcx], yn[bcx], zn[bcx], c=zn[bcx], cmap='Reds');
	ax.scatter3D(xn[bcy], yn[bcy], zn[bcy], c=zn[bcy], cmap='Reds');
	ax.scatter3D(xn[bcz], yn[bcz], zn[bcz], c=zn[bcz], cmap='Blues');
	plt.show()
	# ---------------------------------------------------------------- #
	# material point
	wl   = 0.15*lz0
	x    = np.arange(xB[0]+0.5*h[0]/ni, xB[1], h[0]/ni)
	y    = np.arange(xB[4]+0.5*h[0]/ni, xB[5], h[1]/ni)
	z    = np.arange(xB[2]+0.5*h[2]/ni, lz0  , h[2]/ni)

	xp, yp, zp = np.meshgrid(x, y, z, sparse=False, indexing='ij')
	xp = xp.flatten('F')
	yp = yp.flatten('F')
	zp = zp.flatten('F')
	x  = np.linspace(min(xp),max(yp),200,endpoint=True)
	a  = -1.25
	z  = a*x
	x  = x+0.5*lx
	xf = []
	yf = []
	zf = []
	for mp in range(xp.size):
		for p in range(z.size):
			dx = xp[mp]-x[p]
			dz = zp[mp]-z[p]
			nx = a
			nz = -1
			s  = dx*nx+dz*nz
			if( s>0 ):
				pos = 1
			else:
				pos = 0
			if( zp[mp]<wl ):
				pos = 1	

		if( pos == 1 ):
			xf.append(xp[mp])
			yf.append(yp[mp])
			zf.append(zp[mp])
	xp = np.array(xf)
	yp = np.array(yf)
	zp = np.array(zf)
	"""
	fig, ax = plt.subplots(figsize=(5,3), dpi=80)
	plt.gca().set_aspect('equal', adjustable='box')
	ax = plt.axes(projection='3d')
	#ax.scatter3D(x, y, z, c=z, cmap='Blues');
	ax.scatter3D(xp, yp, zp, c=zp, cmap='Reds');
	plt.show()
	"""
	nmp = xp.size
	lp  = np.ones((nmp,3),dtype=float)*(h[0]/ni)/2
	vp  = np.ones((nmp,1),dtype=float)*(h[0]/ni*h[0]/ni*h[0]/ni)
	mp  = np.ones((nmp,1),dtype=float)*(h[0]/ni*h[0]/ni*h[0]/ni)*rho0
	cohp= np.ones((nmp,1),dtype=float)*coh0
	phip= np.ones((nmp,1),dtype=float)*phi0
	p   = np.where(zp<2*wl) 

	phip[p,0] = phir


	# ---------------------------------------------------------------- #
	# export
	# scalars
	p = np.array([nmp,nn,nno,h[0],h[1],h[2],min(xn),min(yn),min(zn),nnx,nny,nnz])
	np.savetxt("param.txt", p.flatten('F')                                   , fmt="%f", delimiter="\n")
	# mesh-related quantities
	np.savetxt("xn.txt"   , np.hstack([xn,yn,zn]).flatten('F')               , fmt="%f", delimiter="\n")
	np.savetxt("bcs.txt"  , BC.flatten('F')                                  , fmt="%d", delimiter="\n")
	np.savetxt("e2n.txt"  , e2n.flatten('F')                                 , fmt="%d", delimiter="\n")
	# point-related quantities
	np.savetxt("xp.txt"   , np.hstack([xp,yp,zp]).flatten('F')               , fmt="%f", delimiter="\n")
	np.savetxt("lp.txt"   , np.hstack([lp[:,0],lp[:,1],lp[:,2]]).flatten('F'), fmt="%f", delimiter="\n")
	np.savetxt("vol.txt"  , vp                                               , fmt="%f", delimiter="\n")
	np.savetxt("mp.txt"   , mp                                               , fmt="%f", delimiter="\n")
	np.savetxt("cohp.txt" , cohp                                             , fmt="%f", delimiter="\n")
	np.savetxt("phip.txt" , phip                                             , fmt="%f", delimiter="\n")
	
	print("init & export: done")
	return(dx)

