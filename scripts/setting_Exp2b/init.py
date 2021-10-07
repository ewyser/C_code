import math
import numpy as np
import matplotlib.pyplot as plt 
# -----------------------------------------------------------------#
# call external functions in module fun -> fun.py in same dir
import fun as init
# -----------------------------------------------------------------#
## main program for init
# non-dimensional constant
nu   = 0.3
ni   = 2
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
dx   = init.mesh(nel,lx,ly,lz,ni,rho0,coh0,cohr,phi0,phir) #nn,nno,h,min(x),min(y),min(z),nnx,nny,nnz
Hp   = H*dx
p    = np.array([g,rho0,psi0,nu,E,Kc,Gc,cohr,Hp,t,te,tg])
# export physics
np.savetxt("phys.dat",p,fmt="%f",delimiter="\n")

