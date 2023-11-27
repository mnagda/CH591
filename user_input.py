## This file contains user provided data/functions needed for numerically exact quantum dynamics performed in evolve_QD.py

import numpy as np

## DVR grid ('position' basis)
ndvr=1001
mass=2000
v = 0.00547722558
hbar = 1

# Potential parameters
AA = 0.02
CC = 0.00035       #V12       
PP = 0.00035       #V23
QQ = 0.00005       #V13
BB = 0.6
DD = 1.0
F1 = -AA*BB
F2 = 0
F3 = AA*BB
sgnF12 = (F1 - F2)/abs(F1 - F2)
sgnF23 = (F2 - F3)/abs(F2 - F3)
sgnF13 = (F1 - F3)/abs(F1 - F3)
sgnv = v/abs(v)
F12 = abs(F1 - F2)
F23 = abs(F2 - F3)
F13 = abs(F1 - F3)


## Evolution parameters
dt=2000.0
total_time=6000.0
nprint=2            ## Prints rho after every nprint steps
nsave=10            ## Prints snapshots of the wavefunction every nsave steps

## Position grid
xgrid=np.linspace(-50,50,ndvr)

####################################################################################################################
## User defined potential as a function of coordinate x
## Takes as input a real number x, and returns a 3x3 matrix V[3,3], which is the potential energy matrix at the provided x.
def pot(x):

    V=np.zeros((3,3))

    V[0,0]=AA*np.tanh(BB*x)
    ##V[0,0]=0.0405
    V[2,2]=-V[0,0]
    #V[2,2]=0.0305
    #V[1,1]=-AA*np.sinh(0.05*x)
    V[1,1] = 0
    #V[1,1]=AA/5*np.cosh(0.05*x)
    V[0,1]=CC*np.exp(-DD*x*x)
    #V[0,1]=0
    V[1,0]=V[0,1]
    V[1,2]=PP*np.exp(-DD*x*x)
    #V[1,2]=0
    V[2,1]=V[1,2]
    V[0,2]=QQ*np.exp(-DD*x*x)
    #V[0,2]=0.0000
    V[2,0]=V[0,2]
    return(V)

####################################################################################################################
## User defined initial wavefunction
def init_psi(n,xgrid):
    psi=np.zeros(3*n,dtype=complex)
    ## Initial momentum value
    k=np.sqrt(2*mass*(0.03))
    sigma=2.0
    ## psi(t=0)=N exp((x-x0)**2/sigma**2) . e(ikx)
    for i in range(n):
        psi[i+n]=np.exp(-(xgrid[i]+10)**2/(sigma**2)) * np.exp(1.j*k*xgrid[i])
    ## Normalizing the wavefunction
    psi=psi/np.sqrt(np.vdot(psi,psi))
    return(psi)