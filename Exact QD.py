## Code to perform numerically exact quantum dynamics by solving the eigen functions and energies of the full Hamiltonian.
## Requires file user_input.py, that contains simulation parameters and function for the potential energy and initial conditions.
## This code written specifically for a 3 level problem coupled to a coordinate x.

import numpy as np
import matplotlib.pyplot as plt
from user_input import *
import os


def empty(file):
    return os.path.getsize(file) == 0
#################################################################
## Constructs the Hamiltonian matrix in a DVR basis
def construct_Hamiltonian(n,xgrid,mass):
    dq=xgrid[1]-xgrid[0]
    KE_dvr=np.zeros((n,n))
    ## Kinetic energy in DVR basis. See Appendix A of Colbert and Miller JCP 96, 1982 (1991)
    for i in range(n):
        for j in range(n):
            KE_dvr[i,j] = hbar ** 2 / (2 * mass * dq ** 2) * (-1.0) ** (i - j)
            if (i == j):
                KE_dvr[i,j]*= np.pi ** 2 / 3.0
            else:
                KE_dvr[i,j] *= 2.0 / (i - j) ** 2
    Hamil=np.zeros((3*n,3*n))
    Hamil[:n,:n]=KE_dvr
    Hamil[n:2*n,n:2*n]=KE_dvr
    Hamil[2*n:3*n,2*n:3*n]=KE_dvr
    V11=np.zeros(n)
    V22=np.zeros(n)
    V33=np.zeros(n)
    #UU=np.zeros((3,3,n))
    ## Adding potential energy to the diagonal of the Hamiltonian matrix
    for i in range(n):
        V=pot(xgrid[i])
        Hamil[i,i]     += V[0,0]
        Hamil[i+2*n,i+2*n] += V[2,2]
        Hamil[i,i+2*n]   += V[0,2]
        Hamil[i+2*n,i]   += V[2,0]
        Hamil[i,i+n] += V[0,1]
        Hamil[i+n,i] += V[1,0]
        Hamil[i+n,i+2*n] += V[1,2]
        Hamil[i+2*n,i+n] += V[2,1]
        Hamil[i+n,i+n] += V[1,1]
        
        ## Diabatic energies
        V11[i]=Hamil[i,i]
        V22[i]=Hamil[i+n,i+n]
        V33[i]=Hamil[i+2*n,i+2*n]

    # Plots the diabatic potential energy surfaces
    #plt.plot(xgrid,V11*100)
    #plt.plot(xgrid,V22*100)
    #plt.plot(xgrid,V33*100)
    #plt.title("Diabatic Potentials")
    #plt.xlabel("R")
    #plt.xlim([-10,10])
    #plt.ylabel(r'$V_{d} (x 10^{-2})$')
    #plt.tight_layout()
    #plt.savefig("V_ad.png", dpi=res)
    #plt.show()

    return (Hamil)
#################################################################
## Finds eigen-fns and eigen-energies for the full Hamiltonian (electronic+nuclear) H.eig_fn(x,i) = eig_en(i).eig_fn(x,i)
## Uses inbuilt Python library np.linalg.eigh to find the eigen vectors/values.
def find_eigen_soln():
    Hamil=construct_Hamiltonian(ndvr,xgrid,mass)
    eig_en,eig_vec=np.linalg.eigh(Hamil)
    return(eig_en,eig_vec)
#################################################################
## Calculates the expansion of the initial psi in terms of eigen-fns:
## ci[i]=<eig_fn(i)|psi(t=0)>
def compute_ci(psi,eig_fn):
    ci=np.zeros(3*ndvr,dtype=complex)
    for i in range(3*ndvr):
        #ci[i]=np.vdot(eig_fn[:,i],psi)
        ci[i]=sum(psi*eig_fn[:,i])
    return(ci)
#################################################################

## Evolves psi to time t. Note in this method, no time step dt is required. Evolution is numerically exact for any time t.
## psi(x,t)=\sum_i c(i) e^(-i.t.eig_en(i)/hbar) eig_fn(x,i) 
## where eig_fn(x,i) and eig_en(i) are the ith eigenfunctions of the total Hamiltonian: H.eig_fn(x,i) = eig_en(i).eig_fn(x,i)
## c(i) are the quantum coefficients at t=0: c(i)=<eig_fn(i)|psi(t=0)>
def evolve(ci,eig_en,eig_fn,t):
    psi=np.zeros(3*ndvr,dtype=complex)
    for i in range(3*ndvr):
        psi[:]+=ci[i]*eig_fn[:,i]*np.exp(-1.j*eig_en[i]*t/hbar)
    return(psi)
#################################################################

## Calculates diabatic density matrix psi_d and traces over the position to get a 3x3 electronic matrix rho in diabatic basis
def compute_rho(psi):
    rho_ad=np.zeros((3,3),dtype=complex)
    rho_d=np.zeros((3,3),dtype=complex)
    chi_ad=np.zeros(3,dtype=complex)
    chi_d=np.zeros(3,dtype=complex)
    psi_ad=np.zeros(3*ndvr,dtype=complex)
    psi_d=np.zeros(3*ndvr,dtype=complex)
    for i in range(ndvr):
      chi_d[0]=psi[i]
      chi_d[1]=psi[i+ndvr]
      chi_d[2]=psi[i+2*ndvr]
      ## Transforms wavefunction from diabatic basis to adiabatic basis
      #chi_ad=np.transpose(UU[:,:,i])@chi_d
      #psi_ad[i]=chi_ad[0]
      #psi_ad[i+ndvr]=chi_ad[1]
      #psi_ad[i+2*ndvr]=chi_ad[2]

      psi_d[i]=chi_d[0]
      psi_d[i+ndvr]=chi_d[1]
      psi_d[i+2*ndvr]=chi_d[2]

      ## Adiabatic electronic density matrix: rho[i,j] = <\chi_ad[i] | \chi_ad[j]>
      #rho_ad[0,0]+=abs(chi_ad[0])**2
      #rho_ad[1,1]+=abs(chi_ad[1])**2
      #rho_ad[2,2]+=abs(chi_ad[2])**2
      #rho_ad[0,1]+=np.conj(chi_ad[0])*chi_ad[1]
      #rho_ad[1,0]+=np.conj(chi_ad[1])*chi_ad[0]
      #rho_ad[0,2]+=np.conj(chi_ad[0])*chi_ad[2]
      #rho_ad[2,0]+=np.conj(chi_ad[2])*chi_ad[0]
      #rho_ad[1,2]+=np.conj(chi_ad[1])*chi_ad[2]
      #rho_ad[2,1]+=np.conj(chi_ad[2])*chi_ad[1]

      # Diabatic electronic density matrix: rho[i,j] = <\chi_d[i] | \chi_d[j]>
      rho_d[0,0]+=abs(chi_d[0])**2
      rho_d[1,1]+=abs(chi_d[1])**2
      rho_d[2,2]+=abs(chi_d[2])**2
      rho_d[0,1]+=np.conj(chi_d[0])*chi_d[1]
      rho_d[1,0]+=np.conj(chi_d[1])*chi_d[0]
      rho_d[0,2]+=np.conj(chi_d[0])*chi_d[2]
      rho_d[2,0]+=np.conj(chi_d[2])*chi_d[0]
      rho_d[1,2]+=np.conj(chi_d[1])*chi_d[2]
      rho_d[2,1]+=np.conj(chi_d[2])*chi_d[1]
      
    return(rho_d,psi_d)
#################################################################

## Plotting parameters
res=600
SMALL_SIZE = 14
MEDIUM_SIZE = 14
BIGGER_SIZE = 14

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

#################################################################

## Finds the numerically exact eigen functions and eigen values.
eig_en,eig_fn=find_eigen_soln()
#print("first five eigen values = ",eig_en[:5])

## Sets the initial wave function provided by the user in the file user_input.py
psi=init_psi(ndvr,xgrid)
psi_init=psi
## Calculates the expansion coefficients ci: psi(t=0) = psi(x,t=0)=\sum_i ci(i) eig_fn(x,i)
ci=compute_ci(psi,eig_fn)

max_y=np.max(abs(psi))

## Evolve wavefunction for the user defined steps. Parametesr dt, total_time in the file user_input.py
nsteps=int(total_time/dt)
tim=np.zeros(nsteps)
tim_sav=np.zeros(nsteps)
psi_sav=np.zeros((3*ndvr,nsteps),dtype=complex)
rho_d=np.zeros((3,3,nsteps),dtype=complex)
j=0
for i in range(nsteps):
    psi=evolve(ci,eig_en,eig_fn,i*dt)           ## Evolves psi numerically exact to the time t
    rho_d[:,:,i],psi_d=compute_rho(psi)       ## Calculates the diabatic density matrix
    tim[i]=i*dt
    # Every nprint steps, plot the wavefunction as a function of time to create a movie
    if(i%nprint==0):
        plt.clf()
        plt.subplot(3,1,1)
        plt.ylim((0,max_y))
        plt.title(str(i*dt))
        plt.plot(xgrid,abs(psi_d[:ndvr]))
        plt.subplot(3,1,2)
        plt.ylim((0,max_y))
        plt.plot(xgrid,abs(psi_d[ndvr:2*ndvr]))
        plt.subplot(3,1,3)
        plt.ylim((0,max_y))
        plt.plot(xgrid,abs(psi_d[2*ndvr:3*ndvr]))
        plt.pause(0.00001)            ## Time between 2 frames

    ## Every nprint steps, store the wavefunction 
    if(i%nsave==0):
      tim_sav[j]=i*dt
      psi_sav[:,j]=psi_d
      j+=1
plt.show()

###########################################################

## Plotting electronic desnsity matrix (absolute values)
#plt.plot(tim,abs(rho_d[0,0,:]),label=r'$\rho_{el}^{00}$')
#plt.plot(tim,abs(rho_d[1,1,:]),label=r'$\rho_{el}^{11}$')
#plt.plot(tim,abs(rho_d[0,1,:]),label=r'$|\rho_{el}^{01}|$')
#plt.plot(tim,abs(rho_d[2,2,:]),label=r'$\rho_{el}^{22}$')
#plt.plot(tim,abs(rho_d[0,2,:]),label=r'$\rho_{el}^{02}$')
#plt.plot(tim,abs(rho_d[1,2,:]),label=r'$|\rho_{el}^{12}|$')
#plt.xlabel("time (a.u.)")
#plt.ylabel(r'$\rho_{el}$')
#plt.rc('font', size=14)
#plt.legend(loc="center right")
#plt.tight_layout()
#plt.savefig("rho.png",dpi=res)
#plt.show()

# Plots snapshots of the wavefunction after every nsave steps (nsave defined in file user_input.py)
#for i in range(j):
#  plt.subplot(3,j,i+1)
#  plt.ylim((0,max_y))
#  plt.title("t="+str(tim_sav[i]))
#  #plt.xlabel("R (a.u.)")
#  if(i==0):
#    plt.ylabel(r'$\chi_2$')
#    plt.tick_params(left = False, right = False,labelbottom = False,bottom = False)
#  if(i>0):
#    plt.tick_params(left = False, right = False,labelleft = False ,labelbottom = False,bottom = False)
#  plt.plot(xgrid,abs(psi_sav[ndvr:2*ndvr,i]))
#plt.subplots_adjust(wspace=None, hspace=None)

#  plt.subplot(3,j,j+i+1)
#  plt.ylim((0,max_y))
#  #plt.xlabel("R (a.u.)")
#  plt.xlim([-10,20])
#  #plt.xticks(np.arange(-5, 20, step=10))
#  if(i==0):
#    plt.ylabel(r'$\chi_1$')
#  if(i>0):
#    plt.tick_params(left = False, right = False,labelleft=False)
#  plt.plot(xgrid,abs(psi_sav[:ndvr,i]))
#  
#  plt.subplot(3, j, 2*j+i+1)  # Creating the third subplot in the third row
#  plt.ylim((0, max_y))
#  plt.xlabel("R (a.u.)")
#  if(i==0):
#     plt.ylabel(r'$\chi_3$')
#  if(i>0):
#    plt.tick_params(left = False, right = False,labelleft=False)
#  plt.plot(xgrid, abs(psi_sav[2*ndvr:, i]))
  
#plt.subplots_adjust(wspace=0, hspace=0)
##plt.tight_layout()
##plt.savefig("psi_d.png",dpi=res)
#plt.show()
print(np.sqrt(abs(rho_d[0,0,nsteps-1])), np.sqrt(abs(rho_d[1,1,nsteps-1])), np.sqrt(abs(rho_d[2,2,nsteps-1])))
