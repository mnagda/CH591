import numpy as np 
from user_input import *
import matplotlib.pyplot as plt
import os

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

H = np.zeros((3,3), dtype=complex)
init_x = -5
init_c = np.zeros((3), dtype=complex)
init_c[0] = 1
init_c[1] = 0
init_c[2] = 0
total_time = 2100000

mod_c1_sqr = np.zeros(total_time)
mod_c2_sqr = np.zeros(total_time)
mod_c3_sqr = np.zeros(total_time)

plot_x = np.zeros(total_time)

H = pot(init_x)
dt = 0.001
c_dot = np.zeros(3, dtype=complex)
c_dot = -1.j*np.matmul(H,init_c)
c_t = np.zeros(3, dtype=complex)
c_t = init_c
x_t = init_x
for i in range(total_time):
    c_t = c_t + c_dot*dt
    x_t = x_t + v*dt
    H = pot(x_t)
    c_dot = -1.j*np.matmul(H,c_t)
    mod_c1_sqr[i] = np.square(abs(c_t[0]))
    mod_c2_sqr[i] = np.square(abs(c_t[1]))
    mod_c3_sqr[i] = np.square(abs(c_t[2]))
    if x_t > 20: break

print(abs(c_t[0]), abs(c_t[1]), abs(c_t[2]))

#plt.plot(range(total_time),mod_c1_sqr, label=r'$\rho_{el}^{00}$')
#plt.plot(range(total_time),mod_c2_sqr, label=r'$\rho_{el}^{11}$')
#plt.plot(range(total_time),mod_c3_sqr, label=r'$\rho_{el}^{22}$')
#plt.plot(range(total_time),plot_x, label=r'$\rho_{el}^{00}$')
#plt.xlabel("time (a.u.)")
#plt.ylabel(r'$\rho_{el}$')
#plt.rc('font', size=14)
#plt.legend(loc="center right")
#plt.tight_layout()
#plt.savefig("rho.png",dpi=res)
#plt.show()
