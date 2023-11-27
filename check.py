import numpy as np
import os
from user_input import *

def empty(file):
    return os.path.getsize(file) == 0

c3 = np.sqrt(2*np.pi*PP*PP/(hbar*v*F23) + (np.pi*CC*QQ/(hbar*v))**2/(F12*F13) + (2*np.pi)**(3/2)*CC*QQ*PP/hbar**(3/2)*np.sin((np.pi/4)*(sgnF13 - sgnF12 - sgnF23)*sgnv)/np.sqrt((v**3*F12*F13*F23)))
c1 = np.sqrt(2*np.pi*CC*CC/(hbar*v*F12) + (np.pi*QQ*PP/(hbar*v))**2/(F23*F13) + (2*np.pi)**(3/2)*CC*QQ*PP/hbar**(3/2)*np.sin((np.pi/4)*(-sgnF13 + sgnF23 + sgnF12)*sgnv)/np.sqrt((v**3*F12*F13*F23)))
c2 = np.sqrt(1 - c1*c1 - c3*c3)

print(c1, c2, c3)
