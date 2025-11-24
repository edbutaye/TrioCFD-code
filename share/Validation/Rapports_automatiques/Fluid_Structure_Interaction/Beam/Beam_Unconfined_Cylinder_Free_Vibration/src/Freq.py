import numpy as np
from decimal import *

f_air = np.loadtxt('Air/Freq.txt', unpack=True, usecols=[0]) 
f_water=np.loadtxt('Water/Freq.txt', unpack=True, usecols=[0]) 


DataOut = np.column_stack((f_air, f_water))
np.savetxt('Freq.txt', DataOut)


