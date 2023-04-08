#/bin/usr/python3 

import matplotlib.pyplot as plt
import numpy as np
from aux import*

#We define the fitting function

def y(x,sigma,x0,y0):
	y=y0*np.exp(-(x0/sigma)**2*((x-x0)/x0-np.log(x/x0)))
	return y



Ec = 0.020
dEdx = 0.018
E0 = 5
q0 = -1

myshower = emShower(E0,q0,0.)
print("Initial shower particle: ", myshower.shower[0])
myshower.develope(Ec,dEdx)
print(myshower)

Xmax = 50.
Nprof = 200

shprofile=myshower.profile(Nprof,Xmax)
