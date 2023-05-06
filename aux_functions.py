# EM shower class

import numpy as np
import matplotlib.pyplot as plt


# Class for single shower particle

class shParticle:
    'Particle of the EM shower (in 1D)'
    E = 0.0
    q = 0
    xmin = 0.0
    xmax = 0.0
    
    def __init__(self,E,q,xmin):
        self.E = E
        self.q = q
        self.xmin = xmin
        self.xmax = xmin
        
    def __str__(self):
        if self.xmin == self.xmax :
            return 'E = %g   q = % g    created at X = %g [X0]' % (self.E, self.q, self.xmin)
        else :
            return 'E = %g   q = % g    propagated from X = %g to %g  [X0]' % (self.E, self.q, self.xmin, self.xmax)
        
    def setlen(self,shlen):
        self.xmax = self.xmin + shlen
 
    def getlen(self):
        return self.xmax - self.xmin
    
    def getsamp(self,dX):
        i1 = int(self.xmin/dX)
        i2 = int(self.xmax/dX)
        return i2-i1


class emShower:
    'Model of EM shower development'

    def __init__(self,E0,q0,x0):
        # initialize with input particle
        self.shower = []
        self.shower.append(shParticle(E0,q0,x0))
        # only propagated particles are counted
        # so set counters to zero at that point
        self.npar = 0  
        self.nch = 0
        self.chlen = 0.0
        self.eloss = 0.0

    def __str__(self):
        return 'Shower with %d particles, %d charged, total track lenght %g X0, energy loss %g GeV' % (self.npar, self.nch, self.chlen, self.eloss)
    
    def develope(self,Ec,dEdx=0.0):
        # consider increasing number of particles in the shower
        for par in self.shower :
            # conversion/radiation length [X0]
            if par.q == 0 :
                intlen = 9.0/7.0
            else :
                intlen = 1.0
               
            # generate conversion/radiation point
            
            shlen = np.random.exponential(intlen)
            
            # check energy loss of particle
            
            elen = 0
            
            if par.q != 0 and dEdx > 0 :
                elen = shlen * dEdx
                if elen > par.E :
                    shlen = par.E/dEdx
                    elen = par.E
            
            # propagate particle to this point
            
            par.setlen(shlen)
            self.npar+=1
            if par.q != 0 :
                self.nch+=1
                self.chlen+=shlen
                self.eloss+=elen

            # Energy at the final point
            
            Eleft = par.E - elen
            
            # Final point - origin for new particles
            
            xnew = par.xmax 
           
            # If above critical energy:
            #   convert gamma to two photons or radiate photon
            if Eleft > Ec :
                E1 = Eleft * np.random.random(1)    #  Very siplified energy splitting
                E2 = Eleft - E1
                if par.q == 0 :
                    self.shower.append(shParticle(E1,+1,xnew))
                    self.shower.append(shParticle(E2,-1,xnew))
                else:
                    self.shower.append(shParticle(E1,par.q,xnew))
                    self.shower.append(shParticle(E2,0,xnew))
                    

#We modify the ecuation so it makes also the fit
    def profile(self,Nbin,Xmax):
        hprof=[]
	
        for par in self.shower:
            if par.q != 0 :
                nx = par.getsamp(Xmax/Nbin)
                for ix in range(nx):
                    hprof.append(par.xmin+ix*Xmax/Nbin)

        return hprof
    
