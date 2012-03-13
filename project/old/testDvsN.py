# -*- coding: utf-8 -*-
"""

Bootstrapped to test 1d code

A slight modification to transfer1d.py that only takes Nstep steps for
all photons, then stops, so that we can compare rms distance traveled
in N steps to rms step size.

Must provide tau >> Nstep so no photons step off the edge of the nebula.



"""

from matplotlib import *
from pylab import *
from scipy.interpolate import *

c = 3.0e10 # speed of light (cm/s)


def testD(nphotons,Nstep,tau,ngrid=1000):
    ngrid=1000
    density_type = 'constant'
    Rstar = 0.0
    Rinner = 1.0
    Router = 1.0e15
    kappa=0.2
    
    def makeGrid():
        # produces the grid of radii
        # for 1D, goes from -Router to Router
        #
        # for 2D, if we only do pi=0 to pi/2, how
        # do we handle stepping off the edge? maybe they could
        # just reflect back - this is like assuming all 4 sections
        # are identical
        
        rgrid = linspace(Rinner,Router,ngrid)
        return rgrid
    
    def makeDensity(rgrid):
        # produces a mass density
        # options for density type:
            # 'constant' - constant density
            # 'square' - inverse square law
        
        # we want optical depth of roughly one
        # tau = kappa*rho*dr so rho~tau/(kappa*dr)   
        dr = Router-Rinner
        rho_avg=tau/(kappa*dr)
        
        if density_type=='constant':
            rho = ones(len(rgrid))*rho_avg
        elif density_type=='square':
            # let's do a 1/r^2 density, which is caused by high-vel winds
            # (vel >> escape velocity)
            # chosen so that average rho from Rinner to Router is rho_avg
            rho = Rinner*Router/rgrid**2 * rho_avg
    
        return rho
    
    
    def rho2mfp(rho):
        # converts the grid of densities to a grid of mean free paths
        
        mfp = 1/(kappa * rho) 
        return mfp
    
    def make_mfpInterp(rgrid,mfp_grid):
        # creates a spline interpolation function that we can
        # call repeatedly to sample the mean free path when the photon
        # is at a location in between the points of the density grid
        
        # is there a function to extend spline interpolation to 2d?
        
        mfpspline=splrep(rgrid,mfp_grid)
        def mfpInterp(r):
            mfp=splev(r,mfpspline)
            return mfp
        return mfpInterp
    
    class photon:
        # hold the properties of a single photon:
        # current position and direction, history,
        # and whether it's emerged yet   
        
        def __init__(self):
            self.r = Rinner # holds current position
            self.rhist = [Rinner] # list of all past positions
            self.nsteps = 0 # number of steps it's taken
            self.emerged = False # for when photon escapes Router
            self.eaten = False # in case the photon gets eaten by the star
            self.done = False # eaten or emerged
            self.direction = 0. # in 1D, 0 and pi are allowed
                           # 0 is outwards, pi is inwards
            self.time = 0. # time for photon to travel out after being emitted
            self.stephist = []
            self.rmsstep = 0.
        
        def randDir(self):
            self.direction = pi*(rand()>0.5)
        
        def randStep(self,mfp):
            # step to a new location using current direction and
            # a random step size drawn from an exponential distribution
            # with average step size of current mfp
            
            y=rand()
            l=-1. * mfp * log(1-y) # step size
            
            dr = l * cos(self.direction)
            rnew = self.r + dr
            
            self.stephist.append(l)
      
            # if the photon enters the central cavity, approximate the
            # probability of hitting the star and being absorbed
            # as p(absorbed) = Rstar/Rinner
            # if it does not hit the star (not directly possible in 1D
            # which is why we are using a probability) then have it
            # "transmitted" by reflecting back (as if it came from the
            # other half since we are only simulating one half)
            if rnew < Rinner:
                y=rand()
                if y>(Rstar/Rinner):
                    reflectdist = Rinner-rnew
                    rnew = Rinner + reflectdist # reflect rather than absorb
                    self.time += 2*Rinner/c # account for time spent in cavity
            
            return rnew        
            
        def scatter(self,mfpInterp):
            if self.done:
                return
            
            mfp = mfpInterp(self.r)
            rnew = self.randStep(mfp)
    
            dt = (rnew-self.r)/c
            self.time += dt
            
            self.nsteps += 1
            self.r = rnew
            self.rhist.append(rnew)
            self.emerged = abs(rnew) > Router
            self.eaten = rnew < Rinner
            self.done = self.emerged or self.eaten
            
            if self.done:
                self.rmsstep = sqrt(mean(array(self.stephist)**2))
            
            # set direction of photon randomly, except
            # if photon just scattered off inner boundary,
            # then make sure it's pointing outwards
            self.randDir()
            if rnew == Rinner:
                self.direction = 0.0
    
    rgrid=makeGrid()
    rho=makeDensity(rgrid)
    mfp_grid = rho2mfp(rho)
    mfpInterp = make_mfpInterp(rgrid,mfp_grid)
    
    photonlist = [photon() for _ in range(nphotons)]
    print Nstep ## for testing purposes
    for i in range(Nstep):
        [p.scatter(mfpInterp) for p in photonlist]
        if i/100.==floor(i/100.):
            print i
            
    for p in photonlist:
        p.rmsstep = sqrt(mean(array(p.stephist)**2))
            
    # make sure no photons emerged from the star
    ind_done = find([p.done for p in photonlist])
    print 'emerged:',ind_done
    
    return photonlist