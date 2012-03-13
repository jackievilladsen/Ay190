# -*- coding: utf-8 -*-
"""

A library of functions needed to do 1D radiative transfer
problems for isotropic, wavelength-independent scattering.

"""

from matplotlib import *
from pylab import *
from scipy.interpolate import *

c = 3.0e10 # speed of light (cm/s)

class Parameters():
    # default parameters to call run1D with
    def __init__(self):
        self.nphotons=1000
        self.tau=1.0
        self.density_type='constant'
        self.kappa=0.2 #cm^2/g
        self.Rstar=0.0 # cm - effectively no star, so no prob of absorption
        self.Rinner=1.0e14 # cm
        self.Router = 1.0e15 # cm
        self.quitcondition='done'
        self.ngrid=1000
        self.ntimeout=10000 # number of scatterings before a photon times out
        self.lymanalpha = False # for planetary nebulae - 1/2 probability
                                # of automatic escape from anywhere in
                                # nebula by being converted to Balmer photon

class ReturnData():
    # quantities returned by run1D
    def __init__(self):
        self.photonlist=[]
        self.emergedlist=[]
        self.balmerlist=[]
        self.eatenlist=[]
        self.timedoutlist=[]
        self.rgrid=0.0
        self.rho = 0.0

#def run1D(nphotons,tau,density_type,kappa,Rstar,Rinner,Router,quitcondition='done',ngrid=1000):
def run1D(params):
    
    nphotons = params.nphotons
    tau = params.tau
    density_type = params.density_type
    kappa = params.kappa
    Rstar = params.Rstar
    Rinner = params.Rinner
    Router = params.Router
    quitcondition = params.quitcondition
    ngrid = params.ngrid
    ntimeout = params.ntimeout
    lymanalpha=params.lymanalpha
    
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
            # 'power[N]' - power law (i.e 'power2' -> r^-2)
            # 'clumpy[N]' - N sinusoidal clumps ('clumpy3' -> 3 periods)
        # we want optical depth of roughly one
        # tau = kappa*rho*dr so rho~tau/(kappa*dr)   
        dr = Router-Rinner
        rho_avg=tau/(kappa*dr)
        
        if density_type=='constant':
            rho = rho_avg * ones(len(rgrid))

# outdated by 'power2'
#        elif density_type=='square':
#            # let's do a 1/r^2 density, which is caused by high-vel winds
#            # (vel >> escape velocity)
#            # chosen so that average rho from Rinner to Router is rho_avg
#            rho = rho_avg * Rinner*Router/rgrid**2

        elif density_type[0:5]=='power':
            # dont use for powers below -20ish b/c then it rounds to zero
            a = int(density_type[5:]) # the power
            f = Router/Rinner
            A = (a-1)*(f-1) / (1-f**(1-a))
            rho = rho_avg * A * (Rinner/rgrid)**a

        elif density_type[0:6]=='clumpy':
            # approximate clumpy density with a sinusoid
            # for example, density_type='clumpy10' gives 10 clumps            
            nclumps=int(density_type[6:]) # this an input parameter
            omega = 2*pi*nclumps/dr            
            rho = rho_avg * (1-0.5*cos(omega*(rgrid-Rinner)))
            # the 0.5 is to avoid photons escaping entirely in low rho regions
                
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
            self.timedout = False # in case the photon takes too long and it gets boring            
            self.balmer = False # (for planetary nebulae) in case converted to Balmer photon            
            self.done = False # eaten or emerged or timed out
            self.direction = 0. # in 1D, 0 and pi are allowed
                           # 0 is outwards, pi is inwards
            self.time = (Rinner-Rstar)/c # time for photon to travel out after being emitted
            self.stephist = []
            self.rmsstep = 0.
        
        def randDir(self):
            self.direction = pi*(rand()>0.5)
        
        def reflect(self,rnew):
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
        
        def randStep(self,mfp):
            # step to a new location using current direction and
            # a random step size drawn from an exponential distribution
            # with average step size of current mfp
            
            y=rand()
            l=-1. * mfp * log(1-y) # step size
            
            dr = l * cos(self.direction)
            rnew = self.r + dr
            
            # don't let r be > Router
            # also, only count time it takes to get to Router (not beyond)
            if rnew>Router:
                rnew=Router
                self.time+=(Router-self.r)/c
            else:
                # check if the photon hits the inner edge of the cavity
                # if so, reflect it off
                rnew = self.reflect(rnew)
                if rnew>Router: # check if photon has reflected out of nebula
                    rnew=Router
                    self.time+=((Router-Rinner)+(self.r-Rinner))/c
                else:
                    self.time += l/c
            
            self.stephist.append(l)
            
            # set direction of photon randomly
            self.randDir()
            
            return rnew
        
        def calc_rms_step(self):
            self.rmsstep = sqrt(mean(array(p.stephist)**2))
            
        def scatter(self,mfpInterp):
            if self.done:
                return
            elif lymanalpha:
                p_balmer = 0.5 # alpha_A ~= alpha_B, so Lyman alpha
                               # is converted to Balmer photon by
                               # roughly 1/2 of recombinations
                               # (could look up precise numbers:
                               #  p_balmer = alpha_B / (alpha_A + alpha_B) )
                y=rand()
                if y<p_balmer:
                    self.balmer=True
                    rnew = Router
                    self.time += (Router-self.r)/c
                else:
                    mfp = mfpInterp(self.r)
                    rnew = self.randStep(mfp)
            else:        
                mfp = mfpInterp(self.r)
                rnew = self.randStep(mfp)
            
            self.nsteps += 1
            self.r = rnew
            self.rhist.append(rnew)
            self.emerged = (rnew >= Router) and ~self.balmer
            self.eaten = (rnew < Rinner)
            self.timedout = (self.nsteps>ntimeout)
            self.done = (self.emerged or self.eaten or self.timedout or self.balmer)
            
            if self.done:
                self.calc_rms_step()
                
    class PhotonList:
        
        def __init__(self,nphotons,rgrid,rho):
            self.plist = [photon() for _ in range(nphotons)]
            self.allDone = False
            self.nphotons = nphotons
            
            # a memory of where the photons came from            
            self.rgrid = rgrid
            self.rho = rho
            
            # only define these attributes at end to save computation time
            self.nemerged = 0
            self.nbalmer = 0
            self.neaten = 0
            self.ntimedout = 0
            self.emergedlist = []
            self.balmerlist = []
            self.eatenlist = []
            self.timedoutlist = []
            
        
        def scatter_all(self):
            
            
        def separatePhotons(self):
            # separate photons into four categories:
            # emerged by diffusion
            # timed out
            # absorbed by star (eaten)
            # and converted to Balmer photons (for planetary nebulae)
            
            # make a list of the photons that emerged from the nebula
            ind_emerged = find([p.emerged for p in self.plist])
            self.emergedlist = [self.plist[i] for i in ind_emerged]
            self.nemerged = len(ind_emerged)            
            if self.nemerged!=0:
                print self.nemerged, 'photons emerged'         
            
            # make a list of the photons that got eaten by the star
            ind_eaten = find([p.eaten for p in self.plist])
            self.eatenlist = [self.plist[i] for i in ind_eaten]
            self.neaten = len(ind_eaten)       
            if self.neaten!=0:
                print self.neaten, 'photons eaten by star'
            
            # make a list of the photons that got eaten by the star
            ind_timedout = find([p.timedout for p in self.plist])
            self.timedoutlist = [self.plist[i] for i in ind_timedout]
            self.ntimedout = len(ind_timedout)
            if self.ntimedout!=0:
                print self.ntimedout, 'photons timed out after', repr(ntimeout), 'iterations'
            
            # make a list of the photons that were converted to Balmer photons
            ind_balmer = find([p.balmer for p in self.plist])
            self.balmerlist = [self.plist[i] for i in ind_balmer]
            self.nbalmer = len(ind_balmer)
            if self.nbalmer!=0:
                print self.nbalmer, 'Balmer photons'
                
        def nsteps_list(self):
            return array([p.nsteps for p in self.plist])

        def time_list(self):
            return array([p.time for p in self.plist])

#    def separatePhotons(photonlist,returnMe):
#        # separate photons into four categories:
#        # emerged by diffusion
#        # timed out
#        # absorbed by star (eaten)
#        # and converted to Balmer photons (for planetary nebulae)
#                
#        # make a list of the photons that emerged from the nebula
#        ind_emerged = find([p.emerged for p in photonlist])
#        emerged_list = [photonlist[i] for i in ind_emerged]
#        if len(ind_emerged)!=0:
#            print len(ind_emerged), 'photons emerged'         
#        
#        # make a list of the photons that got eaten by the star
#        ind_eaten = find([p.eaten for p in photonlist])
#        eaten_list = [photonlist[i] for i in ind_eaten]
#        if len(ind_eaten)!=0:
#            print 'Warning:', len(ind_eaten), 'photons eaten by star'
#        
#        # make a list of the photons that got eaten by the star
#        ind_timedout = find([p.timedout for p in photonlist])
#        timedout_list = [photonlist[i] for i in ind_timedout]
#        if len(ind_timedout)!=0:
#            print 'Warning:', len(ind_timedout), 'photons timed out after', repr(ntimeout), 'iterations'
#        
#        # make a list of the photons that were converted to Balmer photons
#        ind_balmer = find([p.balmer for p in photonlist])
#        balmer_list = [photonlist[i] for i in ind_balmer]
#        if len(ind_balmer)!=0:
#            print len(ind_balmer), 'Balmer photons'       
#
#        returnMe.photonlist = photonlist
#        returnMe.balmerlist = balmerlist
#        returnMe.eatenlist = eatenlist
#        returnMe.emergedlist = emergedlist
#        returnMe.timedoutlist = timedoutlist
    
    ## here's where we use all the other functions
    
    rgrid=makeGrid()        # make the grid of radii where we will define density
    rho=makeDensity(rgrid)  # calculate the density on that grid
    mfp_grid = rho2mfp(rho) # calculate the mean free path on that grid
    mfpInterp = make_mfpInterp(rgrid,mfp_grid) # make a function to interpolate mfp anywhere
    
    # make a bunch of photons
    photonlist = [photon() for _ in range(nphotons)]
    
    if quitcondition=='done':
        # loop until all photons are absorbed or emerged
        # return only emerged photons        
        
        # iterate through scatterings until photons have all met their end
        allDone = False
        i=0 # counter is just to let me watch progress
        while not allDone:
            [p.scatter(mfpInterp) for p in photonlist]
            allDone = prod([p.done for p in photonlist])  #True if all photons are done
            i+=1
            if i/100.==floor(i/100.):
                print i
        
        returnMe = ReturnData()
        separatePhotons(returnMe)
        returnMe.rgrid = rgrid
        returnMe.rho = rho
        
        return returnMe
        
#        useful_data = [eaten_list,
#                       timedout_list,
#                       rgrid,
#                       rho,
#                       balmer_list]
#            
#        return emerged_list, useful_data
    
    else:
        # this is for testing <D^2> = N*<l^2>        
        # quitcondition, if not 'done', should be a number of scatterings
        # return all photons after that many scatterings
        for _ in range(quitcondition):
                [p.scatter(mfpInterp) for p in photonlist]
        
        [p.calc_rms_step for p in photonlist]
        
#        for p in photonlist:
#            p.rmsstep = sqrt(mean(array(p.stephist)**2))
        
#        # make sure no photons emerged from the star
#        ind_done = find([p.done for p in photonlist])
#        if len(ind_done)!=0:
#            print 'Warning:', len(ind_done), 'photons eaten or emerged'
        
        returnMe = ReturnData()
        separatePhotons(returnMe)
        returnMe.rgrid = rgrid
        returnMe.rho = rho
        
        

