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
        self.kappa=0.2 #cm^2/g # Thomson scattering, for mu=2
        self.Rstar=0.0 # cm - effectively no star, so no prob of absorption
        self.Rinner=1.0e14 # cm
        self.Router = 1.0e15 # cm
        self.quitcondition='done'
        self.ngrid=10000
        self.ntimeout=10000 # number of scatterings before a photon times out
        self.lymanalpha = False # for planetary nebulae - 1/2 probability
                                # of automatic escape from anywhere in
                                # nebula by being converted to Balmer photon
        self.iterator = 100. # how often to print cycle number

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
    iterator = float(params.iterator)
    
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

        elif density_type[0:5]=='power':
            # dont use for powers below -20ish b/c then it rounds to zero
            a = int(density_type[5:]) # the power
            f = Router/Rinner
            A = (a-1)*(f-1) / (1-f**(1-a))
            rho = rho_avg * A * (Rinner/rgrid)**a

        elif density_type[0:6]=='clumpy':
            # approximate clumpy density with a sinusoid
            # for example, density_type='clumpy10' gives 10 clumps  
            # and 'clumpy0' is the same as 'constant'
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
            self.rmsstep = sqrt(mean(array(self.stephist)**2))
        
        def get_last_scatter(self):
            return self.rhist[self.nsteps-1]
           
        def scatter(self,mfpInterp):
            if self.done:
                return
            elif lymanalpha:
                #p_balmer = 2.59/4.18 #  p_balmer = alpha_B / (alpha_A + alpha_B) )
                p_balmer = 0.2                
                y=rand()
                if y<p_balmer and self.stephist!=[]:
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
            self.emerged = (rnew >= Router) and not self.balmer
            self.eaten = (rnew < Rinner)
            self.timedout = (self.nsteps>ntimeout)
            self.done = (self.emerged or self.eaten or self.timedout or self.balmer)
                
    class PhotonList:
        
        def __init__(self,nphotons,rgrid=0.0,rho=0.0):
            self.plist = [photon() for _ in range(nphotons)]
            self.allDone = False
            self.nphotons = nphotons
            
            # a memory of where the photons came from
            # they both remain zero if not defined
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
        
        def makePhotList(self,mylist,rgrid=0.0,rho=0.0):
            # input: mylist is a plain list of photon objects
            # returns: a PhotonList of photon objects
            
            photlist = PhotonList(0,rgrid,rho)
            photlist.plist = mylist
            return photlist   
        
        def scatter_all(self,mfpInterp):
            [p.scatter(mfpInterp) for p in self.plist]
            self.allDone = prod([p.done for p in self.plist])
            
        def finalize(self):
            # calculate final rms step size on photons and            
            # separate photons into four categories:
            #  emerged by diffusion
            #  timed out
            #  absorbed by star (eaten)
            #  and converted to Balmer photons (for planetary nebulae)
            
            [p.calc_rms_step() for p in self.plist]
            
            # make a list of the photons that emerged from the nebula
            ind_emerged = find([p.emerged for p in self.plist])
            emergedlist = [self.plist[i] for i in ind_emerged]
            self.emergedlist = self.makePhotList(emergedlist,self.rgrid,self.rho)            
            self.nemerged = len(ind_emerged)            
            if self.nemerged!=0:
                print self.nemerged, 'photons emerged'         
            
            # make a list of the photons that got eaten by the star
            ind_eaten = find([p.eaten for p in self.plist])
            eatenlist = [self.plist[i] for i in ind_eaten]
            self.eatenlist = self.makePhotList(eatenlist,self.rgrid,self.rho)          
            self.neaten = len(ind_eaten)       
            if self.neaten!=0:
                print self.neaten, 'photons eaten by star'
            
            # make a list of the photons that got eaten by the star
            ind_timedout = find([p.timedout for p in self.plist])
            timedoutlist = [self.plist[i] for i in ind_timedout]
            self.timedoutlist = self.makePhotList(timedoutlist,self.rgrid,self.rho)            
            self.ntimedout = len(ind_timedout)
            if self.ntimedout!=0:
                print self.ntimedout, 'photons timed out after', repr(ntimeout), 'iterations'
            
            # make a list of the photons that were converted to Balmer photons
            ind_balmer = find([p.balmer for p in self.plist])
            balmerlist = [self.plist[i] for i in ind_balmer]            
            self.balmerlist = self.makePhotList(balmerlist,self.rgrid,self.rho)            
            self.nbalmer = len(ind_balmer)
            if self.nbalmer!=0:
                print self.nbalmer, 'Balmer photons'
                
        def nsteps_list(self):
            return array([p.nsteps for p in self.plist])

        def time_list(self):
            return array([p.time for p in self.plist])
        
        def r_list(self):
            return array([p.r for p in self.plist])
        
        def rmsstep_list(self):
            return array([p.rmsstep for p in self.plist])
        
        def last_scatter_list(self):
            return array([p.get_last_scatter() for p in self.plist])

    
    ## here's where we use all the other functions
    
    rgrid=makeGrid()        # make the grid of radii where we will define density
    rho=makeDensity(rgrid)  # calculate the density on that grid
    
    # could move these into the photon or photonlist class for prettiness    
    mfp_grid = rho2mfp(rho) # calculate the mean free path on that grid
    mfpInterp = make_mfpInterp(rgrid,mfp_grid) # make a function to interpolate mfp anywhere
    
    # make a bunch of photons
    photonlist = PhotonList(nphotons,rgrid,rho) 
  
    if quitcondition=='done':       
        # iterate through scatterings until photons have all met their end
        i=0 # counter is just to let me watch progress
        while not photonlist.allDone:
            photonlist.scatter_all(mfpInterp)            
            i+=1
            if i/iterator==floor(i/iterator):
                print i
        
        photonlist.finalize()
    
        return photonlist
    
    else:
        # this is for testing <D^2> = N*<l^2>        
        # quitcondition, if not 'done', should be a number of scatterings
        # return all photons after that many scatterings
        
        for i in range(quitcondition):
            photonlist.scatter_all(mfpInterp)
            if i/iterator==floor(i/iterator):
                print i
       
        photonlist.finalize()
                
        return photonlist
        
def runShockBreakout(edgetype,Rstar,nphotons):    

    param4 = Parameters()
    param4.iterator = 10.
    param4.nphotons = nphotons
    param4.tau = 30. # condition for shock breakout to start
    
    param4.Rstar = Rstar
    dR = 5*Rstar # many times the scale length
    
    if edgetype == 'stellar wind':
        param4.Rinner = param4.Rstar # Rstar is edge of shock
        param4.density_type = 'power2'
    elif edgetype == 'atmosphere':
        param4.Rinner = param4.Rstar * 49/50 #(i.e. one scale height)
        param4.density_type='power50'
    param4.Router = param4.Rinner + dR
        
    photonlist = run1D(param4)
    difft = photonlist.emergedlist.time_list()
           
    dr_light_time = dR/c
    R_light_time = param4.Rstar/c

    #apply light travel time delay
    thetalist = pi*rand(photonlist.nemerged)
    lightt = R_light_time * sin(thetalist)
    totalt = lightt + difft
    
    class returndata:
        def __init__(self,totalt,difft,lightt):
            self.totalt = totalt
            self.difft = difft
            self.lightt = lightt
            self.totalsig = std(totalt)
            self.diffsig = std(difft)
            self.lightsig = std(lightt)
            self.edgetype = edgetype
            self.Rstar = Rstar
            self.mintime = Rstar/c
    
    mydata = returndata(totalt,difft,lightt)
    return mydata
