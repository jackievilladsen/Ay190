# -*- coding: utf-8 -*-
"""

Examples to test the 1D radiative transfer code
(code in transfer1d.py)

"""

from transfer1d import *

"""TEST 1: Simple 1-D random walk for N steps"""
#### purpose: demonstrate that <d^2> = N <l^2>
#### and plot a random walk
#param1=Parameters()
#param1.tau=1000.0
#param1.nphotons=10000
#Nstep = arange(10,110,10)
#param1.Rinner=1.0
#
#ratio = zeros(len(Nstep))
#for j in range(len(Nstep)):
##    photons = run1D(nphotons,tau,density_type,kappa,Rstar,Rinner,Router,Nstep[j])    
#    param1.quitcondition=Nstep[j]    
##    photons = run1D(param1)
#    results = run1D(param1)
#    photons = results.photonlist
#    rlist = array([p.r for p in photons])
#    rmsdist = sqrt(mean(rlist**2))
#    steplist = array([p.rmsstep for p in photons])
#    rmsstep = sqrt(mean(steplist**2))
#    ratio[j] = rmsdist/rmsstep
#
#figure
#plot(Nstep,ratio,'*')
#xlabel('Number of scatterings')
#ylabel('RMS distance from center / RMS step size')
#title('Demonstration of <D^2> = N<l^2>')
#
#hold=True
#ngrid = linspace(0,100,100)
#ygrid = ngrid**(0.5)
#plot(ngrid,ygrid)
#legend(('Simulations','sqrt(n)'),loc='lower right')
#savefig('writeup/simpletest.pdf')


"""TEST 2: Compare density distributions"""
param2 = Parameters()

param2.tau = 10.0
param2.nphotons = 1
density_type = ['constant','clumpy4','power2','power50']
rho=list()
phot=list()
for d in density_type:
    param2.density_type = d
#    [photons,usefuldata] = run1D(param2)
#    rho.append(usefuldata[3])
#    results = run1D(param2)
    photons = run1D(param2)
#    photons = results.emergedlist
    rho.append(photons.rho)
    [phot.append(p) for p in photons.plist]
rgrid=photons.rgrid

# make a plot of different density distributions
figure(figsize=(8,12))
subplot(311)
hold(True)
for r in rho:
    plot(rgrid,r)
legend(('Flat','Clumpy (N=4)','Power law (n=-2)','Power law (n=-50)'))
#maxy=max(max(rho[0]),max(rho[1]),max(rho[2]),max(rho[3]))
axis([0.,max(rgrid),0.,max(rho[0])*5.])
xlabel('Distance r from center of nebula (cm)')
ylabel('Density rho(r)')

# plot example photon paths
subplot(312)
hold(True)
maxn=0
for p in phot:
    plot(p.rhist,range(len(p.rhist)))
    maxn = max(maxn,p.nsteps)
#maxn = min(maxn,200)
axis([0.,max(rgrid),0.,maxn*1.1])
xlabel('Distance r from center of nebula (cm)')
ylabel('# of scatterings')

# now plot time distributions for many photons for each dist
param2.nphotons = 10000
#tlist=list()
#
#for d in density_type:
#    param2.density_type = d
#    photons = run1D(param2)
#    times = get_time(photons.emergedlist)
#    tlist.append(times)

subplot(313)
light_time = (param2.Router-param2.Rstar)/c
est_time = param2.tau
histrange=[0.,est_time] # play with this
nbins=40
for i in range(len(density_type)):
    hold(True)
    times=array(tlist[i])/light_time # time / min time (light travel time for no scatterings)
    hist(times,histtype='step',bins=nbins,range=histrange)
xlabel('Time to emerge / Light-travel time')
ylabel('# of photons')
savefig('writeup/densities.pdf')


"""APPLICATION 1: Planetary nebulae: smooth vs clumpy"""
## purpose: demonstrate the code on typical parameters for a
## planetary nebula - explore effects of clumps

#param3 = Parameters()
#kappa = 0.2 # cm^2/g
#Router = 1.0e15 # cm - outer edge of nebula
#Rinner = 1.0e14 #1.0e14 # cm - inner edge of nebula - doesn't matter here
#Rstar = 6.0e10 # cm - this is ~ Rsun - used to determine whether photon
#               # travels thru central cavity or is absorbed by star
#Rstar = 0 # zero probability of being absorbed by a star
#ngrid = 1000 # number of points in radial grid (inclusive)
#nphotons=1000  # may need to process each photon individually, or in small
#             # groups, if I need a high number here, since really long
#             # lists can take up a lot of memory
#density_type = 'constant' # constant density and mfp throughout
#
#emerged = run1D(nphotons,tau,density_type,kappa,Rstar,Rinner,Router,ngrid=1000)




"""APPLICATION 2: Supernova shock breakout"""
### purpose: use the code to compare times for the flash from
### shock breakout to diffuse out with a 1/r^2 density law and
### a 1/r^2 density law (also compare constant density since it is
### the easiest to understand)
