# -*- coding: utf-8 -*-
"""

Examples to test the 1D radiative transfer code
(code in transfer1d.py)

"""

### Control which tests to run
test1 = False # test <D^2> = N<l^2>
test2 = False # compare diffusion times for varying density profiles
application1 = False # simulate clumpy planetary nebulae
application2 = True # simulate supernova shock breakout


from transfer1d import *

Rsun = 6.955e10 # cm
Msun = 1.99e33 # g
c = 3.0e10 # cm/s

if test1:
    """TEST 1: Simple 1-D random walk for N steps"""
    #### purpose: demonstrate that <d^2> = N <l^2>
    #### and plot a random walk
    
    param1=Parameters()
    param1.tau=1000.0
    param1.nphotons=10000
    Nstep = arange(10,110,10)
    param1.Rinner=1.0
    
    ratio = zeros(len(Nstep))
    for j in range(len(Nstep)):
        param1.quitcondition=Nstep[j]    
        photons = run1D(param1)
        rlist = photons.r_list()
        rmsdist = sqrt(mean(rlist**2))
        steplist = photons.rmsstep_list()
        rmsstep = sqrt(mean(steplist**2))
        ratio[j] = rmsdist/rmsstep
    
    figure
    subplot(211)
    plot(Nstep,ratio,'*')
    ylabel('<D^2> / <l^2>')
    title('Demonstration of <D^2> = N<l^2>')
    
    hold(True)
    ngrid = linspace(0,100,100)
    ygrid = ngrid**(0.5)
    plot(ngrid,ygrid)
    legend(('Simulations','sqrt(n)'),loc='lower right')
    
    subplot(212)
    resid = ratio - Nstep**0.5
    plot(Nstep,resid,'*')
    hold(True)
    plot([0,100],[0,0])
    xlabel('Number of scatterings')
    ylabel('Residuals')
    savefig('writeup/simpletest.pdf')

if test2:
    """TEST 2: Compare density distributions"""
    param2 = Parameters()
    param2.tau = 10.0
    param2.nphotons = 1
    density_type = ['constant','clumpy4','power2','power50']
    rho=list()
    phot=list()
    
    # make a plot of different density distributions
    for d in density_type:
        param2.density_type = d
        photons = run1D(param2)
        rho.append(photons.rho)
        [phot.append(p) for p in photons.plist]
    rgrid=photons.rgrid
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
    tlist=list()
    
    for d in density_type:
        param2.density_type = d
        photons = run1D(param2)
        times = photons.emergedlist.time_list()
        tlist.append(times)
    
    subplot(313)
    light_time = (param2.Router-param2.Rstar)/c
    est_time = param2.tau
    histrange=[0.,est_time] # play with this
    nbins=40
    for i in range(len(density_type)):
        hold(True)
        times=array(tlist[i])/light_time # time / min time (light travel time for no scatterings)
        hist(times,histtype='step',bins=nbins,range=histrange)
        print 'Median time:', median(times)
    xlabel('Time to emerge / Light-travel time')
    ylabel('# of photons')
    savefig('writeup/densities.pdf')

if application1:
    """APPLICATION 1: Planetary nebulae: smooth vs clumpy"""
    # purpose: demonstrate the code on typical parameters for a
    # planetary nebula - explore morphology of clumps in Balmer emission
    
    param3 = Parameters()
    
    param3.Rstar = 0.01*Rsun # typical white dwarf size
    param3.Rinner = 5.0e17 # make inner cavity half the radius of the nebula (somewhat arbitrary)
    param3.Router = 1.0e18 # typical planetary nebula size ~ 1 ly
    param3.tau = 15. # assume optically thick for Lyman alpha photons
                     # but also assume optically thin for Balmer photons
                     # (this is common for ionized gas - for example HII regions)
    param3.lymanalpha = True
    
    param3.iterator=1
    param3.nphotons = 100000
    
    nclumps = [0,1,2,3]
    density_types = ['clumpy' + repr(n) for n in nclumps]

    figure(figsize=(8,12))
    i=1
    for d in density_types:
        param3.density_type = d 
        
        photons = run1D(param3)
        
        subplot(len(nclumps),2,i)
        plot(photons.rgrid,photons.rho)
        axis([min(photons.rgrid),max(photons.rgrid),0,max(photons.rho)*1.1])    
        if i==1:
            title('Density Profile')
        if i==(len(nclumps)*2-1):
            xlabel('r (cm)')
        
        subplot(len(nclumps),2,i+1)
        lastr = photons.balmerlist.last_scatter_list()
        lastr_emerged = photons.emergedlist.last_scatter_list()
        hist(lastr,range=(param3.Rinner,param3.Router),bins=50,log=True) # histtype='step',
        if i==1:
            title('Radial Balmer Intensity')
        if i==(len(nclumps)*2-1):
            xlabel('r (cm)')
        i+=2 
        savefig('writeup/pn.pdf')


if application2:
    """APPLICATION 2: Supernova shock breakout"""
    ### purpose: use the code to compare times for the flash from
    ### shock breakout to diffuse out through different density profiles:
    ###     stellar wind (1/r^2) and
    ###     stellar atmosphere (1/r^50 <--> scale height ~ R*/50)
    ### and for a range of stellar radii

         
    ## 1: compare total time width for radii from
    ##    50 Rsun (BSG) to 1000 Rsun (RBG)
    ##    make plots of Rsun vs t for both edge types
    
    edgetypes = ['atmosphere','stellar wind']
    Rmin = 50 * Rsun # BSG
    Rmax = 1000 * Rsun # RBG
    Rrange = linspace(Rmin,Rmax,10)
    nphotons = 100000

    mydata = [list(),list()]
    sigtlist = [list(),list()]
    nlist = [list(),list()]
    errtlist = [list(),list()]
    i=0
    for edgetype in edgetypes:
        for Rstar in Rrange:
            results = runShockBreakout(edgetype,Rstar,nphotons)
            mydata[i].append(results)
        sigtlist[i] = array([d.totalsig for d in mydata[i]])/3600.
        nlist[i] = array([len(d.totalt) for d in mydata[i]])
        errtlist[i] = sigtlist[i]/sqrt(nlist[i]-1)
        i+=1

    mintlist = array(Rrange/c)/3600.

    figure()    
    hold(True)
    errorbar(Rrange/Rsun,sigtlist[0],errtlist[0])
    errorbar(Rrange/Rsun,sigtlist[1],errtlist[1])
    plot(Rrange/Rsun,mintlist)
    title('Shock Breakout Duration vs. Stellar Radius')
    xlabel('Stellar Radius (in solar radii)')
    ylabel('Shock Breakout Duration (in hours)')
    legend(('Total Time (Atmosphere)','Total Time (Stellar Wind)','Light Travel Time'))
    savefig('writeup/breakout.pdf')
    
    ## 2: make hist of distributions of total time,
    ##    light travel time, and diffusion time for
    ##    an example near Rcrit (transition between
    ##    light travel- & diffusion-dominated)
    
    d = mydata[1][9]
    ttotal = array(d.totalt)/3600.
    tlight = array(d.lightt)/3600.
    tdiff = array(d.difft)/3600.
    
    figure()
    hold(True)
    maxt=10
    hist(tlight,bins=40,range=(0,maxt),histtype='bar',log=True)
    hist(ttotal,bins=40,range=(0,maxt),histtype='bar',log=True)
    hist(tdiff,bins=40,range=(0,maxt),histtype='step',log=True)
    legend(('Diffusion Time','Light Travel Time','Total Time'))
    xlabel('Time (hours)')
    ylabel('Photon Distribution')
    title('Shock Breakout Lightcurve in a Red Supergiant')
    savefig('writeup/rsg.pdf')