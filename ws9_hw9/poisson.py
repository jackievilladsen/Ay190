# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 10:26:07 2012

@author: jackie
"""

from pylab import *
from scipy.interpolate import *

G=6.67e-8

def readData():
    myData=loadtxt('presupernova.dat')
    r=myData[:,2] # ranges from 10^6 to 10^13 - adaptive step size
    rho=myData[:,4] # ranges for 1e10 to 1e-10
    M=myData[:,1] # ranges from 10^30 (initial enclosed mass?) to 10^34
    T=myData[:,3] # ranges from 7e9 to 3e3
    
    return r,rho

def regrid(r,rho,nsteps):
    # set up new, regular r grid
    # and then interpolate rho onto it
    
    rmin=0.
    rmax=1.e9
    rnew=linspace(rmin,rmax,nsteps)
    
    rhospline=splrep(r,rho)
    rhonew=splev(rnew,rhospline)
    
    #deal with r values < min(r)
    def interpLin(rval):
        slope = (rho[1]-rho[0])/(r[1]-r[0])
        intercept = rho[0] - slope*r[0]
        rhoval = slope*rval + intercept
        return rhoval        
    ind=find(rnew<min(r))
    rhonew[ind]=interpLin(rnew[ind])
    
    return rnew,rhonew


def get_dY_dr(r,rho,Y):
    # how to solve with ODEs:
    # z= d(phi)/dr
    # dz/dr + 2/r * z = 4*pi*G*rho
    # (define a vector Y=[phi,z,M])
    z=Y[1]
    dY_dr = zeros(3)
    dY_dr[0] = z   # d(phi)/dr
    dY_dr[1] = 4*pi*G*rho - 2/r * z   # dz/dr
    dY_dr[2] = 4*pi*r**2*rho
    ind=find(isnan(dY_dr))
    dY_dr[ind]=0
    return dY_dr

def stepFwdEuler(r,rho,Y,dr):
    Ynew = Y + dr * get_dY_dr(r,rho,Y)
    return Ynew

def applyBoundary(phi,M,r):
    l=len(r)
    r_outer = r[l-1]
    M_outer = M[l-1]
    phi_outer= phi[l-1]
    phi_outer_true = - G*M_outer / r_outer
    dphi = phi_outer - phi_outer_true
    phi_true = phi - dphi
    return phi_true

def integrateEuler(r,rho):
    
    Y0 = [0,0,0] # assume at center phi=0, dphi/dr=0, Menc=0
    Y=Y0    
    
    l=len(r)
    dr=r[1]-r[0]
    Yarray=zeros((l,3))
    for i in range(l-1):
        Yold=Y
        Y=stepFwdEuler(r[i],rho[i],Yold,dr)
        Yarray[i+1,:]=Y
    
    phi_nobound = Yarray[:,0]
    M=Yarray[:,2]
    phi = applyBoundary(phi_nobound,M,r)
    return phi


[r_uneven,rho_uneven]=readData()
[r,rho]=regrid(r_uneven,rho_uneven,10000)
  
#figure
#hold(True)
#loglog(r_uneven,rho_uneven,linewidth=2)
#loglog(r,rho,'--r',linewidth=2)
#legend(('Model','Interpolated'))
#xlabel('Radius')
#ylabel('Density')
#show()

phi=integrateEuler(r,rho)
figure
plot(r,phi)
show()

