# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 17:20:03 2012

@author: jackie
"""

from pylab import *
from scipy import *
from scipy.special import *

def g(x):
    # returns the form of the number density integral needed for
    # Gauss-Legendre quadrature (i.e. W(x)=1)
    g = x**2 / (1 + exp(x))
    return g
    
def GL1(f,n):
    # returns the integral of f(x) from -1 to 1 using n nodes
    # for Gauss-Legendre integration
    
    [xi,wi]=p_roots(n)
    Q1 = sum(wi*f(xi))
    return Q1

def GaussLegendre(g,n,a,b):
    # returns the integral of g(x) from a to b
    # calculated using Gauss-Legendre quadrature with n nodes

    # transform:
    # x -> t (using t = cx+d => x=(t-d)/c )
    # a -> -1
    # b -> 1
    # g(x)dx -> f(t)dt
    c = 2/(b-a)
    d = -(a+b)/(b-a)
    
    def f(t):
        # transform from g(x)dx to f(t)dt
        # using t = c*x + d
        # so f(t) = g(x) * dx/dt
        
        dx_dt = 1/c
        x=(t-d)/c
        ft = g(x)*dx_dt
        return ft
        
    Q=GL1(f,n)
    return Q

def calcConst():
    # calculates the constant of integration (in cm^-3) on the front
    # of the expression for electron density:
    # n_e(+/-) = A * integral ( x^2 dx / (e^x + 1))
    # where A = 8*pi * (kT)^3 / (2*pi * hbar * c)^3
    # assuming kT = 20 MeV
    
    ergs_per_MeV = 1.60217646e-6 # ergs/MeV
    kT_MeV = 20
    kT = kT_MeV * ergs_per_MeV
    
    hbar = 1.05457148e-27 # cgs
    c = 2.99792458e10 # cm/s
    
    A = 8*pi * (kT/(2*pi*hbar*c))**3
    return A

def calc_ne_ab(Ea,Eb,n):
    # calculates density of electrons with energies between
    # Ea and Eb (in MeV) using n nodes
    # for Gauss-Legendre quadrature
    
    kT = 20. # MeV
    a = Ea/kT
    b = Eb/kT
    
    n_e = calcConst() * GaussLegendre(g,n,a,b)
    return n_e
    
n = 70
dE = 5.
Emax = 250.
nsteps = Emax/dE
Erange = linspace(0.,Emax-dE,nsteps)
l=len(Erange)
n_e = zeros(l)

for i in range(l):
    Ea = Erange[i]
    Eb = Ea + dE
    n_e[i] = calc_ne_ab(Ea,Eb,n)

dn_dE = n_e / dE
 
figure
plot(Erange,dn_dE)
xlabel('Electron Energy (MeV)')
ylabel('dn/dE')
title('Electron Energy Distribution')
savefig('legendre.pdf')
show()

print sum(dn_dE * dE)