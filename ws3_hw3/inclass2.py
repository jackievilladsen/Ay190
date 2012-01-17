# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 17:20:03 2012

@author: jackie
"""

from pylab import *
from scipy import *
from scipy.special import *

def gGL(x):
    # returns the form of the number density integral needed for
    # Gauss-Laguerre quadrature: g(x)e^-x = f(x)
    g = x**2 / (1 + exp(-x))
    return g

def GaussLaguerre(g,n):
    # returns the integral of g(x)*e^-x from 0 to infinity
    # calculated using Gauss-Laguerre quadrature with n nodes
    
    [xi,wi]=l_roots(n)
    Q=sum(wi*g(xi))
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

def calc_ne_total(n):
    # calculates total electron (or positron) density using n nodes
    # for Gauss-Laguerre quadrature
    
    n_e = calcConst() * GaussLaguerre(gGL,n)
    return n_e
    
nmax = 100
nrange = range(2,nmax+1)
n_e = zeros(len(nrange))
i=0

for n in nrange:
    n_e[i] = calc_ne_total(n)
    i=i+1
 
ne_final = n_e[i-1]
percent_error = (n_e-ne_final)/ne_final

print ne_final

figure
semilogy(nrange,abs(percent_error))
xlabel('Number n of nodes in Gauss-Laguerre Quadrature')
ylabel('Fractional Error Compared to n=100')
title('Convergence of Electron Density Integral')
savefig('laguerre.pdf')
show()

