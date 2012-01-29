# -*- coding: utf-8 -*-
"""
Created on Sat Jan 28 16:38:35 2012

@author: jackie
"""

from pylab import *

def prodn(ri):
    # returns the product of (ri-rj) over all j~=i
    # for each i from 0 to n-1
    n=len(ri)
    
    # create matrix whose element [i,j] is r_i - r_j    
    (rHorz,rVert)=meshgrid(ri,ri)
    diffs=rVert-rHorz
    
    # replace diagonal (i=j) with ones
    diffs = diffs + eye(n)
    
    return prod(diffs,1) # product along horizontal axis
    
def findroots(ai):
    # ai is array of coefficients for polynomial ai[0] + ai[1]*x + ...
    # returns all (complex) roots of the polynomial
    # (found using Durand-Kerner method)
    
    # remove any trailing zeros from ai (i.e. 0*x^n + x^(n-1) -> x^(n-1) )   
    ind = find(ai==0)
    if len(ind)>0:
        while max(ind)==(len(ai)-1):
            ai=ai[0:len(ai)-1]
            ind = find(ai==0)
    
    # renormalize ai so that p(x) = x^n + ...
    ai = ai/ai[len(ai)-1]
        
    
    # n = degree of polynomial 
    n=len(ai)-1
    nrange = arange(n+1)        
    
    # define polynomial function
    def f(x):
        if size(x)>1:
            # operates on vector x
            [nmat,xmat]=meshgrid(nrange,x)
            [aimat,xmat]=meshgrid(ai,x)
            return sum(aimat * xmat**nmat,1) #sums along horizontal axis
        else:
            # operates on scalar x
            return sum(ai * x ** nrange)

    # initial guesses for the n roots
    ri = (0.4 + 0.9j)**nrange[0:n]
    
    for k in range(100):  # easier than worrying about fractional error
        rnext = ri - f(ri)/prodn(ri)
        
        # determining convergence is hard :(        
        frac_err = abs(f(rnext)-f(ri)) / abs(f(ri))
        frac_err = frac_err - frac_err*(f(ri)==0)
        frac_err = max(frac_err)
        print frac_err
        
        if frac_err == 0:
            break
        
        # print ri
        ri = rnext
    
    return ri
    
    
        

ai = array([0,0,0,-1,5,3]) # f(x) = 3x^5 + 5x^4 - x^3
#ai = array([-5,3,-3,1]) # - example from wikipedia - works! yay!
#ai = array([0,0,1]) # - yields ~10^-30 for both - yay!
r = findroots(ai)
print r