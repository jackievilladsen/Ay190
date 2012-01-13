# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 15:10:42 2012

@author: jackie
"""

from pylab import *

def Lnj(j,x,xn):
    n=len(xn)-1
    k=arange(0,n+1) 
    k=find(k!=j)
    Lnj = prod((x-xn[k])/(xn[j]-xn[k]))
    return Lnj
    
def interpLagrange(x,xn,fn):
    # x: single point that we want to interpolate f(x) for
    # xn: array of xn for which we know f(xn)
    # fn: array of f(xn)
    #
    # with n data points, interpolates polynomial of degree n-1
    
    n=len(xn)-1
    p=0
    for j in range(n+1):
       p=p+fn[j]*Lnj(j,x,xn)
    return p
    
def interpLin(x,xn,fn):
    # piecewise linear interpolation
    # xn must be linearly increasing (so we don't have to sort)
    
    l=len(xn)
    j = max(find(xn[0:(l-1)]<=x))
    
    slope = (fn[j+1]-fn[j])/(xn[j+1]-xn[j])
    p = fn[j] + slope * (x-xn[j])
    return p
    
def interpQuad(x,xn,fn):
    # piecewise quadratic interpolation
    # uses 2 points <= our point + 1 point > our point
    # if can't do this, uses 3 leftmost points or 3 rightmost points
    # as long as x is somewhere between min and max xn
    
    l=len(xn)
    if (x<xn[0] or x>xn[l-1]):
        return nan
        
    j = max(find(xn[0:(l-1)]<=x))
    if j<1:
        j=1
    elif j>(l-2):
        j=l-2
    
    xtemp=xn[j-1:j+2]
    ftemp=fn[j-1:j+2]
    p=interpLagrange(x,xtemp,ftemp)
    return p

ceph = loadtxt('ceph.dat',comments='#')
time = ceph[:,0]
mag = ceph[:,1]

dt = 0.01
trange = arange(min(time),max(time)+dt,dt)
#trange = array([0.034])
pLag = zeros(len(trange))
pLin = zeros(len(trange))
pQuad = zeros(len(trange))

for i in range(len(trange)):
    pLag[i]=interpLagrange(trange[i],time,mag)
    pLin[i]=interpLin(trange[i],time,mag)
    pQuad[i]=interpQuad(trange[i],time,mag)
    
plot(trange,pLag,'-b')
plot(trange,pLin,'-r')
plot(trange,pQuad,'-g')
plot(time,mag,'k*')
xlabel('Time (Days)')
ylabel('Magnitude')
title('Interpolation of Cepheid Lightcurve')
legend(('n=8 Polynomial','Piecewise Linear','Piecewise Quadratic','Data'),loc='upper center')
axis([0,1,0,2])
savefig('ceph.pdf')
show()