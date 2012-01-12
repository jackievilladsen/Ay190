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
    n=len(xn)-1
    p=0
    for j in range(n+1):
       p=p+fn[j]*Lnj(j,x,xn)
    return p

ceph = loadtxt('ceph.dat',comments='#')
time = ceph[:,0]
mag = ceph[:,1]

plot(time,mag,'r*')

dt = 0.01
trange = arange(min(time),max(time)+dt,dt)
prange = zeros(len(trange))

for i in range(len(trange)):
    prange[i]=interpLagrange(trange[i],time,mag)
    
plot(trange,prange,'-b')
show()