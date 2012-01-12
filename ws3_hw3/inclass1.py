# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 10:30:05 2012

@author: jackie
"""

from matplotlib import *
from pylab import *

def int_trap(xmin,xmax,f,h):
    # implements trapezoid rule to integrate
    #    
    # xmin, xmax - range to integrate over
    # f - function to integrate
    # h - step size (must divide into xmax-xmin)
    
    xi=arange(xmin+h,xmax,h)
    fint = h*(sum(f(xi))+(f(xmin)+f(xmax))/2)
    return fint
    
def int_simp(xmin,xmax,f,h):
    # implements Simpson's rule to integrate
    #
    # same inputs as int_trap
    
    xi=arange(xmin+h,xmax,h)
    xmid=arange(xmin+h/2,xmax+h/2,h)    
    fint=h/3*(sum(f(xi)) + 2*sum(f(xmid)) + (f(xmin)+f(xmax))/2)
    return fint
    
def xsinx(x):
    return x*sin(x)
    
xmin = 0
xmax = pi
f = sin
hrange=(xmax-xmin)/arange(30,1000)
#hrange = arange(0.003,0.1,0.003)

int_h_trap = zeros(len(hrange))
int_h_simp = zeros(len(hrange))

for i in range(len(hrange)):
    int_h_trap[i] = int_trap(xmin,xmax,f,hrange[i])
    int_h_simp[i] = int_simp(xmin,xmax,f,hrange[i])

plot(hrange,int_h_simp,'-b.',label='Simpson\'s Rule')
plot(hrange,int_h_trap,'-r',label='Trapezoidal Rule')
plot([0,.1],[2,2],'-g')
axis([0,1,1.99,2.01])
show()

ind=find(int_h_simp<2)
plot(ind)
show()

f=