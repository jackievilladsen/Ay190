# -*- coding: utf-8 -*-
"""
Created on Sat Jan 28 13:59:03 2012

@author: jackie
"""

from pylab import *

# properties for root finding
max_frac_err = 1e-10

# properties of Earth's orbit
P = 365.25635 # period of Earth in days
omega = 2*pi/P # orbital frequency (per day)
a = 1 # semi-major axis of Earth in AU
a_km = 1.496e6 # semi-major axis of Earth in km
e = 0.0167 # eccentricity of Earth's orbit
#e = 0.99999 # eccentricity for bad crazy Earth orbit
b = a * sqrt(1-e**2) # semi-minor axis of Earth's orbit

def pos(E):
    x = a * cos(E)
    y = b * sin(E)
    return x,y
    
def makef(t):
    # returns functions f(E) and df_dE(E)
    # defined by f(E) = E - omega * t - e sin E = 0
    # for a specific time t
    
    def f(E):
        return E - omega*t - e*sin(E)
    
    def df_dE(E):
        return 1 - e*cos(E)
    
    return f,df_dE
    
def guessE(t):
    # guesses E as an initial value to feed into the root finder
    return omega * t
    
def findroot((f,df_dx),x0):
    # implements Newton's method to find a root of f near x0    
    
    frac_err = max_frac_err + 1 # ensures we will enter the while loop
    
    xn = x0
    fn = f(xn)
    
    niter = 0
    
    while (frac_err > max_frac_err):
        xnext = xn - fn/df_dx(xn)        
        fnext = f(xnext)
        frac_err = abs((fnext-fn)/fn)
        fn = fnext
        xn = xnext
        niter += 1
    
    return xn,niter

def find_xyE(t):
    # takes time t and returns x, y, and E
    
    #(E,niter) = findroot(makef(t),guessE(t))
    (E,niter) = findroot(makef(t),0)

    (x,y) = pos(E)
    
    return x,y,E,niter

# times (in days) at which to evaluate E,x,y
t0 = 0
t1 = 91
t2 = 182
t3 = 273

(x0,y0,E0,n0) = find_xyE(t0)
(x1,y1,E1,n1) = find_xyE(t1)
(x2,y2,E2,n2) = find_xyE(t2)
(x3,y3,E3,n3) = find_xyE(t3)

print n0,n1,n2,n3

figure
hold(True)
plot(x0,y0,'*')
plot(x1,y1,'x')
plot(x2,y2,'o')
plot(x3,y3,'.')
axis([-1.2,1.2,-1.2,1.2])
legend(('0 days','91 days','182 days','273 days'))
xlabel('x (AU)')
ylabel('y (AU)')
title('Orbit of Earth')
filename = 'orbit' + repr(e) + '.pdf'
savefig(filename)
show()