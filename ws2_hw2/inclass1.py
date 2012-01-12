# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
/Users/jackie/.spyder2/.temp.py
"""

from pylab import *

def f(x):
    return x**3 - 5*x**2 + x

def df_dx(x):
    return 3*x**2 - 10*x + 1

def diff_fwd(x,f):
    l=len(x)
    dx=x[1]-x[0]
    df_dx_fwd = (f[1:(l-1)]-f[0:(l-2)])/dx
    x_fwd = x[0:(l-2)]
    return x_fwd,df_dx_fwd

def diff_cen(x,f):
    l=len(x)
    dx=x[1]-x[0]
    df_dx_cen = (f[2:(l-1)]-f[0:(l-3)])/(2*dx)
    x_cen = x[1:(l-2)]
    return x_cen,df_dx_cen
    
def err_fwd(dx):
    x=arange(-2,6,dx)
    [x_fwd,df_dx_fwd]=diff_fwd(x,f(x))
    err_fwd = df_dx_fwd - df_dx(x_fwd)
    return x_fwd,err_fwd
    
def err_cen(dx):
    x=arange(-2,6,dx)
    [x_cen,df_dx_cen]=diff_cen(x,f(x))
    err_cen = df_dx_cen - df_dx(x_cen)
    return x_cen,err_cen

[x_fwd2,err_fwd2] = err_fwd(0.2)
[x_cen2,err_cen2] = err_cen(0.2)

[x_fwd1,err_fwd1] = err_fwd(0.1)
[x_cen1,err_cen1] = err_cen(0.1)

figure()
plot(x_cen2,err_cen2,'g-',label='dx=0.2')
plot(x_cen1,err_cen1,'r-',label='dx=0.1')
xlabel('x')
ylabel('Error on df/dx')
title('Error Using Central Differencing')
legend()
savefig('cen.pdf')

figure()
plot(x_fwd2,err_fwd2,'g-',label='dx=0.2')
plot(x_fwd1,err_fwd1,'r-',label='dx=0.1')
xlabel('x')
ylabel('Error on df/dx')
title('Error Using Forward Differencing')
legend()
savefig('fwd.pdf')


