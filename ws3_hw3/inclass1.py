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

def plot_result(xmin,xmax,f,hrange,ans,plotname,filename):
    int_h_trap = zeros(len(hrange))
    int_h_simp = zeros(len(hrange))
    
    for i in range(len(hrange)):
        int_h_trap[i] = int_trap(xmin,xmax,f,hrange[i])
        int_h_simp[i] = int_simp(xmin,xmax,f,hrange[i])
        
    # what happens to the error when h -> h/2?
    err_trap_ratio = (int_h_trap[0]-ans)/(int_h_trap[3]-ans)
    err_simp_ratio = (int_h_simp[0]-ans)/(int_h_simp[3]-ans)
    
    print 'When h decreases by a factor of', hrange[0]/hrange[3]
    
    print 'Error decreases by a factor of', err_trap_ratio, '(trapezoid rule)'
    print err_simp_ratio, '(Simpson\'s rule)'

    plot(hrange,int_h_simp,'-b',label='Simpson\'s Rule')
    plot(hrange,int_h_trap,'-r',label='Trapezoidal Rule')
    plot([0,max(hrange)],[ans,ans],'-g',label='Analytic Answer')
    axis([0,max(hrange),ans-0.05,ans+0.01])
    title(plotname)
    xlabel('Step size h')
    ylabel('Result of integral')
    legend(loc='lower right')
    savefig(filename)

    
xmin = 0
xmax = pi
hrange=(xmax-xmin)/arange(3,1000)

f = sin
ans = 2
plotname='Integration of sin x from 0 to pi'
filename='sinx.pdf'
plot_result(xmin,xmax,f,hrange,ans,plotname,filename)

f=xsinx
ans=pi
plotname = 'Integration of x*sin x from 0 to pi'
filename='xsinx.pdf'
plot_result(xmin,xmax,f,hrange,ans,plotname,filename)

