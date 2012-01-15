# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 15:10:42 2012

@author: jackie
"""

from pylab import *
from scipy.interpolate import *

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
        # or could wrap around since lightcurve is periodic
        
    j = max(find(xn[0:(l-1)]<=x))
    if j == l:
        j=0
    
    xtemp=array([xn[j-1], xn[j], xn[j+1]])
    ftemp = array([fn[j-1], fn[j], fn[j+1]])
    p=interpLagrange(x,xtemp,ftemp)
    return p
    
def df_dx_cen(i,xn,fn):
    # takes the derivative of f(x) at point x_i
    # assumes a periodic function (i.e. f(xn) = f(x0) --> wraps around
    
    n=len(fn)-1
    a=i-1
    b=i+1
    if b>n:
        b=b-n
        period = xn[n]-xn[0]
        dx = xn[b]+period-xn[a]
    else:
        dx=xn[b]-xn[a]
    
    df_dx = (fn[b]-fn[a])/dx
    return df_dx
    
def psi0(z):
    return 2 * z**3 - 3 * z**2 + 1
    
def psi1(z):
    return z**3 - 2 * z**2 + z

def interpHerm(x,xn,fn):
    # implements piecewise cubic Hermite interpolation
    
    l=len(xn)
    i=max(find(xn[0:(l-1)]<=x))
    
    dx= xn[i+1]-xn[i]
    z=(x-xn[i])/dx
    
    H3 = ( fn[i] * psi0(z) + fn[i+1]*psi0(1-z)
         + df_dx_cen(i,xn,fn) * dx * psi1(z)
         - df_dx_cen(i+1,xn,fn) * dx * psi1(1-z) )
    return H3

ceph = loadtxt('ceph.dat',comments='#')
time = ceph[:,0]
mag = ceph[:,1]

dt = 0.01
trange = arange(min(time),max(time)+dt,dt)
#trange = array([0.034])
pLag = zeros(len(trange))
pLin = zeros(len(trange))
pQuad = zeros(len(trange))
pHerm = zeros(len(trange))
pSpline = zeros(len(trange))

l=len(time)
period=time[l-1]-time[0]
time2=zeros(3*l-2)
mag2=zeros(3*l-2)
time2[0:l]=time-period
mag2[0:l]=mag
time2[l-1:2*l-1]=time
mag2[l-1:2*l-1]=mag
time2[(2*l-2):(3*l-2)]=time+period
mag2[2*l-2:3*l-2]=mag
print time2
print mag2
interpSpline = interp1d(time2,mag2,kind=3)

for i in range(len(trange)):
    pLag[i]=interpLagrange(trange[i],time,mag)
    pLin[i]=interpLin(trange[i],time2,mag2)
    pQuad[i]=interpQuad(trange[i],time2,mag2)
    pHerm[i]=interpHerm(trange[i],time2,mag2)
    pSpline[i]=interpSpline(trange[i])

figure
subplot(311)    
plot(trange,pLag,'-r',label='Lagrange (n=8 polynomial)')
plot(time,mag,'k*')
title('Interpolation of Cepheid Lightcurve')
legend(loc='upper center')
axis([0,1,0,3])
subplot(312)
plot(trange,pLin,'-r',label='Piecewise Linear')
plot(time,mag,'k*')
axis([0,1,0,1])
legend(loc='upper left')
ylabel('Magnitude')
subplot(313)
plot(trange,pQuad,'-r',label='Piecewise Quadratic')
plot(time,mag,'k*')
axis([0,1,0,1])
legend(loc='upper left')
xlabel('Time (Days)')
savefig('ceph.pdf')
show()


#subplot(514)
#plot(trange,pHerm,'-r',label='Piecewise Cubic Hermite')
#plot(time,mag,'k*')
#axis([0,1,0,1])
#legend(loc='upper left')
#subplot(515)
#plot(trange,pSpline,'-r',label='Cubic Spline')
#plot(time,mag,'k*')
#xlabel('Time (Days)')
#axis([0,1,0,1])
#legend(loc='upper left')
#savefig('ceph.pdf')
#show()

figure
plot(trange,pHerm,'-b')
plot(trange,pSpline,'-r')
plot(time,mag,'k*')
xlabel('Time (Days)')
ylabel('Magnitude')
axis([0,1,0,1])
title('Interpolation of Cepheid Lightcurve')
legend(('Piecewise Cubic Hermite','Cubic Spline'),loc='upper left')
savefig('ceph2.pdf')
show()