# -*- coding: utf-8 -*-
"""
Created on Sat Jan 28 16:38:35 2012

@author: jackie
"""

'''general Runge-Kutta stuff'''

from pylab import *

def RK1(old_data,x,dx,RHS):
    # uses first-order Runge-Kutta (forward Euler)
    # to step forward one step in the integration
    k = dx * RHS(x,old_data)
    new_data = old_data + k
    return new_data

def RK2(old_data,x,dx,RHS):
    # uses 2nd-order Runge-Kutta to step forward one step
    # in the integration
    k1 = dx * RHS(x,old_data)
    k2 = dx * RHS(x+dx/2,old_data+k1/2)    
    new_data = old_data + k2
    return new_data
    
def RK4(old_data,x,dx,RHS):
    # uses 4th-order Runge-Kutta to step forward one step
    # in the integration
    k1 = dx * RHS(x,old_data)
    k2 = dx * RHS(x+dx/2,old_data+k1/2)
    k3 = dx * RHS(x+dx/2,old_data+k2/2)
    k4 = dx * RHS(x+dx,old_data+k3)
    new_data = old_data + (k1 + 2*k2 + 2*k3 + k4)/6
    return new_data
    
def backEuler(old_data,dx):
    # uses backwards Euler to step forward one step in the
    # integration
    
    # note: this function is specific to this problem,
    # but other integrating functions aren't
    
    A=matrix([[1.,-99.],[-1.,99.]])
    thingToInvert = eye(2)-dx*A
    multiplyBy = asarray(thingToInvert.I)
    new_data = dot(multiplyBy,old_data)
    return new_data
    
def makeGrid(x0,xf,dx):
    # makes grid from x0 to xf (inclusive) with spacing dx
    
    n_points = (xf-x0)/dx + 1
    xgrid=linspace(x0,xf,n_points)
    return xgrid

def RKintegrate(x0,xf,dx,y0,dy_dx,order):
    # uses RK of order "order" with step size dx to integrate
    # dy_dx from x0 to xf, with initial condition y(x0)=f0
    
    xgrid=makeGrid(x0,xf,dx)
    
    def RK(old_data,x,dx,RHS):
        if order==1:
            return RK1(old_data,x,dx,RHS)
        if order==2:
            return RK2(old_data,x,dx,RHS)
        if order==4:
            return RK4(old_data,x,dx,RHS)
        if order=='1b':
            return backEuler(old_data,dx)
    
    # store initial condition in grid of solution
    ygrid=zeros((len(xgrid),len(y0)))
    ygrid[0]=y0
    
    # iterate to fill in the rest of the gridded solution
    nsteps = len(xgrid)-1
    y=y0
    for i in range(nsteps):
        new_y=RK(y,xgrid[i],dx,dy_dx)
        ygrid[i+1,:]=new_y
        y=new_y
    
    return xgrid,ygrid

'''specific to this problem'''

def dY_dt(t,Y):
    # RHS function for the system of equations:
    # dY1/dt = Y1 - 99 Y2
    # dY2/dt = -Y1 + 99 Y2
    # where the input Y is the vector [Y1,Y2]
    
    A=[[1,-99],[-1,99]]
    return dot(A,Y)

def Yanalytic(t):
    # returns the analytic solution to the system of eqns at time t
    Y1 = (exp(100*t) + 99.)/100.
    Y2 = (-exp(100*t) + 1.)/100.
    return Y1,Y2

t0=0.
tf=4.
Y0=[1,0]

#dt=0.001
#tgrid=makeGrid(t0,tf,dt)
#[Y1analytic,Y2analytic]=Yanalytic(tgrid)
[Y1analytic,Y2analytic]=Yanalytic(4.)

#figure
#hold(True)
#subplot(211)
#semilogy(tgrid,Y1analytic)
#subplot(212)
#semilogy(tgrid,-Y2analytic)

order_range=[1,2,4,'1b']
#order_range=['1b']
nsteps=floor(10**arange(2.5,6,0.5))
dt_grid=(tf-t0)/nsteps

n_orders = len(order_range)
n_stepsizes = len(dt_grid)

frac_err=zeros((n_stepsizes,n_orders))
Y1final=zeros((n_stepsizes,n_orders))
Y2final=zeros((n_stepsizes,n_orders))

for j in range(n_orders):
    for k in range(n_stepsizes):
        dt=dt_grid[k]
        [tgrid,Ygrid]=RKintegrate(t0,tf,dt,Y0,dY_dt,order_range[j])
        Y1=Ygrid[:,0]
        Y2=Ygrid[:,1]
        l=len(Y1)
        Y1final[k,j]=Y1[l-1]
        Y2final[k,j]=Y2[l-1]
        print j,k
    
#    subplot(211)
#    semilogy(tgrid,Y1)
#    subplot(212)
#    semilogy(tgrid,-Y2)
#    legend(('analytic','numerical'))
#    
#    frac_err1 = abs((Y1-Y1analytic)/Y1)
#    frac_err2 = abs((Y2-Y2analytic)/Y2)
#    l=len(frac_err1)
#    print 'order:', order_range[j]
#    print 'fractional error on Y1(t=4):', frac_err1[l-1]
#    print 'fractional error on Y2(t=4):', frac_err2[l-1]
#    frac_err[j]=frac_err1[l-1]

frac_err = abs((Y1final-Y1analytic)/Y1final)
frac_err_v2 = abs((Y1final-Y1analytic)/Y1analytic)
figure
subplot(211)
hold(True)
loglog(dt_grid,frac_err[:,0])
loglog(dt_grid,frac_err[:,1])
loglog(dt_grid,frac_err[:,2])
legend(('Euler','RK2','RK4'))
ylabel('Fractional Error')
title('Convergence of Stiff ODE Integrators')
subplot(212)
loglog(dt_grid,frac_err_v2[:,3])
xlabel('Step Size dt')
ylabel('Fractional Error')
title('Backwards Euler')
savefig('stiff.pdf')
show()