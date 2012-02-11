import sys,math
from pylab import *

sigma=sqrt(15.)
x0=30.
v=0.1

def apply_bcs(x,y):
    # apply boundary conditions    
    l=len(y)
    y[0]=y[1]
    y[l-1]=y[l-2]
    return y

def step_upwind(yold,yolder,v,dx,dt):
    # implements upwind method to step forward by dt
    # where yold is the array of values of y at all x,
    # returns ynew, the array of new y values at all x
    
    l=len(yold)
    ydiff = zeros(l)
    ydiff[1:l-1] = yold[1:l-1]-yold[0:l-2]
    ynew = yold - v*dt/dx * ydiff
    return ynew
    
def step_downwind(yold,yolder,v,dx,dt):
    # implements downwind method to step forward by dt
    # where yold is the array of values of y at all x,
    # returns ynew, the array of new y values at all x
    
    l=len(yold)
    ydiff = zeros(l)
    ydiff[0:l-2] = yold[1:l-1]-yold[0:l-2]
    ynew = yold - v*dt/dx * ydiff
    return ynew
    
def step_FTCS(yold,yolder,v,dx,dt):
    # implements FTCS method to step forward by dt
    # where yold is the array of values of y at all x,
    # returns ynew, the array of new y values at all x
    
    l=len(yold)
    ydiff = zeros(l)
    ydiff[1:l-2] = yold[2:l-1]-yold[0:l-3]
    ynew = yold - v*dt/dx * ydiff
    return ynew
    
def step_LaxFriedrich(yold,yolder,v,dx,dt):
    # implements Lax-Friedrich method to step forward by dt
    # where yold is the array of values of y at all x,
    # returns ynew, the array of new y values at all x
    
    l=len(yold)
    ydiff = zeros(l)
    ydiff[1:l-2] = yold[2:l-1]-yold[0:l-3]
    ymean=zeros(l)
    ymean[1:l-2] = (yold[2:l-1]+yold[0:l-3])/2
    ynew = ymean - v*dt/dx * ydiff/2
    return ynew
    
def step_leapfrog(yold,yolder,v,dx,dt):
    # implements leapfrog method to step forward by dt
    # where yold is the array of values of y at all x,
    # returns ynew, the array of new y values at all x
    
    l=len(yold)
    ydiff = zeros(l)
    ydiff[1:l-2] = yold[2:l-1]-yold[0:l-3]
    ynew = yolder - v*dt/dx * ydiff
    return ynew

def step_LaxWendroff(yold,yolder,v,dx,dt):
    # implements Lax-Wendroff method to step forward by dt
    # where yold is the array of values of y at all x,
    # returns ynew, the array of new y values at all x
    
    l=len(yold)
    ydiff = zeros(l)
    ydiff[1:l-2] = yold[2:l-1]-yold[0:l-3]
    ymean=zeros(l)
    ymean[1:l-2] = (yold[2:l-1]+yold[0:l-3])/2
    ynew = yold - v*dt/dx * ydiff/2 + (v*dt/dx)**2 * (ymean - yold)
    return ynew

def analytic(t,x):
    # returns y(x) evaluated at time t
    # if x is an array, then y is an array too
    
    y=normpdf(x-v*t,x0,sigma)
    return y

def makeGrid(x0,xf,dx):
    # makes grid from x0 to xf (inclusive) with spacing dx
    
    n_points = (xf-x0)/dx + 1
    xgrid=linspace(x0,xf,n_points)
    return xgrid

def get_err_vs_t(step_function,dx,alpha):
    # parameters
    v = 0.1
    
    # set up the grid
    xmin=0
    xmax=100
    x = makeGrid(xmin,xmax,dx)
    n = len(x)
    y = zeros(n)
    cfl = alpha
    dt = cfl * dx/v # satisfy upwind CFL stability condition
    ntmax = int((xmax-10.-xmin)/(v*dt))
    #ntmax=5
    t = 0.0
    
    #set up initial conditions
    y = analytic(t,x)
    
    # evolve (and show evolution)
    ion()
    figure()
    plot(x,y,'x-') # numerical data
    plot(x,y,'r-') # analytic data
    show()
    maxy=max(y)*1.1
    miny=min(y)-0.1*max(y)
    maxx=100
    minx=0
    axis([minx,maxx,miny,maxy])
    
    yold2 = y
    yold = y
    err_vs_t = zeros(ntmax)
    for it in range(ntmax):
        t = t+dt
        # save previous and previous previous data
        yold2 = yold
        yold = y

        # get new data
        y = step_function(yold,yold2,v,dx,dt)
    
        # after update, apply boundary conditions
        y = apply_bcs(x,y) 
    
        # get analytic result for time t
        yana = analytic(t,x)
        # compute error estimate
        err = sum(abs(y-yana))/sum(yana)
        print "it = ",it, err
        err_vs_t[it]=err
        clf()
        # plot numerical result
        plot(x,y,'x-')
        # plot analytic results
        plot(x,yana,'r-')    
        axis([minx,maxx,miny,maxy])
        draw()
    #savefig('laxfried' + repr(ntmax) + '.pdf')
    show()
    t=linspace(dt,dt*ntmax,ntmax)
    return t, err_vs_t



dx=1. # recommended: 0.1
alpha=0.9
[t,errlw]=get_err_vs_t(step_LaxWendroff,dx,alpha)
[t2,errlw2]=get_err_vs_t(step_LaxWendroff,dx/2,alpha)
[tleap,errleap]=get_err_vs_t(step_leapfrog,dx,alpha)

figure()
hold(True)
semilogy(t,errlw)
semilogy(tleap,errleap)
xlabel('Time')
ylabel('Error')
title('Error for alpha=0.5, dx=1')
legend(('Lax-Wendroff','Leapfrog'),loc='lower right')
savefig('err_lw_vs_leap.pdf')
show()

conv=abs((log(errlw2[range(0,200,2)])-log(errlw))/(log(dx)-log(dx/2)))
figure()
plot(t,conv)
xlabel('Time')
ylabel('Calculated Convergence Rate')
title('Convergence Rate of Lax-Wendroff Method')
savefig('conv.pdf')
show()