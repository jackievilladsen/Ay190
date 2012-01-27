#!/opt/local/bin/python

import sys
from scipy import *
from pylab import *

# everything is in cgs

# global constants
ggrav = 6.67e-8
clite = 3.0e10
msun = 1.99e33


# EOS
# neutron stars:
# polyG = 2.0
# polyK = 100.0 * 5.55e38/6.1755e17**polyG

# EOS for 
# white dwarfs:
polyG = 4.0/3.0
polyK = 1.244e15*0.5**polyG

# central values
rhoc = 1.0e10

# polytropic relation between rho and P
def P_of_rho(rho):
    return polyK * rho ** polyG
def rho_of_P(P):
    return (P/polyK)**(1/polyG)

# minimum pressure
rhomin = 1.e-10 * rhoc
min_press = P_of_rho(rhomin)

# define internal energy
def int_energy(P,rho):
    return P / ((polyG-1) * rho)


# grid
rmax = 1.0e9


def set_grid(rmax,nzones):
    # set up the grid and return the
    # radius array and dr
    
    # uniform linear grid

    rad = linspace(0.,rmax,nzones) # max r is actually rmax-dr
    dr = rad[1]-rad[0]
    
    return (rad,dr)

def tov_RHS(r,data):
    # evaluates the right-hand side of the stellar structure equations
    # rhs[0] = dP/dr = -G * M(r) * rho(r) / r^2
    # rhs[1] = dM/dr = 4*pi*r^2 * rho(r)

    rhs = zeros(2)

    mass = data[1]
    press = max(data[0],min_press)

    rho = rho_of_P(press)

    # hydrostatic equilibrium
    if abs(r)<=1e-10:  # i.e. at r=0
        rhs[0] = 0     # lim (r->0) dP/dr = 0
    else:
        rhs[0] = - ggrav * mass * rho / r**2
        
    # mass conservation
    rhs[1] = 4*pi*r**2 * rho

    return rhs

def tov_RK1(old_data,r,dr):
    # uses first-order Runge-Kutta (forward Euler)
    # to step forward one step in the integration
    k = dr * tov_RHS(r,old_data)
    new_data = old_data + k
    return new_data

def tov_RK2(old_data,r,dr):
    # uses 2nd-order Runge-Kutta to step forward one step
    # in the integration
    k1 = dr * tov_RHS(r,old_data)
    k2 = dr * tov_RHS(r+dr/2,old_data+k1/2)    
    new_data = old_data + k2
    return new_data
    
def tov_RK3(old_data,r,dr):
    # uses 3rd-order Runge-Kutta to step forward one step
    # in the integration
    k1 = dr * tov_RHS(r,old_data)
    k2 = dr * tov_RHS(r+dr/2,old_data+k1/2)
    k3 = dr * tov_RHS(r+dr,old_data-k1+2*k2)    
    new_data = old_data + (k1 + 4*k2 + k3)/6
    return new_data

def tov_RK4(old_data,r,dr):
    # uses 4th-order Runge-Kutta to step forward one step
    # in the integration
    k1 = dr * tov_RHS(r,old_data)
    k2 = dr * tov_RHS(r+dr/2,old_data+k1/2)
    k3 = dr * tov_RHS(r+dr/2,old_data+k2/2)
    k4 = dr * tov_RHS(r+dr,old_data+k3)
    new_data = old_data + (k1 + 2*k2 + 2*k3 + k4)/6
    return new_data

def tov_integrate(rmax,nzones,RKn):

    # set up grid
    (rad,dr) = set_grid(rmax,nzones)

    # initialize some variables
    tovdata = zeros((nzones,2))
    # 0 -- press
    # 1 -- mbary

    tovout = zeros((nzones,4))
    # 0 -- rho
    # 1 -- press
    # 2 -- eps
    # 3 -- mass
     
    # central values
    Pc = P_of_rho(rhoc)
    tovdata[0,0] = Pc
    tovdata[0,1] = 0.0
    
    tovout[0,0] = rhoc
    tovout[0,1] = Pc
    tovout[0,2] = int_energy(Pc,rhoc)
    tovout[0,3] = 0.0


    # you will need to track the surface (where press <= press_min)
    isurf = 0
    for i in range(nzones-1):

        # integrate one step using RK1 (can change this to
        # RK2 or RK3 or RK4)

        if RKn==1:
            tovdata[i+1,:] = tov_RK1(tovdata[i,:],rad[i],dr)

        elif RKn==2:
            tovdata[i+1,:] = tov_RK2(tovdata[i,:],rad[i],dr)

        elif RKn==3:
            tovdata[i+1,:] = tov_RK3(tovdata[i,:],rad[i],dr)

        elif RKn==4:
            tovdata[i+1,:] = tov_RK4(tovdata[i,:],rad[i],dr)

        # check if press below 0
        if(tovdata[i+1,0] <= min_press):
            if isurf==0:
                isurf = i
            tovdata[i+1,0] = min_press # prevent negative pressures

        # add pressure to output array tovout
        P = tovdata[i+1,0]
        tovout[i+1,1] = P

        # add mass to output array tovout
        if (i+1 > isurf and isurf > 0):
            tovout[i+1,3] = tovdata[isurf,1]
        else:
            tovout[i+1,3] = tovdata[i+1,1]
            
        # compute density and add to output array tovout
        rho = rho_of_P(P)        
        tovout[i+1,0] = rho
        
        # compute eps (internal energy) and add to output array tovout
        tovout[i+1,2] = int_energy(P,rho)

    return (tovout,isurf,dr)



# for convergence: 
# number of points
na = array([10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000])
# to store masses
masses = zeros((len(na),4))
# to store the drs
radii = zeros((len(na),4))

for i in range(len(na)):
    print na[i]
    for j in range(4):
        print j+1
        (tov_star,isurf,dr) = tov_integrate(rmax,na[i],j+1)
    
        # tov_star (n x 4 array):
            # 0 -- rho
            # 1 -- press
            # 2 -- eps
            # 3 -- mass
           
        # record mass and dr for this number of points    
        masses[i,j] = tov_star[isurf,3]
#        print masses[i,j]/msun # mass in M_Sun
        radii[i,j] = dr*isurf

# stuff for plotting radial profile of stellar properties
#    # save stuff into variables with names I understand, for easy display
#    rho = tov_star[0:isurf,0]    
#    P = tov_star[0:isurf,1]
#    eps = tov_star[0:isurf,2]
#    M = tov_star[0:isurf,3]
#    
#    # define an array of the radius for plotting purposes
#    (rad,dr) = set_grid(rmax,na[i])
#    rad=rad[0:isurf]
#    print rad[isurf-1]/(1.e5) # radius in km
#
#    # plot interesting quantities (normalized so max is 1)
#    plot(rad,P/P_of_rho(rhoc),'b-')
#    plot(rad,rho/rhoc,'g-')
#    plot(rad,M/masses[i],'r-')
#    plot(rad,eps/eps[0],'k-')
#    legend(('Pressure','Density','Enclosed Mass','Internal Energy'),loc='center right')
#    xlabel('Radius (cm)')
#    ylabel('Quantity/Max(Quantity)')
#    title('Radial Profile of White Dwarf')
#    savefig('radial.pdf')    
#    show()

l=len(na)
ans = masses[l-1,3]
frac_err = abs((masses-ans)/ans)

figure
hold(True)
loglog(na,frac_err[:,0],'-*')
loglog(na,frac_err[:,1],'-*')
loglog(na,frac_err[:,2],'-*')
loglog(na,frac_err[:,3],'-*')
legend(('RK1','RK2','RK3','RK4'))
xlabel('Number of Integration Steps')
ylabel('Fractional Error')
title('Convergence of Runge-Kutta Integration')
savefig('rk.pdf')
show()

na_ratio = na[1:l-1]/na[0:l-2]
err_ratio = frac_err[1:l-1]/frac_err[0:l-2]
convergence_power1 = log(err_ratio[:,0])/log(na_ratio)
print convergence_power1

convergence_power2 = log(err_ratio[:,1])/log(na_ratio)
print convergence_power2

convergence_power3 = log(err_ratio[:,2])/log(na_ratio)
print convergence_power3

convergence_power4 = log(err_ratio[:,3])/log(na_ratio)
print convergence_power4
