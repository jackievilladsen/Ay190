# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 14:50:08 2012

@author: jackie
"""

from copy import deepcopy
from pylab import *
from time import *
    
######### code from online: http://snippets.dzone.com/posts/show/5645
    
# this function, swapRows, was adapted from
# Numerical Methods Engineering with Python, Jean Kiusalaas
def swapRows(v,i,j):
	"""Swaps rows i and j of vector or matrix [v]."""
	if len(v) == 1:
		v[i],v[j] = v[j],v[i]
	else:
		temp = v[i].copy()
		v[i] = v[j]
		v[j] = temp
	
def pivoting(a, b):
	"""changes matrix A by pivoting"""
	
	n = len(b)
	
	for k in range(0, n-1):
		p = int(argmax(abs(a[k:n, k]))) + k
		if (p != k):
			swapRows(b, k, p)
			swapRows(a,k,p)

def gauss(a, b, t=1.0e-9, verbose=False):
	""" Solves [a|b] by gauss elimination"""
	
	n = len(b)
	
	# make copies of a and b so as not to change the values in the arguments
	tempa = deepcopy(a)
	tempb = deepcopy(b)
	
	# check if matrix is singular
	if abs(linalg.det(tempa)) < t:
		print "asn"
		return -1
	
	pivoting(tempa, tempb)

	for k in range(0,n-1):	
		for i in range(k+1, n):
			if tempa[i,k] != 0.0:
				m = tempa[i,k]/tempa[k,k]
				if verbose:
				    print "m =", m
				tempa[i,k+1:n] = tempa[i,k+1:n] - m * tempa[k,k+1:n]
				tempb[i] = tempb[i] - m * tempb[k]
	
	# Back substitution
	for k in range(n-1,-1,-1):
		tempb[k] = (tempb[k] - dot(tempa[k,k+1:n], tempb[k+1:n]))/tempa[k,k]

	return tempb
 
############### end code from online, begin my code:
    
def load_Ab(i):
    # loads matrix A and b for system linear equations Ax = b
    # i = 1 through 5
    
    Afile = 'LSE' + repr(i) + '_m.dat'
    bfile = 'LSE' + repr(i) + '_bvec.dat'
    
    A = loadtxt(Afile)
    b = loadtxt(bfile)

    ndims=len(b)    
    
    print 'number of equations in system', i, ':', ndims  
    print 'det(A):', det(A)
    
    return A,b

def test_gauss():
    # to test the Gaussian elimination function I got online
    A=matrix([[0.,1.,0.],[1.,0.,0.],[0.,0.,2.]])
    b=array([2.,1.,3.])
    x=gauss(A,b)
    print x # yields [1. 2. 1.5], the correct answer

def compare_times(A,b):
    # returns times to get solutions for Gaussian elimination
    # vs. built-in solver (which uses LU)
    t1=time()
    x=solve(A,b)
    t2=time()
    x2=gauss(A,b)
    t3=time()
    t_builtin = t2-t1
    t_Gauss = t3-t2
    print 'For system', i, 'Gaussian elimination took', t_Gauss, 'seconds and the built-in solver took',t_builtin, 'seconds'
    return t_Gauss, t_builtin

#test_gauss()
irange=range(1,6)
tG = zeros(len(irange))
tBI = zeros(len(irange))
n_eqns = zeros(len(irange))

for i in irange:
    [A,b]=load_Ab(i)
    n_eqns[i-1]=len(b)
    [tG[i-1],tBI[i-1]]=compare_times(A,b)
    
figure    
loglog(n_eqns,tG,'-*')
hold(True)
loglog(n_eqns,tBI,'-*')
xlabel('Number of Equations')
ylabel('Time (s)')
title('Efficiency of Methods for Solving Systems of Linear Equations')
legend(['Gauss Elimination','Built-In Solver'],loc='upper left')
savefig('linsystimes.pdf')
show()