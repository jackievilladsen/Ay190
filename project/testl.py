# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 14:29:12 2012

@author: jackie
"""

from matplotlib import *
from pylab import *

def testl(mfp):
    y=rand()
    l=-1. * mfp * log(1-y) # step size
    return l

mfp = 1
ntest=100000
l_array=zeros(ntest)
for i in range(ntest):
    l_array[i]=testl(mfp)

figure()
hist(l_array,bins=50,histtype='step',log=True)
print 'mean:', mean(l_array)
print 'mean l^2:', mean(l_array**2)
print 'rms l:', sqrt(mean(l_array**2))
xlabel('Step size / Mean free path')
ylabel('# Occurrences out of 100,000')
title('Test of Photon Step Size Monte Carlo')
savefig('writeup/mctest.pdf')