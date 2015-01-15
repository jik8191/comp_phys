#!/usr/bin/python
from pylab import *
import matplotlib.pyplot as plt

data= loadtxt("poissontest.dat", 'float')
plot(data[:,0], data[:,1], 'k-', label='$\rho (r)$')
xlabel("$r$")
ylabel("$S(r)$")
savefig('Sofr_test.png')
#savefig('numint.pdf')
show()

plot(data[:,0], data[:,2], 'k-', label='$v(r)$')
xlabel("$r$")
ylabel("$v(r)$")
savefig('vofr_test.png')
#savefig('numint.pdf')
show()


