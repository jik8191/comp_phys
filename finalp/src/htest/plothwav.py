#!/usr/bin/python
from pylab import *
import matplotlib.pyplot as plt

data= loadtxt("hf1d_h.dat", 'float')
plot(data[:,0], data[:,1], 'b-', label='$R(r)$')
plot(data[:,0], data[:,2], 'r-', label='$R_{num}(r)$')
xlabel("$r$")
ylabel("$R(r)$")
rc('text', usetex=True)
rc('font', family='serif')
savefig('../../img/Rofr_h.png')
legend(loc=4)
#savefig('numint.pdf')
show()
