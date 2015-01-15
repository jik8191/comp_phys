#!/usr/bin/python
from pylab import *
import matplotlib.pyplot as plt

data= loadtxt("hf1d_he.dat", 'float')
plot(data[:,0], data[:,1], 'b-', label='$R(r)$')
xlabel("$r$")
ylabel("$R(r)$")
rc('text', usetex=True)
rc('font', family='serif')
savefig('../../img/Rinitofr_he.png')
#savefig('numint.pdf')
show()
