#!/usr/bin/python
from pylab import *
import matplotlib.pyplot as plt

datai= loadtxt("heinit.dat", 'float')
plot(datai[:,0], datai[:,1], 'b-', label='$R(r)$')
xlabel("$r$")
ylabel("$R(r)$")
rc('text', usetex=True)
rc('font', family='serif')
savefig('../img/Rinitofr_he.png')
#savefig('numint.pdf')
show()

numiter= '6'
data0= loadtxt("data/u_0.dat", 'float')
data= loadtxt("data/u_" + numiter + ".dat", 'float')
plot(data0[:,0], data0[:,1], 'b-', label='$u_i(r)$')
plot(data[:,0], data[:,1], 'k-', label='$u_f(r)$')
xlabel("$r$")
ylabel("$u_(r)$")
xlim(0.0,3.0)
rc('text', usetex=True)
rc('font', family='serif')
legend(loc=4)
savefig('../img/uofr_he.png')
show()


data= loadtxt("hf1d_he.dat", 'float')
plot(datai[:,0], datai[:,1], 'b-', label='$R_i(r)$')
plot(data[:,0], data[:,1], 'k-', label='$R_f(r)$')
xlabel("$r$")
ylabel("$R(r)$")
xlim(0.0,3.0)
rc('text', usetex=True)
rc('font', family='serif')
legend(loc=4)
savefig('../img/Rofr_he.png')
show()
