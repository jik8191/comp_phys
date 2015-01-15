/* AEP 4380 Final Project

   Calculation of the ground state energy of helium
   using a self-consistent loop.

   Austin Liu
   ayl42
   2014-12-15

   run on a core i3 with g++ 4.7.2 on the 
   Debian 7.6 Linux distribution.
*/


#include <cstdlib>
#include <cmath>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

const double pi = 3.14159265358979;
const double a = 0.5299; // Bohr radius in Angstroms a = hbar^2/(me^2)
const double hbsq_m = 7.6359; // hbar^2/m in eV-Angstroms^2
const double esq = 14.409; // e^2 in eV-Angstroms
const double Z = 2.00;
const int N = 2001;
const double ri = 0.0;
const double rf = 4.0; // in Angstroms
const double h = (rf-ri)/(N-1);
int i;
int numiter; // counter
double r [N];
double R [N]; // for R(r)
double u [N]; // for u(r) = r*R(r)
double usq [N];
double unew [N];
double unewsq [N];
double w [N]; // to be used as a weighing function for Numerov's method
double dudr [N];
double dudrsq [N];
double eptmp [N]; // for holding the potential energy intergrand
double V [N]; // for the potential energy V(r)
double v [N]; // v(r) = rV(r) - easier to solve for this first.
// double rho [N]; // for the charge density distribution
double f [N]; // for the function on the right hand side of
              // Poisson's equation,  i.e. -4*pi*r*(rho(r))
double blank [N]; // to handle optional arguments into Numerov
double Ek;
double Ep;
double E;
double Eold;
  
/* A trial radial part of the wavefunction
   If R(r) = u(r)/((sqrt(4pi)r) then this is u(r) */
double uradial(double r){
  double t = (Z/a);
  return 2*sqrt(t)*t*r*exp(-t*r);  
}

/* A trial radial part of the wavefunction
   If R(r) = u(r)/((sqrt(4pi)r) then this is R(r) */
double Rradial(double r){
  double t = (Z/a);
  return sqrt(t/pi)*t*exp(-t*r);  
}

/* Initialize values in r and w */
void initialize(){
  for(i=0; i<N; i++){
    r[i]= h*i;
    u[i]= uradial(r[i]);
    R[i]= Rradial(r[i]);
    usq[i]= u[i]*u[i];
    blank[i]= 0;
  }  
}

/* Update usq after determining a new u(r) */
void updateusq(){
  for(i=0; i<N; i++){
    usq[i]= u[i]*u[i];
  }  
}

/* Computes integral using the trapezoid rule */
double trapint(double f [N]){
  double sum;
  sum= f[0]/2;
  for(i=1; i<N-1; i++){
    sum+= f[i];
  }
  sum += f[N-1]/2;
  sum *= h;
  return sum;
}

/* Compute integral using Simpson's rule, N must be odd */
double simpint(double f [N]){
  double sum;
  sum= f[0];
  for(i=1; i<N-1; i+=2){
    sum+= 4*f[i];
  }
  for(i=2; i<N-2; i+=2){
    sum+= 2*f[i];
  }
  sum += f[N-1];
  sum *= h/3;
  return sum;
}

/* Updates the derivative of a function using numerical differentiation
   Assumes f is a function of r. */
void updatederiv(double f [N], double (&dfdr) [N]){  
  // have to use forward difference
  dfdr[0]= (f[1]-f[0])/h;
  // use central difference for middle ones
  for (i=1; i< N-1; i++){
    dfdr[i]= 0.5*(f[i+1]-f[i-1])/h;
  }
  // use backward difference
  dfdr[N-1]= (f[N-2]-f[N-1])/h;
}

/* Compute the kinetic energy of a wavefunction */
double calcEk(double u [N]){
  double Ek;
  updatederiv(u, dudr);
  for (i=1; i<N; i++){
    dudrsq[i]= dudr[i]*dudr[i];  
  }
  Ek = hbsq_m*0.5*simpint(dudrsq);
  return Ek;
}

/* Compute the potential energy of a wavefunction using a 
   predetermined potential*/
double calcEp(double u [N]){
  for (i=1; i<N; i++){
    eptmp[i]= (-Z*esq/r[i] + 0.5*V[i])*usq[i];
  }
  return simpint(eptmp);
}

/* Use the current wavefunction to set the right hand side of 
   Poisson's equation. */
void setf(){
  updateusq();
  f[0]= 0.0; // boundary condition
  for (i=1; i<N; i++){
    f[i]= -esq*usq[i]/r[i];
  }
}


/* Uses the Numerov method to update the (k+1)th term of u
   k must be in range (2, N-2). Assumes terms
   u[0] through u[k] are already initialized and
   w is initialized. */
void numerovf(int k, double (&u) [N], double w [N], double s [N]){
  double b = h*h/12.0;
  u[k+1]= ( 2*(1-5*b*w[k])*u[k] - (1+ b*w[k-1])*u[k-1]
	    + b*(s[k+1]+10*s[k]+s[k-1]) )/(1+b*w[k+1]);
}

/* Uses the Numerov method to update the (k-1)th term of u
   k must be in range (2, N-2). Assumes terms
   u[0] through u[k] are already initialized and
   w is initialized. */
void numerovb(int k, double (&u) [N], double w [N], double s [N]){
  double b = h*h/12.0;
  u[k-1]= ( 2*(1-5*b*w[k])*u[k] - (1+ b*w[k+1])*u[k+1]
	    + b*(s[k+1]+10*s[k]+s[k-1]) )/(1+b*w[k-1]);
}

/* Solves Poisson's equation del^2 v = f*/
void calcv(){
  // initial conditions
  v[0]= 0.0; // require v[0] = 0 since v(r)/r = V(r)
  v[1]= 0.0; // v[1] is roughly 0.0
  for (i=1; i<N-1; i++){
    numerovf(i, v, blank, f);
  }
  // Need correction using Konnin's method p.54
  double m = (v[N-1]-v[N-11])/(h*10);
  // cout << "m = " << m << endl;
  double b = v[N-1]-m*(N-1)*h;
  // cout << "b = " << b << endl;
  for (i=0; i<N; i++){
    v[i]= v[i]- m*h*i;
  }
}

/* Calculate V(r) for use in the potential energy calculation. */
void calcV(){
  for(i=1; i<N; i++){
    V[i]= v[i]/r[i];
  }
  // linear approximation for V(r=0)
  V[0]= (V[1]-V[2])+V[1];
}

// write u(r) output to a file
void writeutofile(){
  stringstream fname;
  // cout << t << endl;
  fname << "data/u_" << numiter << ".dat";
  ofstream fp;
  fp.open(fname.str().c_str());
  if ( fp.fail() ){
    cout << "cannot open file" << endl;
  }
  for (i=0; i<N; i++){
    fp << r[i] << setw(15) << u[i] << setw(15) << endl;
    // << w[i] << setw(15)
    // << dudr[i] << setw(15) << R[i] << setw(15)
    // << f[i] << setw(15) << v[i] << 
      
  }
  fp.close();
}

// write R(r) output to file last
void writeRtofile(){
  ofstream fp;
  fp.open("hf1d_he.dat");
  if ( fp.fail() ){
    cout << "cannot open file" << endl;
  }
  for (i=0; i<N; i++){
    fp << r[i] << setw(15) << R[i] << setw(15) << endl;
    // << w[i] << setw(15)
    // << dudr[i] << setw(15) << R[i] << setw(15)
    // << f[i] << setw(15) << v[i] <<       
  }
  fp.close();
}


/* Determine the weighing function for use in calculating
   unew using Numerov's method. */
void calcw(double E){
  for(i=1; i<N; i++){
    w[i]= (2/hbsq_m)*(Z*esq/r[i] - 0.5*V[i] + E);
  }
  w[0]= 2*w[1]-w[2]; // linear extrapolation
}

/* Given a value of the total energy, use Numerov's method
   to trace a new wavefunction and normalize it accordingly. */
void calcunew(){
  unew[N-1]= 0.1*u[N-1];
  unew[N-2]= 0.1*u[N-2]; // initial conditions
  // for (i=1; i<N/10; i++){
  //   numerovf(i, unew, w, blank);
  // }
  for (i=N-2; i>1; i--){
    numerovb(i, unew, w, blank);
  }
  // normalization required
  for (i=0; i<N; i++){
    unewsq[i]= unew[i]*unew[i];
  }
  double m = sqrt(simpint(unewsq));
  // cout << "norm. const. = " << m << endl;
  for (i=0; i<N; i++){
    unew[i]/= m;
    unewsq[i]= unew[i]*unew[i];
  }
  // cout << "normalized to: " << simpint(unewsq) << endl;
}

void updateR(){
  for (i=1; i<N; i++){
    R[i]= u[i]/(sqrt(4*pi)*r[i]);
   }
}

void scfiter(){
  setf(); // del^2 v = f (Poisson's equation)
  calcv();
  calcV();
  Ep= calcEp(u);
  Ek= calcEk(u);
  E= Ep+Ek;
  cout << "kinetic energy is: " << Ek << endl;
  cout << "potential energy is: " << Ep << endl;
  cout << "total energy is: " << E << endl;
  calcw(E); // weighing function to solve for new u
  calcunew();
  // update old values for comparison
  for (i=0; i<N; i++){
    u[i]= unew[i];
  }
  numiter++;
  cout << "iteration " << numiter << " completed." << endl;
  writeutofile();
}

int main(){
  // normalization test
  initialize();
  writeutofile(); // write initial u
  Eold= -50.0; // start E at some value
  E= 0.0;
  while (abs(E-Eold) > 0.001 && numiter < 15){
    Eold= E;
    scfiter();
  }
  Ep= calcEp(u);
  Ek= calcEk(u);
  E= Ep+Ek;
  cout << "final kinetic energy is: " << Ek << endl;
  cout << "final potential energy is: " << Ep << endl;
  cout << "Ek/Ep is: " << Ek/Ep << endl;
  cout << "final total energy is: " << E << endl;
  cout << "stepsize was: " << h << endl;
  cout << "cutoff radius was: " << rf << endl;
  writeRtofile();  
  return(EXIT_SUCCESS);
} 


