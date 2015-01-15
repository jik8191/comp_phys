/* AEP 4380 Final Project

   Test of Numerov's method and
   various methods to calculate kinetic and
   potential energy for a hydrogen 1s wavefunction.

   Austin Liu
   ayl42
   2014-12-15

   run on a core i3 with g++ 4.7.2 on the Debian 7.6 Linux distribution.
*/

#include <cstdlib>
#include <cmath>

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

const double pi= 3.14159265358979;
const int N = 10001;
const double a = 0.5299; // Bohr radius in Angstroms a = hbar^2/(me^2)
const double esq = 14.409; // e^2 in eV-Angstroms
const double hbsq_m = 7.6359; // hbar^2/m in eV-Angstroms^2
const double ri = 0.0;
const double rf = 3.0;
const double h = (rf-ri)/(N-1);
const double Z = 1.0;
int i;
double Ep, Ek, E;
double r [N];
double u [N];
double usq [N];
double dudr [N];
double dudrsq [N];
double w [N];
double R [N];
double Rnum [N];
double eptmp [N];
  
/* The hydrogenic 1s wavefunction */
double f(double r){
  double t = (Z/a);
  return sqrt(t/pi)*t*exp(-t*r);
}

/* Weighting function to solve for the
   wavefunction numerically, i.e. w(r) in the 
   expression for Numerov's method. */
double g(double r){
  return (2/hbsq_m)*((Z*esq/r) - 13.605);   
}

/* Sets the initial conditions*/
void set_ics(double (&u)[N]){
  u[0]= 0.0;
  u[1]= 2*sqrt(Z/a)*(Z/a)*h*exp(-Z*h/a) ; 
  // try exact values
  u[N-1]= 2*sqrt(Z/a)*(Z/a)*h*(N-1)*exp(-Z*h*(N-1)/a);
  u[N-2]= 2*sqrt(Z/a)*(Z/a)*h*(N-2)*exp(-Z*h*(N-2)/a);
}

/* Update usq after determining a new u(r) */
void updateusq(){
  for(i=0; i<N; i++){
    usq[i]= u[i]*u[i];
  }  
}

/* Uses the Numerov method to update the (k+1)th term 
   k must be in range (0, N-1). Assumes terms
   u[0] through u[k] are already initialized and
   w is initialized. */
void numerovf(int k, double (&u) [N], double w [N]){
  double a = h*h/12.0;
  u[k+1]= ( 2*u[k]-u[k-1]-a*(10*w[k]*u[k] + w[k-1]*u[k-1]) )/(1+a*w[k+1]);
}

/* Uses the Numerov method to update the (k-1)th term 
   k must be in range (0, N-1). Assumes terms
   u[N-1] through u[k] are already initialized and
   w is initialized. */
void numerovb(int k, double (&u) [N], double w [N]){
  double a = h*h/12.0;
  u[k-1]= ( 2*u[k]-u[k+1]-a*(10*w[k]*u[k] + w[k+1]*u[k+1]) )/(1+a*w[k-1]);
}

/* Initialize values in r and w */
void initialize(double (&r) [N], double (&w) [N]){
  for(i=0; i<N; i++){
    r[i]= h*i;
    w[i]= g(r[i]);
    R[i]= f(r[i]);
  }  
  // use linear extrapolation because w[0] goes to infinity
  // due to 1/r term in denominator.
  w[0]= 2*w[1]-w[2]; 
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
    eptmp[i]= (-Z*esq/r[i])*usq[i];
  }
  return simpint(eptmp);
}


void writetofile(){
  ofstream fp;
  fp.open("hf1d_h.dat");
  if (fp.fail()){
    cout << "cannot open file" << endl;
  }
  for (i=0; i<N; i++){
    fp << r[i] << setw(15) << R[i] << setw(15) 
       << Rnum[i] << endl;
  }
}

/* Solve for the hydrogenic wavefunction numerically */
void calcu(){
  // trace forward for first half
  for (i= 1; i<N/2; i++){
    numerovf(i, u, w);
  }
  // trace backward for second half
  for (i= N-2; i>N/2; i--){
    numerovb(i, u, w);
  }
  // normalization required
  for (i=0; i<N; i++){
    usq[i]= u[i]*u[i];
  }
  double m = sqrt(simpint(usq));
  cout << "norm. const. = " << m << endl;
  for (i=0; i<N; i++){
    u[i]/= m;
    usq[i]= u[i]*u[i];
  }
  cout << "normalized to: " << simpint(usq) << endl;
}

int main(){

  initialize(r,w);
  set_ics(u);

  // for (i= 1; i<N-1; i++){
  //   numerovf(i, u, w);
  // }
  calcu();
  updateusq();
  Ep= calcEp(u);
  Ek= calcEk(u);
  E= Ep+Ek;
  cout << "kinetic energy is: " << Ek << endl;
  cout << "potential energy is: " << Ep << endl;
  cout << "total energy is: " << E << endl;

  // generate Rnum values for comparison
  for (i= 1; i<N; i++){
    Rnum[i] = u[i]/(r[i]*sqrt(4*pi));
  }  
  Rnum[0] = 2*Rnum[1]-Rnum[2];

  writetofile();
  return (EXIT_SUCCESS);
} 
