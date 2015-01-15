/* AEP 4380 Final Project

   Testing solution of Poisson's equation
   using Numerov's method.
   Austin Liu
   ayl42
   2014-10-30
*/

#include <cstdlib>
#include <cmath>

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

const double pi= 3.14159265358979;
const double a = 0.5299; // Bohr radius in Angstroms a = hbar^2/(me^2)
const double ri = 0.0;
const double rf = 4.0;
const int N = 1001;
const double h = (rf-ri)/(N-1);
const double esq = 14.409; // e^2 in eV-Angstroms
double rho [N]; // for the charge density distribution
double S [N]; // for Numerov
double r [N];
double v [N];
double V [N];
int i;

double rhoofr(double r){
  return (1/(pi*a*a*a))*exp(-2*r/a);
}

void initialize(){
  for(i=0; i<N; i++){
    r[i]= h*i;
    rho[i]= rhoofr(r[i]);
    S[i]= -4*pi*esq*r[i]*rho[i];
  }  
}

void set_ics(){
  v[0]= 0;
  v[1]= 0;
  //v[1]= 1-0.5*(h+2)*exp(-h); // test exact solution
}


/* Uses the Numerov method to update the (k+1)th term of u
   k must be in range (2, N-2). Assumes terms
   u[0] through u[k] are already initialized. */
void numerovf(int k, double (&u) [N], double s [N]){
  double b = h*h/12.0;
  u[k+1]= 2*u[k] - u[k-1]+ b*(s[k+1]+10*s[k]+s[k-1]);
}

void calcv(){
  for (i=1;i<N-1;i++){
    numerovf(i, v, S);
  }
  // testing correction using Konnin's method
  double m = (v[N-1]-v[N-11])/(h*10);
  cout << "m = " << m << endl;
  double b = v[N-1]-m*(N-1)*h;
  cout << "b = " << b << endl;
  for (i=0; i<N; i++){
    v[i]= v[i]- m*h*i;
  }
}

void writetofile(){
  // stringstream fname;
  // cout << t << endl;
  // fname << "wave2d_t_" << t*1000 << "_ms.dat";
  ofstream fp;
  // fp.open(fname.str().c_str());
  fp.open("poissontest.dat");
  if ( fp.fail() ){
    cout << "cannot open file" << endl;
  }
  for (i=0; i<N; i++){
    fp << r[i] << setw(15) << S[i] << setw(15) 
       << v[i] <<  endl;
  }
  fp.close();
}

int main(){
  initialize(); // initialize r, rho
  set_ics();
  calcv();
  writetofile();
  return(EXIT_SUCCESS);
}
