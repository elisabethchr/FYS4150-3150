#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <cctype>
#include <string>
#include <math.h>

using namespace std;
using namespace arma;
ofstream ofile;


inline double f_array(double x){
    //vec f;
    //f = 100*exp(-10*x);
    return 100*exp(-10*x);
}
/*
vec exact(double x, int n) {
  vec fexact;
  fexact = 1.0-(1-exp(-10))*x - exp(-10*x);
  return fexact;
}*/
inline double exact(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);}

double step_value(int n){
  double h;
  h = 1/((double) n+1);
  cout << "h: " << h << endl;
  return h;
}

int main(int argc, char * argv[])
{
  //Task 1b)
  string outfilename;
  int m;

  if (argc <= 1){
    cout << "Type also: output file name and n; for max 10^n" << endl;
    exit(1);
  }
  else{
    outfilename = argv[1];
    m = atoi(argv[2]);
  }

  // forloop...? loop over power of 10
  for (int j = 1; j <= m; j++){
    int i; double h;
    int n = (int) pow(10.0,j);

    string outfile = outfilename;
    string argument = to_string(j);
    outfile.append(argument); // need the four loop to run
    outfile.append(".txt");

    cout << "n: " << n << ", outfile name: " << outfile << endl;

    // Define the arrays:
    vec a = ones<vec>(n); a.fill(-1);
    vec b = ones<vec>(n); b.fill(2);
    vec c = ones<vec>(n); c.fill(-1);

    vec x = linspace<vec>(0,1,n);
    vec f = zeros<vec>(n);
    vec btilde = zeros<vec>(n);
    vec v = zeros<vec>(n);    // Solution vector

    // Set up the array for f and Boundary conditions:
    x(0) = 0.0; x(n-1) = 1.0;
    btilde(0) = b(0);
    btilde(n-1) = b(n-1);
    //cout << "lol" << endl;
    for(int i=0; i<n; i++){
        f(i) = 100*exp(-10*x(i));
        }
    h = step_value(n);
    double hh = h*h;

    vec ftilde = hh*f;     // Set up ftilde vector

    // Forward subst:
    for (int i=1; i<n; i++){
      btilde(i) = b(i) - a(i-1)*a(i-1)/((double) btilde(i-1));
      ftilde(i) = ftilde(i) + ftilde(i-1)/((double) btilde(i-1));
    }

    // Backward subst:
    v(n-1) = ftilde(n-1)/btilde(n-1);
    for (int i=n-1; i>1; i--){
      v(i-1) = (ftilde(i-1) + v(i))/((double) btilde(i-1));
    }

    ofile.open(outfile);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << "     x:              approx:       exact:        relative error:        step length:" << endl;
    for (int i=1; i<n; i++){
      double eps_rel = fabs((v(i)-exact(x(i)))/exact(x(i)));
      ofile << setw(15) << setprecision(8) << x(i);
      ofile << setw(15) << setprecision(8) << v(i);
      ofile << setw(15) << setprecision(8) << exact(x(i));
      ofile << setw(15) << setprecision(8) << log10(eps_rel);
      ofile << setw(15) << setprecision(8) << h << endl;
    }
    ofile.close();
  }
  //btilde.print("btilde:");
  //ftilde.print("ftilde:");
  //v.print("v:");


  // For tridiagnoal
  /*int k;    // length of vector
  int n = atoi(argv[1]);    // Input parameter for number of rows and columns
  int a[n], b[n], c[n];     // Differential vectors
  // Define the elements in the vectors
  for (k=0; k<n+1; k++){
    a[k] = 1;
    b[k] = -2;
    c[k] = 1;
  }

  // Define the matrix
  mat A(n,n); A.zeros();
  int i, j; // Index parameters
  // Loop over the rows and columns
  for (i=0; i<=n-1; i++){
    for (j=0; j<=n-1; j++){
      // Set diagonal elements
      if (i<n-1){
        if (i==j){
          A(i,j) = b[i];
          A(i+1,j) = a[i];
          A(i,j+1) = c[i];
        }
      }
      // Set the last element in position i,j = n,n
      else if (i==n-1){
        if (i==j){A(i,j) = b[i];}
      }
    }
  }

  cout << size(A) << endl;
  A.print();
  */


  //return 0;
}
