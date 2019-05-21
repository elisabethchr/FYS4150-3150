#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include "time.h"
#include <string>
#include <math.h>

using namespace std;
using namespace arma;
ofstream ofile;

// f(x) function
inline double f_array(double x){
    return 100*exp(-10*x);
}

// Exact solution of u''(x) = f(x)
inline double exact(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);}

// Calcultate the step lenght
double step_value(int n){
  double h;
  h = 1/((double) n+1);
  cout << "h: " << h << endl;
  return h;
}

int main(int argc, char * argv[])
{
  //Task 1c) With LU decomposition
  string outfilename;
  int m; int aa; int bb; int cc;
  // Check for input arguments
  if (argc <= 1){
    cout << "Type also: output file name, n; for max 10^n, and the diagnal elements a, b, c" << endl;
    exit(1);
  }
  else{
    outfilename = argv[1];
    m = atoi(argv[2]);
    aa = atoi(argv[3]);
    bb = atoi(argv[4]);
    cc = atoi(argv[5]);
  }

  //Loop over the power of 10: for each 10^n
  for (int j = 1; j <= m; j++){
    int i; double h;
    int n = (int) pow(10.0,j);

    // Define the output filem name
    string outfile = outfilename;
    string argument = to_string(j);
    outfile.append(argument);
    outfile.append(".txt");

    cout << "n: " << n << ", outfile name: " << outfile << endl;

    // Define the arrays:
    vec a = ones<vec>(n); a.fill(aa);
    vec b = ones<vec>(n); b.fill(bb);
    vec c = ones<vec>(n); c.fill(cc);

    vec x = linspace<vec>(0,1,n);
    vec f = zeros<vec>(n);
    vec v = zeros<vec>(n);    // Solution vector

    // Set up the array for f and Boundary conditions:
    x(0) = 0.0; x(n-1) = 1.0;

    // Calculate f(x)
    for(int i=0; i<n; i++){
        f(i) = 100*exp(-10*x(i));
        }

    h = step_value(n);      // Call for the step length   !!!
    double hh = h*h;        //      !!
    vec ftilde = hh*f;      // Set up ftilde;     !!

    // Timing of the aglorithm start
    clock_t start, finish;
    start = clock();

    // Define the tridiagnal matrix
    mat A(n,n); A.zeros(); int k;
    for(i=0; i<=n-1; i++){
        for(k=0; k<=n-1; k++){
          if(i<n-1){
            if(i==k){
            A(i+1,k) = a(i);
            A(i, k) = b(i);
            A(i, k+1) = c(i);
            }
          }
          else if(i==n-1){
            if(i==k){A(i, k) = b(i);}
          }
        }
    }

    //A.print("A=");
    mat Up; mat Low;
    lu(Low,Up,A);     // Lower upper decomposition
    vec Solution = solve(A,ftilde);    // Solve Ax=b
    //Solution.print("Solution to Ax=b:");

    finish = clock();     // Stop timing
    double time_used = (double) (finish - start)/(CLOCKS_PER_SEC );
    // Print time used
    cout << setprecision(10) << "Time used: " << time_used << " s at n=" << n << endl;

    // Write to file, with values of x, approximation, exact solution, relative error and step length
    ofile.open(outfile);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << "     x:              approx:       exact:        relative error:        step length:" << endl;
    for (int i=1; i<n; i++){
      double eps_rel = fabs((Solution(i)-exact(x(i)))/exact(x(i)));
      ofile << setw(15) << setprecision(8) << x(i);
      ofile << setw(15) << setprecision(8) << Solution(i);
      ofile << setw(15) << setprecision(8) << exact(x(i));
      ofile << setw(15) << setprecision(8) << log10(eps_rel);
      ofile << setw(15) << setprecision(8) << h << endl;
    }
    ofile.close();

  }
//return 0;
}
