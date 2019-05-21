#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cctype>
#include <math.h>
#include <armadillo>

using namespace std;
using namespace arma;

ofstream ofile;         //create file for output

inline double f_array(double x){
       return 100*exp(-10.0*x);
}


double step_value(int n){
    double h;
    h = 1/((double) n);
    return h;
}

inline double exact(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);}
/*
vec exact(double x){
    vec u;
    u = 1 - (1-exp(-10.0))*x - exp(-10.0*x);
    return u;
}
*/

int main(int argc, char* argv[])
{
//task 1b)
    int power;
    string filename;        //allow for changeable filenames (for output)

    if(argc<=1){
    cout << "Not enough arguments, need output filename and exponential of 10^n" << endl;
    exit(1);
    }
    else{
    filename=argv[1];    //argument vector
    power = atoi(argv[2]);        //ascii to integer of argument vector
    }

//create ouput files for each 10^n power:
    for(int i=1; i<=power; i++){
    int n = pow(10, i);
//declare new filename
    string outputdata = filename;
//convert power 10^i to string
    string argument = to_string(i);
//allow each filename to be changed by adding argument i, e.g. outputdata1, outputdata2, etc.
    outputdata.append(argument);
    outputdata.append(".txt");

    double h;
    vec a = ones<vec>(n+1); a.fill(-1);
    vec d = ones<vec>(n+1); d.fill(2);        //i.e. b-vector as defined in project 1
    vec x = linspace<vec>(0, 1, n+1);
    vec f = zeros<vec>(n);

    for(int i=1; i<n; i++){
       f[i] = 100*exp(-10*x[i]);
    }

    h = step_value(n);
    cout << "Step value: " << h << endl;

    vec v = zeros<vec>(n+1);
    vec b_tilde = pow(h, 2)*f;
    vec d_tilde = zeros<vec>(n+1);
    d_tilde[0] = d[0];

    for(int i=1; i<n; i++){
        d_tilde[i] =(i+1.0)/((double) i);
    }

/*
 * Gaussian elimination:
 * Start first with a forward substitution (to create upper triangular matrix),
 * then a backward substitution (to create lower triangular matrix)
 */

//forward substitution:
    for (int i=2; i<n; i++){
//        d_tilde[i] = d[i] - a[i-1]/((double) d_tilde[i-1]);
        b_tilde[i] = b_tilde[i] + b_tilde[i-1]/(d_tilde[i-1]);
        }
    v[n-1] = b_tilde[n-1]/d_tilde[n-1];

//backward substitution:
    for(int i=n-2; i>0; i--){
        v[i] =(b_tilde[i] + v[i+1])/((double) d_tilde[i]);
    }

//open file and write out results
    ofile.open(outputdata);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << "          x:          approx:         exact:    relative error:     h:" << endl;

    for(int i=1; i<n; i++){
    double RelativeError = fabs((exact(x[i]) - v[i])/(exact(x[i])));
    ofile << setw(15) << setprecision(8) << x[i];
    ofile << setw(15) << setprecision(8) << v[i];
    ofile << setw(15) << setprecision(8) << exact(x[i]);
    ofile << setw(15) << setprecision(8) << log10(RelativeError);
    ofile << setw(15) << setprecision(8) << h << endl;
    }
    ofile.close();
    }
return 0;
}
