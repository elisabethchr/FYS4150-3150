#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
#include <string>

#include "system.h"
#include "metropolis.h"

using namespace  std;
using namespace arma;

// output file

void System::writefile(int n_spin, int MCs, double Temp, vec expValues, string filename)
{
    double norm = 1.0/((double) MCs);
    double averageE = expValues(0)*norm;
    double averageE2 = expValues(1)*norm;
    double averageM = expValues(2)*norm;
    double averageM2 = expValues(3)*norm;
    double averageMabs = expValues(4)*norm;

    //double varE = (averageE2 - averageE*averageE)/n_spin/n_spin;
    //double varM = (averageM2 - averageMabs*averageMabs)/n_spin/n_spin;
    double Cv = (averageE2 - averageE*averageE)/n_spin/n_spin;
    double chi = (averageM2 - averageM*averageM)/n_spin/n_spin;

    ofstream ofile;
    string output = filename;
    string arg = to_string(n_spin);
    output.append(arg);
    output.append(".txt");
    ofile.open(output, ios::app);

    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << setprecision(8) << Temp;
    ofile << setw(15) << setprecision(8) << averageE/n_spin/n_spin;
    ofile << setw(15) << setprecision(8) << Cv/Temp/Temp;
    ofile << setw(15) << setprecision(8) << averageM/n_spin/n_spin;
    ofile << setw(15) << setprecision(8) << chi/Temp;
    ofile << setw(15) << setprecision(8) << averageMabs/n_spin/n_spin << "\n";

   // ofile.close();
}


void System::initialize(int n_spin, double Temp)
{
    /*
 * Initialize lattice matrices containing spins
 */
    spin_matrix.zeros(n_spin,n_spin);
    for (int y=0; y<n_spin; y++){
        for (int x=0; x<n_spin; x++){
            //      if (Temp < 1.5){
            spin_matrix(y,x) = 1;
            //      }
            M += (double) spin_matrix(y,x);
        }
    }
    for (int y=0; y<n_spin; y++){
        for (int x=0; x<n_spin; x++){
            E -= (double) (spin_matrix(y,x)*(spin_matrix(periodic(y,n_spin,-1), x) + spin_matrix(y,periodic(x,n_spin,-1))));
        }
    }
}

int System::periodic(int i, int n_spin, int add)
{
    return (i+n_spin+add) % (n_spin);       //integer remainder
}

arma::mat System::Lattice(){
return spin_matrix;
}

double System::Energy(){
return E;
}

double System::MagneticMoment(){
return M;
}
