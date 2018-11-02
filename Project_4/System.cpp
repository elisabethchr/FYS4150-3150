#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
#include <string>

#include "System.h"
#include "Metropolis.h"

using namespace  std;
using namespace arma;

// output file
ofstream ofile;

int System::periodic(int i, int n_spin, int add)
{
  return (i+n_spin+add) % (n_spin);
}

void System::writefile(int n_spin, int MCs, double Temp, vec average)
{
  double norm = 1.0/((double) MCs);
  double averageE = average(0)*norm;
  double averageE2 = average(1)*norm;
  double averageM = average(2)*norm;
  double averageM2 = average(3)*norm;
  double averageMabs = average(4)*norm;

  //double varE = (averageE2 - averageE*averageE)/n_spin/n_spin;
  //double varM = (averageM2 - averageMabs*averageMabs)/n_spin/n_spin;
  double Cv = (averageE2 - averageE*averageE)/n_spin/n_spin;
  double chi = (averageM2 - averageM*averageM)/n_spin/n_spin;

  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setw(15) << setprecision(8) << Temp;
  ofile << setw(15) << setprecision(8) << averageE/n_spin/n_spin;
  ofile << setw(15) << setprecision(8) << Cv/Temp/Temp;
  ofile << setw(15) << setprecision(8) << averageM/n_spin/n_spin;
  ofile << setw(15) << setprecision(8) << chi/Temp;
  ofile << setw(15) << setprecision(8) << averageMabs/n_spin/n_spin << endl;
}

void System::initialize(int n_spin, double Temp, mat spin_matrix, double &E, double &M)
{
  for (int y=0; y<n_spin; y++){
    for (int x=0; x<n_spin; x++){
      if (Temp < 1.5){
        spin_matrix(y,x) = 1;
      }
      M += (double) spin_matrix(y,x);
    }
  }
  for (int y=0; y<n_spin; y++){
    for (int x=0; x<n_spin; x++){
      E -= (double) (spin_matrix(y,x)*(spin_matrix(periodic(y,n_spin,-1), x) + spin_matrix(y,periodic(x,n_spin,-1)));
    }
  }
}
