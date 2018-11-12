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

int System::periodic(int i, int n_spin, int add)
{
  return (i+n_spin+add) % (n_spin);
}

void System::writefile(int n_spin, int MCs, double Temp, vec average, string filename, int cycle, int counter)
{
  int Cycles = cycle;
  int N_accepted = counter;
  double norm = 1.0/((double) Cycles);

  double averageE = average(0)*norm;
  double averageE2 = average(1)*norm;
  double averageM = average(2)*norm;
  double averageM2 = average(3)*norm;
  double averageMabs = average(4)*norm;
  //cout << N_accepted << endl;

  double varE = (averageE2 - averageE*averageE)/(n_spin*n_spin);
  double varM = (averageM2 - averageMabs*averageMabs)/(n_spin*n_spin);
  double Cv = (averageE2 - averageE*averageE)/(Temp*Temp*n_spin*n_spin);
  double chi = (averageM2 - averageM*averageM)/(Temp*n_spin);

  ofstream ofile;
  string output = filename;
  string arg = to_string(MCs);
  output.append(arg);
  output.append(".txt");
  ofile.open(output, ios::app | ios::out);

  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setw(10) << setprecision(8) << cycle;
  ofile << setw(15) << setprecision(8) << Temp;
  ofile << setw(15) << setprecision(8) << averageE/(n_spin*n_spin);//;
  ofile << setw(15) << setprecision(8) << Cv;// /(n_spin*n_spin);
  ofile << setw(15) << setprecision(8) << averageM/(n_spin*n_spin);//;
  ofile << setw(15) << setprecision(8) << chi;// /(n_spin*n_spin);
  ofile << setw(15) << setprecision(8) << averageMabs/(n_spin*n_spin);
  ofile << setw(10) << setprecision(8) << N_accepted << "\n";// << "\n";
  ofile.close();
}

void System::initialize(int n_spin, double Temp, double Ein, double Min, int choise)//, double &E, double &M)
{

  E = Ein; M = Min;
  spin_matrix.zeros(n_spin, n_spin);
  random_device rd;
  mt19937_64 gen(rd());
  uniform_real_distribution<double> RNG(0.0,1.0);

  if (choise==0){
    for (int y=0; y<n_spin; y++){
      for (int x=0; x<n_spin; x++){
        spin_matrix(y,x) = 1.0;
        M += (double) spin_matrix(y,x);
      }
    }
  }
  else if (choise==1){
    for (int y=0; y<n_spin; y++){
      for (int x=0; x<n_spin; x++){
        if (RNG(gen) <= 0.5) spin_matrix(y,x) = -1.0;
        else if (RNG(gen) > 0.5) spin_matrix(x,y) = 1.0;

        M += (double) spin_matrix(y,x);
      }
    }
  }

  for (int y=0; y<n_spin; y++){
    for (int x=0; x<n_spin; x++){
      E -= (double) (spin_matrix(y,x)*(spin_matrix(periodic(y,n_spin,-1), x) + spin_matrix(y,periodic(x,n_spin,-1))));
    }
  }
}

mat System::Lattice()
{
  return spin_matrix;
}

double System::Energy()
{
  return E;
}

double System::MagneticMoment()
{
  return M;
}
