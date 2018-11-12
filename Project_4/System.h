#ifndef SYSTEM_H
#define SYSTEM_H

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
#include <string>

using namespace arma;

class System
{
public:
  int periodic(int i, int n_spin, int add);
  void writefile(int n_spin, int mcs, double Temp, vec average, std::string filename, int cycles, int N_accepted);
  void initialize(int n_spin, double Temp, double Ein, double Min, int choise);//, double &E, double &M);
  //mat spin_matrix;
  mat Lattice();
  double MagneticMoment();
  double Energy();

private:
  mat spin_matrix;
  double E;
  double M;
};
#endif
