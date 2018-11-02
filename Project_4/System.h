#ifndef METROPOLIS_H
#define METROPOLIS_H

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
  void writefile(int n_spin, int mcs, double Temp, vec average);
  void initialize(int n_spin, double Temp, mat spin_matrix, double &E, double &M);

};
#endif
