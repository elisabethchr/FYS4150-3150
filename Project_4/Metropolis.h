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
#include <mpi.h>

using namespace arma;

class Metropolis
{
public:
  void metropolis(int n_spin, int MCs, double Temp, vec ExpValues, std::string filename, int choise);
  void metropolisMPI(int n_spin, int MCs, double Temp, vec ExpValues, std::string filename, int choise);
};
#endif // METROPOLIS_H
