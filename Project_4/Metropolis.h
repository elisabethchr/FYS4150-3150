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

class Metropolis
{
public:
  void metropolis(int n_spin, int MCs, double Temp, vec ExpValues, std::string filename, int choise);
};
#endif // METROPOLIS_H
