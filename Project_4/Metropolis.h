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
  void Metropolis(int n_spin, int MCs, double Temp, vec ExpValues);
};
#endif // METROPOLIS_H
