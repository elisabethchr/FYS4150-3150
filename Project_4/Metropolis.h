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
  void metropolisMPI(int n_spin, int myLoopBegin, int myLoopEnd, int MCs, double Temp,
                      vec ExpValues, std::string filename, int choise, int myRank);
  //void metropolisMPI(int n_spin, int MCs, double Temp, int nIntervalls, std::string filename,
  //                      int choise, int myRank, vec &ExpValues);
  vec ExpValues;
};
#endif // METROPOLIS_H
