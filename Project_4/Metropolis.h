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
#include "System.h"

using namespace arma;

class Metropolis
{
public:
    //System *Ising;

    void metropolis(int n_spin, int MCs, double Temp, vec ExpValues, int nTemp, std::string filename);//, System &Ising);
};
#endif // METROPOLIS_H
