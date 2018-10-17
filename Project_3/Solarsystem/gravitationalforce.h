//class for GravitationalForce --> this we can push into solarsystem when we have a general algortihm for the force
#ifndef GRAVITATIONALFORCE_H
#define GRAVITATIONALFORCE_H

#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <math.h>

using namespace std;
using namespace arma;

class GravitationalForce
{
public:
    static double Force(double pos, double r);
};

#endif // GRAVITATIONALFORCE_H
