//class for GravitationalForce --> this we can push into solarsystem when we have a general algortihm for the force
#ifndef GRAVITATIONALFORCE_H
#define GRAVITATIONALFORCE_H

#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "celestialobject.h"
#include "vec3.h"

//using namespace std;
//using namespace arma;

class GravitationalForce
{
    CelestialObject* object_a;
    CelestialObject* object_b;
    double G;

public:
    static double Force(double pos, double r);
    double Gravity(CelestialObject *a, CelestialObject *b);
};

#endif // GRAVITATIONALFORCE_H
