//class for CelestialObject
#ifndef CELESTIALOBJECT_H
#define CELESTIALOBJECT_H

#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "vec3.h"

//using namespace arma;

class CelestialObject
{
//private:
//    std::string myname;
//    double mass;
//    mat position;
//    mat velocity;


public:
    vec3 position;
    vec3 velocity;
    vec3 force;
    double mass;

    CelestialObject(vec3 pos, vec3 vel, double mass_);
    void addForce(vec3 addF);
    void resetForce();
    
};
#endif // CELESTIALOBJECT_H
