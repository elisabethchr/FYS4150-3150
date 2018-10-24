//definitions for functions within CelestialObject

#include <iostream>
#include <armadillo>
#include "celestialobject.h"
#include "vec3.h"

//using namespace arma;

CelestialObject::CelestialObject(vec3 pos, vec3 vel, double mass_) {
    position = pos;
    velocity = vel;
    mass = mass_;
}

void CelestialObject::addForce(vec3 addF)
{
    // Adding force to an object.

    force += addF;
}

void CelestialObject::resetForce()
{
    force.zeros();
}
