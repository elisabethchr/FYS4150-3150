//definitions for functions within CelestialObject

#include <iostream>
#include <armadillo>
#include "celestialobject.h"

using namespace arma;

CelestialObject::CelestialObject(mat pos, mat vel, double mass_) {
    position = pos;
    velocity = vel;
    mass = mass_;
}
