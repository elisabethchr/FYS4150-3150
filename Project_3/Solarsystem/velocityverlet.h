//class for VelocityVerlet
#ifndef VELOCITYVERLET_H
#define VELOCITYVERLET_H


#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "vec3.h"
#include "solarsystem.h"


class VelocityVerlet
{

//private:
public:
    SolarSystem *solar;
    arma::mat vel, pos;
    vec3 acc, acc_next;
    vec3 vel1, pos1;
    double mass;
    double m_dt;
    double h;

    void Verlet(double dt);

    void Integrate(SolarSystem &solar, double h, int n_years);

};

#endif // VELOCITYVERLET_H
