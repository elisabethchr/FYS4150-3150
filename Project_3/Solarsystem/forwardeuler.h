//class for ForwardEuler
#ifndef FORWARDEULER_H
#define FORWARDEULER_H

#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "vec3.h"
#include "solarsystem.h"


class ForwardEuler
{

//private:
public:
    SolarSystem *solar;
    arma::mat vel, pos;
    vec3 acc;
    vec3 vel1, pos1;
    double mass;


    //public:
    double m_dt;
    void Euler(double dt);

    arma::mat Integrate(SolarSystem &solar, double h, int n_years);
    
};
#endif // FORWARDEULER_H
