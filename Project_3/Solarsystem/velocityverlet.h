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

//using namespace std;
//using namespace arma;

class VelocityVerlet
{

//private:
public:
    SolarSystem *solar;
    arma::mat vel, pos;
    vec3 acc, acc_next;
    vec3 vel1, pos1;
    double mass;
//    double x0, y0, z0;
//    double vx0, vy0, vz0;

//public:
    double m_dt;
    void Verlet(double dt);

    arma::mat Integrate(SolarSystem &solar);
//    double step (double stepsPrYear){return 1.0/stepsPrYear;}
    //arma::mat InitialPosition(arma::mat pos);
    //arma::mat InitialVelocity(arma::mat vel);
    //arma::mat getPos(){return pos;}
    //arma::mat getVel(){return vel;}
    //arma::mat getAcc(){return acc;}
};

#endif // VELOCITYVERLET_H
