//class for ForwardEuler
#ifndef FORWARDEULER_H
#define FORWARDEULER_H

//#include "system.h"
//#include "Method_Earth_sun.hpp"
#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <math.h>

//using namespace std;
//using namespace arma;

class ForwardEuler
{

//private:
public:
    arma::mat vel, pos, acc;
//    double x0, y0, z0;
//    double vx0, vy0, vz0;

//public:
    double m_dt;
    void Euler(double dt);
    arma::mat Integrate(int dim, int N, std::string filename, double eps, double dt);
//    double step (double stepsPrYear){return 1.0/stepsPrYear;}
    arma::mat InitialPosition(arma::mat pos);
    arma::mat InitialVelocity(arma::mat vel);
    arma::mat getPos(arma::mat pos){return pos;}
    arma::mat getVel(arma::mat vel){return vel;}
    arma::mat getAcc(arma::mat acc){return acc;}
};
#endif // FORWARDEULER_H
