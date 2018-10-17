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

using namespace std;
using namespace arma;

class ForwardEuler
{

//private:
public:
    mat vel, pos, acc;
//    double x0, y0, z0;
//    double vx0, vy0, vz0;

//public:
    double m_dt;
    void Euler(double dt);
    mat Integrate(int dim, int N, string filename, double eps, double dt);
//    double step (double stepsPrYear){return 1.0/stepsPrYear;}
    mat getPos(mat pos){return pos;}
    mat getVel(mat vel){return vel;}
    mat getAcc(mat acc){return acc;}
};
#endif // FORWARDEULER_H
