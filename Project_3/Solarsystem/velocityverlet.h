//class for VelocityVerlet
#ifndef VELOCITYVERLET_H
#define VELOCITYVERLET_H


#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <math.h>

using namespace std;
using namespace arma;

class VelocityVerlet
{

//private:
public:
    mat vel, pos, acc;
//    double x0, y0, z0;
//    double vx0, vy0, vz0;

//public:
    double m_dt;
    void Verlet(double dt);
    void Integrate(int dim, int N, string filename, double eps, double dt);
//    double step (double stepsPrYear){return 1.0/stepsPrYear;}
    mat InitialPosition(mat pos);
    mat InitialVelocity(mat vel);
    mat getPos(){return pos;}
    mat getVel(){return vel;}
    mat getAcc(){return acc;}
};

#endif // VELOCITYVERLET_H
