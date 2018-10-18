//class for VelocityVerlet
#ifndef VELOCITYVERLET_H
#define VELOCITYVERLET_H


#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <math.h>

//using namespace std;
//using namespace arma;

class VelocityVerlet
{

//private:
public:
    arma::mat vel, pos, acc;
//    double x0, y0, z0;
//    double vx0, vy0, vz0;

//public:
    double m_dt;
    void Verlet(double dt);
    arma::mat Integrate(int dim, int N, std::string filename, double eps, double dt);
//    double step (double stepsPrYear){return 1.0/stepsPrYear;}
    arma::mat InitialPosition(arma::mat pos);
    arma::mat InitialVelocity(arma::mat vel);
    arma::mat getPos(){return pos;}
    arma::mat getVel(){return vel;}
    arma::mat getAcc(){return acc;}
};

#endif // VELOCITYVERLET_H
