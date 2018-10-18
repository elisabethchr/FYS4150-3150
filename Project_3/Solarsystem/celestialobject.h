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
    //: myname(name), m(mass), x(X), y(Y), z(Z), vx(vX), vy(vY), vz(vZ){} //X, Y, Z, vX, vY, vZ = initial conditions


/*
    //declare functions:
    void getPosition();
    void getVelocity();
    double getMass();
*/
/*
    double getX(double x){return x;}
    double getY(double y){return y;}
    double getZ(double z){return z;}
    double getvX(double vx){return vx;}
    double getvY(double vy){return vy;}
    double getvZ(double vz){return vz;}
*/
/*
    friend std::ostream& operator<<(std::ostream& out, const CelestialObject& obj) {
        return out << obj.myname << '\t' << obj.x << '\t' << obj.y
                   << '\t' << obj.vx << '\t' << obj.vy << std::endl;
    }
*/
};
#endif // CELESTIALOBJECT_H
