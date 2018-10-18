//definitions for functions withn class GravitationalForce
#include "gravitationalforce.h"
#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <math.h>
//#include "vec3.h"

//using namespace std;
// namespace arma;

double GravitationalForce::Force(double pos, double r)
{
/*
  Compute the Gravitational force, and retun F/M_planet
*/
  double G = 6.67e-11;      // N m^2 kg^2 Gravitational constant
  double Msun = 1.99e+30;   // kg
  double pi = acos(-1);
  double GtimesMsun = 4*pi*pi; //AU^3/yr^2
  double F = -GtimesMsun * pos/((double)r*(double)r*(double)r);

  return F;

}

//alternative way:
    double GravitationalForce::Gravity(CelestialObject *a, CelestialObject *b)
    {
        object_a = a;
        object_b = b;
        double pi = acos(-1);
        double G = 4*pi*pi;

        vec3 F1;
        vec3 F2;
        vec3 r_a = object_a -> position;
        vec3 r_b = object_b -> position;

        double distance = (r_b - r_a).length();

//Calculating the gravitational force
    F1 = G*object_a->mass * object_b->mass*(r_b-r_a)/((double) pow(distance, 3));
    F2 = (-1)*F1;
    object_a -> addForce(F1);
    object_b -> addForce(F2);
    }
