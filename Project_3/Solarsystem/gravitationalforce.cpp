//definitions for functions withn class GravitationalForce
#include "gravitationalforce.h"
#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <math.h>

using namespace std;
using namespace arma;

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
