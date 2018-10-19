#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
//#include "Method_Earth_sun.hpp"
#include "forwardeuler.h"
//#include "velocityverlet.h"
#include "celestialobject.h"
#include "solarsystem.h"
#include "vec3.h"
//#include "initialize.h"
#include "Readfile.h"
//#include "readfile.hpp"

//using namespace std;
//using namespace arma;


int main(int nArgs, char **arguments)
{
  int numTimesteps = 366*2; // 2 years minimum
  if (nArgs >= 2) numTimesteps = atoi(arguments[1])*366; // n years
  int dim = 3;

  double stepsPrYear = 365.25;
  double epsilon = 1e-5;
  double dt = 0.001;
  arma::mat pos_euler;
  arma::mat pos_verlet;


  SolarSystem sol("Initialposition.txt", "Initialvelocity.txt");
  ForwardEuler SolvePositions;
  SolvePositions.Integrate(sol);

//  SolarSystem solarsystem;
// CelestialObject &sun = solarsystem.createCelestialBody( vec3(0,0,0), vec3(0,0,0), 1.0 );
//  std::string obj = "sun";
//  ForwardEuler integrate_euler;     //need object of the class (integrate_euler) to use the member function Integrate within the class
  /*
  pos_euler = integrate_euler.Integrate(numTimesteps, dim, obj, epsilon, dt);
  std::cout  << "pos_euler: " << std::endl;
  std::cout << pos_euler << std::endl;

  VelocityVerlet integrate_verlet;  //need object of the class (integrate_verlet) to use the member function Integrate within the class
  pos_verlet = integrate_verlet.Integrate(numTimesteps, dim, obj, epsilon, dt);
*/

  // Call initial values: 0-CoM, 1-Mercery,2-venus,3-Earth,4-mars,5-jupiter,6-saturn,
  // 7-uranus,8-neptune,9-pluto. Both position and velocity is listed: x, y, z direction

/* -----------> this is put in the solarsystem.cpp file
  arma::mat pos0 = Readfile("Initialposition.txt");     // AU
  arma::mat vel0 = Readfile("Initialvelocity.txt");     // AU/yr
  //pos0.print("initpos");
  //vel0.print("initvel");
  arma::vec pos0_CoM = arma::zeros(3); arma::vec vel0_CoM = arma::zeros(3);           // Center of mass
  arma::vec pos0_Mercery = arma::zeros(3); arma::vec vel0_Mercery = arma::zeros(3);
  arma::vec pos0_Venus = arma::zeros(3); arma::vec vel0_Venus = arma::zeros(3);
  arma::vec pos0_Earth = arma::zeros(3); arma::vec vel0_Earth = arma::zeros(3);
  arma::vec pos0_Mars = arma::zeros(3); arma::vec vel0_Mars = arma::zeros(3);
  arma::vec pos0_Jupiter = arma::zeros(3); arma::vec vel0_Jupiter = arma::zeros(3);
  arma::vec pos0_Saturn = arma::zeros(3); arma::vec vel0_Saturn = arma::zeros(3);
  arma::vec pos0_Uranus = arma::zeros(3); arma::vec vel0_Uranus = arma::zeros(3);
  arma::vec pos0_Neptune = arma::zeros(3); arma::vec vel0_Neptune = arma::zeros(3);
  arma::vec pos0_Pluto = arma::zeros(3); arma::vec vel0_Pluto = arma::zeros(3);

  for (int j=0; j<3;j++){
    pos0_CoM(j) = pos0(0,j); vel0_CoM(j) = vel0(0,j);
    pos0_Mercery(j) = pos0(1,j); vel0_Mercery(j)=vel0(1,j);
    pos0_Venus(j) = pos0(2,j); vel0_Venus(j) = vel0(2,j);
    pos0_Earth(j) = pos0(3,j); vel0_Earth(j) = vel0(3,j);
    pos0_Mars(j) = pos0(4,j); vel0_Mars(j) = vel0(4,j);
    pos0_Jupiter(j) = pos0(5,j); vel0_Jupiter(j) = vel0(5,j);
    pos0_Saturn(j) = pos0(6,j); vel0_Saturn(j) = vel0(6,j);
    pos0_Uranus(j) = pos0(7,j); vel0_Uranus(j) = vel0(7,j);
    pos0_Neptune(j) = pos0(8,j); vel0_Neptune(j) = vel0(8,j);
    pos0_Pluto(j) = pos0(9,j); vel0_Pluto(j) = vel0(9,j);
    }
<------------- */


//    ForwardEuler integrate_euler;
/*
    CelestialObject &sun = SolarSystem.createCelestialObject( vec3(0,0,0), vec3(0,0,0), 1.0 );

    // We don't need to store the reference, but just call the function without a left hand side
    double x0 = 9.528047055398201E-01;       // AU
    double y0 = 3.053612869840809E-01;       // AU
    double z0 = -9.272902073041313E-05;      // AU
    double vx0 = -5.428888690270241E-03*365.25;      // AU/yr
    double vy0 = 1.636353485946535E-02*365.25;       // AU/yr
    double vz0 = -4.491683144318728E-07*365.25;      // AU/yr

    CelestialObject &earth = SolarSystem.createCelestialObject(vec3 (pos0_Earth(0), pos0_Earth(1), pos0_Earth(2)), vec3 (vel0_Earth(0), vel0_Earth(1), vel0_Earth(2)), 3e-6 );

    // To get a list (a reference, not copy) of all the bodies in the solar system, we use the .objects() function
*/
//    ForwardEuler integrate_euler;
//   integrate_euler.Integrate(&sun, &earth);

  //SolarSystem = SolarSystem;
  //ForwardEuler Euler;
  //double &a = Euler.Force(1,1);
  //cout << a << endl;
  //Euler(dim, numTimesteps, "test.txt");

  return 0;
}
