#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
//#include "Method_Earth_sun.hpp"
#include "forwardeuler.h"
//#include "initialize.h"
#include "Readfile.h"
//#include "readfile.hpp"
using namespace std;
using namespace arma;


int main(int nArgs, char **arguments)
{
  int numTimesteps = 366*2; // 2 years minimum
  if (nArgs >= 2) numTimesteps = atoi(arguments[1])*366; // n years
  int dim = 3;

  // Call initial values: 0-CoM, 1-Mercery,2-venus,3-Earth,4-mars,5-jupiter,6-saturn,
  // 7-uranus,8-neptune,9-pluto. Both position and velocity is listed: x, y, z direction
  mat pos0 = Readfile("Initialposition.txt");     // AU
  mat vel0 = Readfile("Initialvelocity.txt");     // AU/yr
  pos0.print("initpos");
  vel0.print("initvel");
  vec pos0_CoM = zeros(3); vec vel0_CoM = zeros(3);           // Center of mass
  vec pos0_Mercery = zeros(3); vec vel0_Mercery = zeros(3);
  vec pos0_Venus = zeros(3); vec vel0_Venus = zeros(3);
  vec pos0_Earth = zeros(3); vec vel0_Earth = zeros(3);
  vec pos0_Mars = zeros(3); vec vel0_Mars = zeros(3);
  vec pos0_Jupiter = zeros(3); vec vel0_Jupiter = zeros(3);
  vec pos0_Saturn = zeros(3); vec vel0_Saturn = zeros(3);
  vec pos0_Uranus = zeros(3); vec vel0_Uranus = zeros(3);
  vec pos0_Neptune = zeros(3); vec vel0_Neptune = zeros(3);
  vec pos0_Pluto = zeros(3); vec vel0_Pluto = zeros(3);

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


  //SolarSystem = SolarSystem;
  //ForwardEuler Euler;
  //double &a = Euler.Force(1,1);
  //cout << a << endl;
  //Euler(dim, numTimesteps, "test.txt");

  return 0;
}
