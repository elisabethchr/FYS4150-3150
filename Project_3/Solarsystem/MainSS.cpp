#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "include_headers.h"

/*
#include "forwardeuler.h"
//#include "velocityverlet.h"
#include "celestialobject.h"
#include "solarsystem.h"
#include "vec3.h"
#include "Readfile.h"
*/
//using namespace std;
//using namespace arma;


int main(int nArgs, char **arguments)
{
  int numTimesteps = 366*2; // 2 years minimum
  if (nArgs >= 2) numTimesteps = atoi(arguments[1])*366; // n years
  int dim = 3;

//  double stepsPrYear = 365.25;
//  double epsilon = 1e-5;
  double dt = 0.001;


  SolarSystem sol("Initialposition.txt", "Initialvelocity.txt", "masses.txt");
//  ForwardEuler SolvePositions_euler;
//  SolvePositions_euler.Integrate(sol);
  VelocityVerlet SolvePositions_verlet;
  SolvePositions_verlet.Integrate(sol);
  return 0;
}
