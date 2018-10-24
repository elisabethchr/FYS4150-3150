#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "include_headers.h"

//to run in terminal:
//-> c++ -o run.x MainSS.cpp velocityverlet.cpp forwardeuler.cpp solarsystem.cpp celestialobject.cpp writetofile.cpp readfile_test.cpp vec3.cpp -larmadillo

int main(int nArgs, char **arguments)
{
    int numTimesteps = 366*2;    // 2 years minimum
    if (nArgs >= 2) numTimesteps = atoi(arguments[1])*366; // n years

    //  double epsilon = 1e-5;
    double h = 1/(365.25);     //declaring step value --> if you want the step value to reflect each day, then <double h = 1.0/((double) 365.25)>;
    int n_years = 20;               //declaring number of years

    //Insert the textfiles including the initial positions, velocities and masses for the different objects,
    //i.e. Solarsystem sol(<position_file>, <velocity_file>, <masses_file>)
    SolarSystem sol("Initialposition.txt", "Initialvelocity.txt", "masses.txt");

    ForwardEuler SolvePositions_euler;
    SolvePositions_euler.Integrate(sol, h, n_years);

    VelocityVerlet SolvePositions_verlet;
    SolvePositions_verlet.Integrate(sol, h, n_years);
    return 0;

//Energy is conserved for planet0, dE = 3.44198e-006 --> for all planets after 248yr
}
