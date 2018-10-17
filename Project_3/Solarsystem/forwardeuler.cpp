//definitions for functions within class ForwardEuler
#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <math.h>
//#include "WriteToFile.hpp"
//#include "EnergyTest.hpp"
#include "forwardeuler.h"
#include "initialize.h"
#include "gravitationalforce.h"

using namespace std;
using namespace arma;

void ForwardEuler::Euler(double dt)// :m_dt(dt)
{
    m_dt = dt;
}

mat ForwardEuler::Integrate(int N, int dim, string filename, double eps, double dt){
  /*
  Compute the position of the planet using forward Euler method.
  */
  double h = 1.0/((double) 365.25)*dt; //1.0/((double) 365.25);
  mat vel = zeros(dim, N); mat pos = zeros(dim, N); mat acc = zeros(dim, N);

  Initialize initialize;
  vel = initialize.Initialvelocity(vel, dim, N);
  pos = initialize.Initialposition(pos, dim, N);

  vec t = zeros(N);
  double r0 = sqrt(pow(pos(0,0),2) + pow(pos(1,0),2) + pow(pos(2,0),2));
  double v20 = vel(0,0)*vel(0,0) + vel(1,0)*vel(1,0) + vel(1,0)*vel(1,0);

  clock_t start, finish;
  start = clock();  // start timing

  // Integrating loop:
  for (int i=0; i<N-1; i++){
    double r = sqrt(pow(pos(0, i),2) + pow(pos(1, i),2) + pow(pos(2, i),2));
    double v2 = vel(0,i)*vel(0,i) + vel(1,i)*vel(1,i) + vel(1,i)*vel(1,i);
    t(i+1) = t(i) + h;
    for (int j=0; j<dim; j++){

      acc(j,i) = GravitationalForce::Force(pos(j,i), r);
      vel(j, i+1) = vel(j, i) + h*acc(j,i);
      pos(j, i+1) = pos(j, i) + h*vel(j,i);

    }

  }
  finish = clock();   // end timing
  double time_used = (double)(finish - start)/(CLOCKS_PER_SEC );
  cout << setprecision(10) << "Time used: " << time_used << " s at " << N/365 <<" yr" << endl;
  cout << "Euler works!" << endl;

//  getPos(pos);
//  getVel(vel);
  return pos;
}
  mat getPos(mat pos){
    return pos;
  }

  mat getVel(mat vel){
    return vel;
  }
