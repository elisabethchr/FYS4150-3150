//definitions for functions within VelocityVerlet
#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <math.h>
//#include "WriteToFile.hpp"
//#include "EnergyTest.hpp"
#include "velocityverlet.h"
#include "initialize.h"
#include "gravitationalforce.h"

using namespace std;
using namespace arma;

void VelocityVerlet::Verlet(double dt)// :m_dt(dt)
{
    m_dt = dt;
}

void VelocityVerlet::Integrate(int N, int dim, string filename, double eps, double dt){
  /*
  Compute the position of the planets using velocity verlet method.
  */
  double r2;
  double h = 1/((double) 365.25)*dt;//1.0/(365.25);
  double hh_half = h*h/2.0;
  double h_half = h/2.0;
  mat acc = zeros(dim, N); mat vel = zeros(dim, N); mat pos = zeros(dim, N);
  vel = Initialize::Initialvelocity(vel, dim, N);
  pos = Initialize::Initialposition(pos, dim, N);

  vec t = zeros(N);

  double r0 = sqrt(pow(pos(0,0),2) + pow(pos(1,0),2) + pow(pos(2,0),2));
  double v20 = vel(0,0)*vel(0,0) + vel(1,0)*vel(1,0) + vel(1,0)*vel(1,0);

  clock_t start, finish;
  start = clock();  // start timing

  for (int i=0; i<N-1; i++){
    double r = sqrt(pos(0, i)*pos(0, i) + pos(1, i)*pos(1, i) + pos(2, i)*pos(2, i));
    acc(0,i) = acc(1,i) = acc(2,i) = 0.0;
    for (int j=0; j<dim; j++){
      acc(j,i) = GravitationalForce::Force(pos(j,i), r);
      pos(j, i+1) = pos(j, i) + h*vel(j, i) + hh_half*acc(j,i);


    }
    double rnew = sqrt(pow(pos(0,i+1),2) + pow(pos(1,i+1),2) + pow(pos(2,i+1),2));
    for (int j=0; j<dim; j++){
      acc(j,i+1) = GravitationalForce::Force(pos(j,i+1), rnew);
      vel(j, i+1) = vel(j,i) + h_half*(acc(j,i+1) + acc(j,i));

    }
    //cout << pos(0,i+1) << " " << pos(1,i+1) << " " << pos(2, i+1) << endl;
  }
  finish = clock();   // end timing
  double time_used = (double)(finish - start)/(CLOCKS_PER_SEC );
  cout << setprecision(10) << "Time used: " << time_used << " s at " << N/365 <<" yr" << endl;

  cout << t(N-1) << " " << r2<< endl;

  double rN = sqrt(pow(pos(0, N-1),2) + pow(pos(1, N-1),2) + pow(pos(2, N-1),2));
  double v2N = vel(0,N-1)*vel(0,N-1) + vel(1,N-1)*vel(1,N-1) + vel(1,N-1)*vel(1,N-1);
  cout << "--> UnitTesting -->" << endl;
  cout << "Verlet works!" << endl;
}
