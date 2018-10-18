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
//#include "initialize.h"
#include "gravitationalforce.h"

void ForwardEuler::Euler(double dt)// :m_dt(dt)
{
    m_dt = dt;
}

arma::mat ForwardEuler::Integrate(vec3 position, vec3 velocity, double mass){ // int N, int dim, std::string obj){
  /*
  Compute the position of the planet using forward Euler method.
  */
//  double dt = 0.001;
  int N = 367;
  int dim = 3;
  double h = 1.0/((double) 365.25); //1.0/((double) 365.25);
  arma::mat vel = arma::zeros(dim, N);
  arma::mat pos = arma::zeros(dim, N);
  arma::mat acc = arma::zeros(dim, N);

//  Initialize initialize;

  vel(0, 0) = velocity(0); vel(1, 0) = velocity(1); vel(2, 0) = velocity(2);
  pos(0, 0) = position(0); pos(1, 0) = position(1); pos(2, 0) = position(2);

//  vel = InitialVelocity(vel);
//  pos = InitialPosition(pos);

  arma::vec t = arma::zeros(N);

  std::clock_t start, finish;
  start = std::clock();  // start timing

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
  finish = std::clock();   // end timing
  double time_used = (double)(finish - start)/(CLOCKS_PER_SEC );
  std::cout << std::setprecision(10) << "Time used: " << time_used << " s at " << N/365 <<" yr" << std::endl;
  std::cout << "Euler works!" << std::endl;

//  getPos(pos);
//  getVel(vel);
  return pos;
}
/*
  mat ForwardEuler::getPos(mat pos){
    return pos;
  }

  mat ForwardEuler::getVel(mat vel){
    return vel;
  }
*/
/*
  arma::mat ForwardEuler::InitialPosition(arma::mat pos){
    double x0 = pos(0,0) = 9.528047055398201E-01;       // AU
    double y0 = pos(1,0) = 3.053612869840809E-01;       // AU
    double z0 = pos(2,0) = -9.272902073041313E-05;      // AU
    std::cout <<"Initial position in x, y, z direction:" << std::endl;
    std::cout << x0 << " " << y0 <<" " << z0 << std::endl;
    return pos;
  }

  arma::mat ForwardEuler::InitialVelocity(arma::mat vel){
  // JPL values are AU/day so multiply with 365.25
    double vx0 = vel(0,0) = -5.428888690270241E-03*365.25;      // AU/yr
    double vy0 = vel(1,0) = 1.636353485946535E-02*365.25;       // AU/yr
    double vz0 = vel(2,0) = -4.491683144318728E-07*365.25;      // AU/yr
    std::cout <<"Initial velocity in x, y, z direction:" << std::endl;
    std::cout << vx0 << " " << vy0 <<" " << vz0 << std::endl;
    return vel;
  }
*/
