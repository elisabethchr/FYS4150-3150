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
#include "celestialobject.h"
#include "solarsystem.h"

void ForwardEuler::Euler(double dt)// :m_dt(dt)
{
    m_dt = dt;
}

arma::mat ForwardEuler::Integrate(SolarSystem &input_system){//vec3 position, vec3 velocity, double mass){ // int N, int dim, std::string obj){
  /*
  Compute the position of the planet using forward Euler method.
  */
//  double dt = 0.001;
  std::cout << "euler?" << std::endl;
  solar = &input_system;
  int  NumberOfObjects = solar->numberOfObjects();

  int N = 367;
  int dim = 3;
  double h = 1.0/((double) 365.25); //1.0/((double) 365.25);
  arma::mat vel = arma::zeros(dim, N);
  arma::mat pos = arma::zeros(dim, N);


//  Obtain the initial velocities and positions for planets in m_objects;
  //vec3 posit1, veloc1;
  //vec3 posit2, veloc2;
    std::cout << "before extracting bodies" << std::endl;
    std::vector<CelestialObject> &bodies = solar->objects();
    arma::mat posit = arma::zeros(dim+1, bodies.size());
    arma::mat veloc = arma::zeros(dim, bodies.size());
    arma::vec mass = arma::zeros(bodies.size());
    for (int i=0; i<bodies.size(); i++){
      CelestialObject &obj = bodies[i];
        std::cout << "i = " << i << std::endl;
      for (int j=0; j < dim; j++){
          veloc(j, i) = obj.velocity(j);
        }
      for(int j=0; j < dim+1; j++){
          posit(j, i) = obj.position(j);
          }
    mass(i) = obj.mass;
    }
    std::cout << "mass" << mass << std::endl;


  vel(0, 0) = veloc(0); vel(1, 0) = veloc(1); vel(2, 0) = veloc(2);
  pos(0, 0) = posit(0); pos(1, 0) = posit(1); pos(2, 0) = posit(2);

//  vel = InitialVelocity(vel);
//  pos = InitialPosition(pos);

  arma::vec t = arma::zeros(N);

  std::clock_t start, finish;
  start = std::clock();  // start timing

  // Integrating loop:
//for(int k=0; k < objects.size(); k++){
//    CelestialObject &obj1 = objects[k];
//    for(int l=k+1; l < objects.size(); l++){
//        CelestialObject &obj2 = objects[l];
  for(int k=0; k<NumberOfObjects; k++){
  std::cout << "printing?" << std::endl;
    for(int i=0; i<N-1; i++){
        for (int j=0; j<dim; j++){
            vec3 acc = bodies[k].force/bodies[k].mass;
            vel(j, i+1) = vel(j, i) + h*acc(j);
            pos(j, i+1) = pos(j, i) + h*vel(j,i);
        }
    }
    std::cout << "pos for planet " << k << std::endl;
    pos.print();
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
