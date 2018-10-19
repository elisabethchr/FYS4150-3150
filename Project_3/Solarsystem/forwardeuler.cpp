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
#include "writetofile.h"

void ForwardEuler::Euler(double dt)// :m_dt(dt)
{
    m_dt = dt;
}

arma::mat ForwardEuler::Integrate(SolarSystem &input_system){//vec3 position, vec3 velocity, double mass){ // int N, int dim, std::string obj){
  /*
  Compute the position of the planet using forward Euler method.
  */
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
    std::vector<CelestialObject> &bodies = solar->objects();
    arma::mat posit = arma::zeros(dim, bodies.size());
    arma::mat veloc = arma::zeros(dim, bodies.size());
    arma::vec mass = arma::zeros(bodies.size());

    for (int i=0; i<bodies.size(); i++){
      CelestialObject &obj = bodies[i];
      for (int j=0; j < dim; j++){
          veloc(j, i) = obj.velocity(j);
        }
      for(int j=0; j < dim; j++){
          posit(j, i) = obj.position(j);
          }
    mass(i) = obj.mass;
    }
//std::cout << "veloc =" << veloc << std::endl;

  arma::vec t = arma::zeros(N);

  std::clock_t start, finish;
  start = std::clock();  // start timing



//forward euler

  for(int k=0; k<NumberOfObjects; k++){
//initializing to position and velociy matrices used in forward euler
    vel(0, 0) = veloc(0, k); vel(1, 0) = veloc(1, k); vel(2, 0) = veloc(2, k);
    pos(0, 0) = posit(0, k); pos(1, 0) = posit(1, k); pos(2, 0) = posit(2, k);
    std::cout << "acceleration for each iteration for each object" << std::endl;
        for(int i=0; i<N-1; i++){
            solar->calculateForcesAndEnergy();
            for (int j=0; j<dim; j++){
//            std::cout << bodies[k].force << std::endl;
                vec3 acc = bodies[k].force/bodies[k].mass;
                vel(j, i+1) = vel(j, i) + h*acc(j);
                pos(j, i+1) = pos(j, i) + h*vel(j,i);
                bodies[k].resetForce();
             }
          std::string filename = "euler_pos";
          std::string arg = std::to_string(k);
          filename.append(arg);
          filename.append(".txt");
          WriteToFile write;
          write.WritetoFile(filename, pos);
         }
}


  finish = std::clock();   // end timing
  double time_used = (double)(finish - start)/(CLOCKS_PER_SEC );
  std::cout << std::setprecision(10) << "Time used: " << time_used << " s at " << N/365 <<" yr" << std::endl;
  std::cout << "Write positions to file" << std::endl;
//  getPos(pos);
//  getVel(vel);
  return pos;
}
