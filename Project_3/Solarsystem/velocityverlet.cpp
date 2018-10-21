//definitions for functions within class ForwardEuler
#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <math.h>
//#include "EnergyTest.hpp"
#include "velocityverlet.h"
#include "gravitationalforce.h"
#include "celestialobject.h"
#include "solarsystem.h"
#include "writetofile.h"

void VelocityVerlet::Verlet(double dt)// :m_dt(dt)
{
    m_dt = dt;
}

arma::mat VelocityVerlet::Integrate(SolarSystem &input_system){
    /*
 * Compute the position of the planet using forward Euler method.
 */
    solar = &input_system;
    int  NumberOfObjects = solar->numberOfObjects();

    int N = 367*2;
    int dim = 3;
    double h = 1.0/((double) 365.25); //1.0/((double) 365.25);
    double hh2 = h*h*0.5; double epsilon = 1e-5;
    arma::mat vel = arma::zeros(dim, N);
    arma::mat pos = arma::zeros(dim, N);


    //Obtain the initial velocities and positions for planets in m_objects;
    std::vector<CelestialObject> &bodies = solar->objects();
    arma::mat pos_init = arma::zeros(dim, bodies.size());
    arma::mat vel_init = arma::zeros(dim, bodies.size());
    arma::vec mass = arma::zeros(bodies.size());

    for (int i=0; i<bodies.size(); i++){
        CelestialObject &obj = bodies[i];
        for (int j=0; j < dim; j++){
            vel_init(j, i) = obj.velocity(j);
        }
        for(int j=0; j < dim; j++){
            pos_init(j, i) = obj.position(j);
        }
        mass(i) = obj.mass;
    }


    //Start timing
    arma::vec t = arma::zeros(N);
    std::clock_t start, finish;
    start = std::clock();

    /*
    * Velocity Verlet calculation
    */

    //arma::cube pos = arma::zeros(dim, N, NumberOfObjects);
    for(int i=0; i<N-1; i++){
      arma::mat acc = arma::zeros(dim, NumberOfObjects);
      solar->calculateForcesAndEnergy();

      //std::cout <<"before planet"<<" " << bodies[1].position << std::endl;
      for(int k=0; k<NumberOfObjects; k++){
        if (i==0){
          //initializing to position and velociy matrices
          vel(0, 0) = vel_init(0, k); vel(1, 0) = vel_init(1, k); vel(2, 0) = vel_init(2, k);
          pos(0, 0) = pos_init(0, k); pos(1, 0) = pos_init(1, k); pos(2, 0) = pos_init(2, k);
        }

        for (int j=0; j<dim; j++){
          acc(j,k) = bodies[k].force(j)/bodies[k].mass;
          bodies[k].position(j) += h*bodies[k].velocity(j) + (h*h*0.5)*acc(j,k);
        }
        //write to file here
        std::string filename = "Verlet_pos";
        std::string arg = std::to_string(k);
        std::string arg2 = std::to_string(h);
        std::string arg3 = std::to_string((int) round(N/365.25));
        filename.append(arg);
        filename.append("_dt");
        filename.append(arg2);
        filename.append("_yr");
        filename.append(arg3);
        filename.append(".txt");
        WriteToFile write;
        //write.WritetoFile(filename, bodies[k].position(0), bodies[k].position(1), bodies[k].position(2));

        bodies[k].resetForce();
      }

      solar->calculateForcesAndEnergy(); // Update gravitational force
      for (int k=0; k<NumberOfObjects; k++){
        vec3 acc_next = bodies[k].force/bodies[k].mass;
        for (int j=0; j<dim; j++){
          bodies[k].velocity(j) += (0.5*h)*(acc_next(j) + acc(j,k));
        }

        // Energy testing and angular momentum:
        double E0, E_end, dl;
        vec3 l0, l_end;
        if (i == 0){
          E0 = solar->totalEnergy();
          //l0 = solar->angularMomentum();
          //std::cout << "Energy at start time = " << E0 << std::endl;
          //l0.print("angular Momentum start:");
        }
        if (i == N-2){
          E_end = solar->totalEnergy();
          //l_end = solar->angularMomentum();
          //std::cout << "Energy at end time = " << E_end << std::endl;
          //l_end.print("angular momentum end:");
          if (fabs(E0 - E_end)< epsilon){
            std::cout <<"Energy is conserved for planet" << k <<", dE = " << E0 - E_end << std::endl;
          }
          else{
            std::cout <<"Energy is not conserved for planet"<< k <<", dE = " << E0 - E_end<< std::endl;
          }
        }
        // angularMomentum
        l0 = solar->angularMomentum();
        l_end = solar->angularMomentum();
        //l0.print("angular Momentum start:");
        //l_end.print("angular Momentum end:");
        dl = (l_end - l0).length();
        //std::cout <<dl << std::endl;

        if (fabs(dl) <= epsilon){
          std::cout <<"Angular momentum is conserved for planet" << k << std::endl;
        }
        else{
          std::cout <<"Angular momentum is not conserved for planet" << k << std::endl;
          std:: cout << "Difference in angular momentum: dl = " << dl << std::endl;
        }

        bodies[k].resetForce();
      }
      //std::cout <<"after  planet" <<" " << bodies[2].position << std::endl;
    }

    //stop timing
    finish = std::clock();   // end timing
    double time_used = (double)(finish - start)/(CLOCKS_PER_SEC );
    std::cout << std::setprecision(10) << "Time used: " << time_used << " s at " << N/365.25 <<" yr" << std::endl;
    return pos;
}
