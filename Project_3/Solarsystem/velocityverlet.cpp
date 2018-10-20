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

arma::mat VelocityVerlet::Integrate(SolarSystem &input_system){//vec3 position, vec3 velocity, double mass){ // int N, int dim, std::string obj){
    /*
 * Compute the position of the planet using forward Euler method.
 */
    solar = &input_system;
    int  NumberOfObjects = solar->numberOfObjects();

    int N = 367*2;
    int dim = 3;
    double h = 1.0/((double) 365.25); //1.0/((double) 365.25);
    double hh2 = h*h*0.5;
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

        std::cout <<"before planet1"<<" " << bodies[1].position << std::endl;
        for(int k=0; k<NumberOfObjects; k++){

            //write to file here
            std::string filename = "Verlet_pos";
            std::string arg = std::to_string(k);
            filename.append(arg);
            filename.append(".txt");
            WriteToFile write;
            write.WritetoFile(filename, bodies[k].position(0), bodies[k].position(1), bodies[k].position(2));

            for (int j=0; j<dim; j++){
                acc(j,k) = bodies[k].force(j)/bodies[k].mass;
                bodies[k].position(j) += h*bodies[k].velocity(j) + (h*h*0.5)*acc(j,k);
            }

            bodies[k].resetForce();
        }

        solar->calculateForcesAndEnergy(); // Update gravitational force
        for (int k=0; k<NumberOfObjects; k++){
            vec3 acc_next = bodies[k].force/bodies[k].mass;

            for (int j=0; j<dim; j++){
                bodies[k].velocity(j) += (0.5*h)*(acc_next(j) + acc(j,k));
            }
            bodies[k].resetForce();

            // Put positions into a 3d matrix for further
            //for (int j=0; j<3; j++){
            //  pos(j,i,k) = bodies[k].position(j);
            //}
        }
        std::cout <<"after  planet1" <<" " << bodies[1].position << std::endl;

    }
    /*
    // Write to file

    for (int k=0; k<NumberOfObjects; k++){
      std::cout << "Plantet"<< k << std::endl;

      string filename = "Verlet_pos";
      string arg = to_string(k);
      filename.append(arg);
      filename.append(".txt")

      for (int i=0; i<N; i++){
        std::cout << pos(0,i,k) << "  "<< pos(1,i,k)<< "  "<< pos(1,i,k) << std::endl;
      }


    }
    */
    //stop timing
    finish = std::clock();   // end timing
    double time_used = (double)(finish - start)/(CLOCKS_PER_SEC );
    //  std::cout << std::setprecision(10) << "Time used: " << time_used << " s at " << N/365 <<" yr" << std::endl;
    return pos;
}
