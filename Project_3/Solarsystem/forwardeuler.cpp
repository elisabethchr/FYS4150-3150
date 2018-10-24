//definitions for functions within class ForwardEuler
#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "forwardeuler.h"
#include "celestialobject.h"
#include "solarsystem.h"
#include "writetofile.h"

void ForwardEuler::Euler(double dt)// :m_dt(dt)
{
    m_dt = dt;
}

arma::mat ForwardEuler::Integrate(SolarSystem &input_system, double h, int n_years){
    /*
 * Compute the position of the planet using forward Euler method.
 */
    solar = &input_system;
    int  NumberOfObjects = solar->numberOfObjects();

    int N = 367*n_years;
    int dim = 3;

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
    * Forward Euler calculation
    */

    for(int k=0; k<NumberOfObjects; k++){
        //initializing to position and velociy matrices used in forward euler
        vel(0, 0) = vel_init(0, k); vel(1, 0) = vel_init(1, k); vel(2, 0) = vel_init(2, k);
        pos(0, 0) = pos_init(0, k); pos(1, 0) = pos_init(1, k); pos(2, 0) = pos_init(2, k);

        //write positions of each object/planet k to file
        std::string filename = "euler_pos";
        std::string arg = std::to_string(k);
        filename.append(arg);
        filename.append(".txt");
        WriteToFile write;

        for(int i=0; i<N-1; i++){
            bodies[k].resetForce();
            solar->calculateForcesAndEnergy();
            bodies[k].position += h*bodies[k].velocity;
            bodies[k].velocity += h*bodies[k].force/bodies[k].mass;

            for (int j=0; j<dim; j++){
                pos(j, i+1) = bodies[k].position(j);
            }
            write.WritetoFile(filename, pos(0,i+1), pos(1,i+1), pos(2,i+1));
        }

        std::cout << "size of pos for k = " << k << " is " << sizeof(pos) << std::endl;
    }

    //stop timing
    finish = std::clock();   // end timing
    double time_used = (double)(finish - start)/(CLOCKS_PER_SEC );
    std::cout << std::setprecision(10) << "Time used: " << time_used << " s at " << N/365 <<" yr" << std::endl;
    return pos;
}
