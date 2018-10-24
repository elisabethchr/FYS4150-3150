//definitions for functions within class ForwardEuler
#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <math.h>

#include "velocityverlet.h"
#include "celestialobject.h"
#include "solarsystem.h"
#include "writetofile.h"

void VelocityVerlet::Verlet(double dt)
{
    m_dt = dt;
}

void VelocityVerlet::Integrate(SolarSystem &input_system, double h, int n_years){
    /*
    Compute the position of the planet using forward Euler method.
    */
    solar = &input_system;
    int  NumberOfObjects = solar->numberOfObjects();

    int dim = 3;
    int N = n_years/h;
    double epsilon = 1e-5;

    //Obtain the initial velocities and positions for planets in m_objects;
    std::vector<CelestialObject> &bodies = solar->objects();
    arma::vec time = arma::zeros(N);
    //Start timing
    arma::vec t = arma::zeros(N);
    std::clock_t start, finish;
    start = std::clock();

    /*
    * Velocity Verlet calculation
    */

    for(int i=0; i<N-1; i++){
        arma::mat acc = arma::zeros(dim, NumberOfObjects);
        //solar->calculateForcesAndEnergy();
        solar->relativisticForcesAndEnergy();
        //std::cout <<"before planet"<<" " << bodies[1].position << std::endl;
        for(int k=0; k<NumberOfObjects; k++){
            double E=solar->totalEnergy();

            //write to file here
            std::string filename1 = "Verlet_pos";
            std::string filename2 = "Verlet_energies";
            std::string arg = std::to_string(k);
            std::string arg2 = std::to_string(h);
            std::string arg3 = std::to_string(n_years);
            //std::string arg4 = std::to_string(beta(b));
            filename1.append(arg);
            filename1.append("_dt");
            filename1.append(arg2);
            filename1.append("_yr");
            filename1.append(arg3);
            //filename.append("_beta");
            //filename.append(arg4);
            filename1.append(".txt");
            filename2.append(arg);
            filename2.append(".txt");
            WriteToFile write;
            write.WritetoFile(filename1, bodies[k].position(0), bodies[k].position(1), bodies[k].position(2));
            write.WritetoFile_Energy_AngMom(filename2, E, time(i));

            for (int j=0; j<dim; j++){
                acc(j,k) = bodies[k].force(j)/bodies[k].mass;
                bodies[k].position(j) += h*bodies[k].velocity(j) + (h*h*0.5)*acc(j,k);
            }

            bodies[k].resetForce();
        }

        //solar->calculateForcesAndEnergy(); // Update gravitational force
        solar->relativisticForcesAndEnergy();
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
                l0 = solar->angularMomentum();
                //std::cout << "Energy at start time = " << E0 << std::endl;
                std::cout << "angular momentum start = " << l0.length() << std::endl;
            }
            if (i == N-2){
                E_end = solar->totalEnergy();
                l_end = solar->angularMomentum();
                //std::cout << "Energy at end time = " << E_end << std::endl;
                std::cout << "angular momentum end = " << l_end.length() << std::endl;

                if (fabs(E0 - E_end)< epsilon){
                    std::cout <<"Energy is conserved for planet" << k <<", dE = " << E0 - E_end << std::endl;
                }
                else{
                    std::cout <<"Energy is not conserved for planet"<< k <<", dE = " << E0 - E_end<< std::endl;
                }
            }

            // testing angular momentum
            dl = (l_end - l0).length();
            if (i == N-2){
                if (fabs(dl) <= epsilon){
                    std::cout <<"Angular momentum is conserved for planet" << k << std::endl;
                    std:: cout << "Difference in angular momentum: dl = " << dl << std::endl;
                }
                else{
                    std::cout <<"Angular momentum is not conserved for planet" << k << std::endl;
                    std:: cout << "Difference in angular momentum: dl = " << dl << std::endl;
                }
            }

            bodies[k].resetForce();
        }
        //std::cout <<"after  planet" <<" " << bodies[2].position << std::endl;
    time(i+1) = time(i)+h;
    }

    //stop timing
    finish = std::clock();   // end timing
    double time_used = (double)(finish - start)/(CLOCKS_PER_SEC );
    std::cout << std::setprecision(10) << "Time used: " << time_used << " s at " << n_years <<" yr" << std::endl;
}
