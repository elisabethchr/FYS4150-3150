#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>

#include "Metropolis.h"
#include "System.h"

using namespace std;
using namespace arma;

void Metropolis::metropolis(int n_spin, int MCs, double Temp, vec ExpValues, int nTemp, string filename)
{
    random_device rd;
    mt19937_64 gen(rd());     //Mersenne-Twister random number generator -> period = 2^19937
    uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);   // Set up uniform distribution from 0 to 1

    System sys;
    double E = 0.; double M = 0.;
    sys.initialize(n_spin, Temp, E, M);
    double Energy = sys.Energy();
    double MagneticMoment = sys.MagneticMoment();
    mat spin_matrix = sys.Lattice();
    vec w = zeros(pow(2, n_spin*n_spin)+1);

    //initialize energy difference matrix
    for (int dE=-8; dE<=8;dE+=4){
        w(dE+8) = exp(-dE/Temp);
    }

    /*
  Monte Carlo algorithm
 */
    int counter = 0;
    for (int cycle=1; cycle <=MCs; cycle++){
        for (int y=0; y<n_spin; y++){
            for (int x=0; x<n_spin; x++){
                int ix = (int) (RandomNumberGenerator(gen) * (double)n_spin);
                int iy = (int) (RandomNumberGenerator(gen) * (double)n_spin);
                int deltaE = 2*spin_matrix(iy,ix) * (spin_matrix(iy,sys.periodic(ix,n_spin,-1))
                                                     + spin_matrix(sys.periodic(iy,n_spin,-1), ix)
                                                     + spin_matrix(iy, sys.periodic(ix,n_spin, 1))
                                                     + spin_matrix(sys.periodic(iy,n_spin,1), ix));

                //cout << "RNG= "<<RandomNumberGenerator(gen) << " w = " << w(deltaE+8)<< " dE ="<<deltaE<< endl;
                if (RandomNumberGenerator(gen) <= w(deltaE+8)){
                    spin_matrix(iy,ix) *= -1;      // Flip one spin
                    MagneticMoment += (double) 2*spin_matrix(iy,ix);
                    Energy += (double) deltaE;
                    counter +=1;
                    //cout << " E= "<< Energy << " RNG= "<<RandomNumberGenerator(gen) << " w = " << w(deltaE+8)<< endl;
                    //cout << " cycle: "<< cycle<<" dE = " << deltaE << endl;
                }
            }
        }
        //cout << "Mc "<< cycle << " E = " << Energy << endl;
        ExpValues(0) += Energy;
        ExpValues(1) += Energy*Energy;
        ExpValues(2) += MagneticMoment;
        ExpValues(3) += MagneticMoment*MagneticMoment;
        ExpValues(4) += fabs(MagneticMoment);
        //ExpValues.print(" ");
        //        sys.writefile(n_spin, MCs, Temp, ExpValues, cycle, nTemp, filename);
        //if ((cycle%100==0) && (cycle != 0)){
            //cout << "i = " << cycle << endl;
            sys.writefile(n_spin, MCs, Temp, ExpValues, cycle, nTemp, filename);
        //  }

    }
    cout << "Number of accepted runs: " << counter << endl;
    //ofile.close();
}
