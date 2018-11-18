#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
#include <mpi.h>

#include "Metropolis.h"
#include "System.h"

using namespace std;
using namespace arma;

void Metropolis::metropolis(int n_spin, int MCs, double Temp, vec ExpValues, string filename, int choise)
{
  random_device rd;
  mt19937_64 gen(rd());
  uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);   // Set up uniform distribution from 0 to 1
  System sys;


  double Ein; double Min;
  sys.initialize(n_spin, Temp, Ein, Min, choise);
  mat spin_matrix = sys.Lattice();

  double Energy = sys.Energy();
  double MagneticMoment = sys.MagneticMoment();
  vec w = zeros(2*fabs(Energy)+1);
  cout << "E0=" << Energy << " M0=" << MagneticMoment << " at T=" << Temp << endl;


  for (int dE=Energy; dE<=fabs(Energy); dE+=4){
    w(dE+fabs(Energy)) = exp(-dE/Temp);
  }

  // Monte Carlo
  for (int cycle=1; cycle<=MCs; cycle++){

    int counter = 0;
    for (int x=0; x<n_spin; x++){
      for (int y=0; y<n_spin; y++){
        int ix = (int) (RandomNumberGenerator(gen) * (double)n_spin);
        int iy = (int) (RandomNumberGenerator(gen) * (double)n_spin);

        int deltaE = 2*spin_matrix(iy,ix) * (spin_matrix(iy,sys.periodic(ix,n_spin,-1))
                                          + spin_matrix(sys.periodic(iy,n_spin,-1), ix)
                                          + spin_matrix(iy, sys.periodic(ix,n_spin, 1))
                                          + spin_matrix(sys.periodic(iy,n_spin,1), ix));

        // Metropolis algorithm
        if (RandomNumberGenerator(gen) <= w(deltaE+8)){
          spin_matrix(iy,ix) *= -1.0;      // Flip one spin
          MagneticMoment += (double) 2*spin_matrix(iy,ix);
          Energy += (double) deltaE;
          counter += 1;

        }
      }
    }

    ExpValues(0) += Energy;
    ExpValues(1) += Energy*Energy;
    ExpValues(2) += MagneticMoment;
    ExpValues(3) += MagneticMoment*MagneticMoment;
    ExpValues(4) += fabs(MagneticMoment);

    // write only every 100 value
    if (cycle >= MCs/10) {
      if ((cycle%100==0) && (cycle != 0)){
        //cout << "counter = " << counter << endl;
        sys.writefile(n_spin, MCs, Temp, ExpValues, filename, cycle, counter);
      }
    }
  }
}

void Metropolis::metropolisMPI(int n_spin, int myLoopBegin, int myLoopEnd, int MCs, double Temp,
                                vec ExpValues, string filename, int choise, int myRank)
{
  random_device rd;
  mt19937_64 gen(rd());
  uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);   // Set up uniform distribution from 0 to 1
  System sys;

  double Ein; double Min;
  sys.initialize(n_spin, Temp, Ein, Min, choise);
  mat spin_matrix = sys.Lattice();

  double Energy = sys.Energy();
  double MagneticMoment = sys.MagneticMoment();
  vec w = zeros(17);

  for (int dE=-8; dE<=8; dE+=4){
    w(dE+8) = exp(-dE/Temp);
  }

  // Monte Carlo
  for (int cycle=myLoopBegin; cycle<=myLoopEnd; cycle++){
    int counter = 0;
    for (int x=0; x<n_spin; x++){
      for (int y=0; y<n_spin; y++){
        int ix = (int) (RandomNumberGenerator(gen) * (double)n_spin);
        int iy = (int) (RandomNumberGenerator(gen) * (double)n_spin);

        int deltaE = 2*spin_matrix(iy,ix) * (spin_matrix(iy,sys.periodic(ix,n_spin,-1))
                                          + spin_matrix(sys.periodic(iy,n_spin,-1), ix)
                                          + spin_matrix(iy, sys.periodic(ix,n_spin, 1))
                                          + spin_matrix(sys.periodic(iy,n_spin,1), ix));

        // Metropolis algorithm
        if (RandomNumberGenerator(gen) <= w(deltaE+8)){
          spin_matrix(iy,ix) *= -1.0;      // Flip one spin
          MagneticMoment += (double) 2*spin_matrix(iy,ix);
          Energy += (double) deltaE;
          counter += 1;

        }
      }
    }
    //if (cycle >= MCs/5){

      ExpValues(0) += Energy;
      ExpValues(1) += Energy*Energy;
      ExpValues(2) += MagneticMoment;
      ExpValues(3) += MagneticMoment*MagneticMoment;
      ExpValues(4) += fabs(MagneticMoment);
    //}
  }

  // Mpi stuff
  vec totExpValues = zeros(5);
  for (int i=0; i<5; i++){
    // Gather the values from the cores
    MPI_Reduce(&ExpValues[i], &totExpValues[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  }

  if (myRank==0){
    // Write to file
    sys.writefileMPI(n_spin, MCs, Temp, totExpValues, filename);

  }

}
