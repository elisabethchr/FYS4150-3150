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

void Metropolis::Metropolis(int n_spin, int MCs, double Temp, vec ExpValues)
{
  double Energy = 0.0; double MagneticMoment = 0.0;
  mat spin_matrix = zeros(n_spin, n_spin);
  vec w = zeros(17);
  System sys;
  sys.initialize(n_spin, Temp, spin_matrix, Energy, MagneticMoment);
  for (int dE=-8; dE<=8;dE+=4){
    w(dE+8) = exp(-dE/Temp);
  }
  // Monte Carlo
  for (int y=0; y<n_spin; y++){
    for (int x=0; x<n_spin; x++){
      int ix = (int) (ran1(MCs)*(double)n_spin);
      int iy = (int) (ran1(MCs)*(double)n_spin);

      int deltaE = 2*spin_matrix(iy,ix) * (spin_matrix(iy,periodic(ix,n_spin,-1))
                  + spin_matrix(periodic(iy,n_spin,-1),ix) + spin_matrix(iy,periodic(ix,n_spin,1))
                  + spin_matrix(periodic(iy,n_spin,1),ix));
      if (ran1(MCs) <= w(deltaE+8)){
        spin_matrix(iy,ix) *= -1;      // Flip one spin
        M += (double) 2*spin_matrix(iy,ix);
        E += (double) deltaE;
      }
    }
  }
}
