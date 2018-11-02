/*
Compute the 2D ising model.
*/

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

ofstream ofile;
int main(int argc, char* argv[])
{
  string filename;
  int n_spin, MCs;      // number of spins and number of Monte carlo simulations
  double Temp0, Temp_final, Temp_step;
  if (argc <= 1){
    cout << "Bad usage: type filename for output" << endl;
    exit(1);
  }
  else{
    filename = argv[1];
  }
  n_spin = 2; MCs = 100; Temp0 = 1.2; Temp_final = 1.7; Temp_step = 0.05;

  string output = filename;
  string arg = to_string(n_spin);
  output.append(arg);
  output.append(".txt");
  ofile.open(output);

  Metropolis metroplis;
  System sys;
  for (double T = Temp0; T<=Temp_final; T+=Temp_step){
    vec ExpValues = zeros(5);

    metroplis.Metropolis(n_spin, MCs, T, ExpValues);

    sys.writefile(n_spin, MCs, T, ExpValues);

  }
  ofile.close();

  return 0;
}
