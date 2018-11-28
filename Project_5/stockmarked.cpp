#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
#include <mpi.h>

#include "stockmarked.h"

using namespace std;
using namespace arma;

void StockMarked::Model(int Nagents, int transactions, string filename, double m0, vec outvalues)
{
  // Set up uniform distribution from 0 to 1
  random_device rd;
  mt19937_64 gen(rd());
  uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);

  // initialize
  initialize(Nagents, m0);
  vec agents = Agents();
  vec agents_post = zeros(Nagents);
  vec Gibbs = zeros(Nagents);
  double beta = 1.0/m0;
  double dm = 0.01;
  double var_m, exp_m, prev_exp_m;
  prev_exp_m = 1e8; // some big number.


  // run transactions
  for (int trans=0; trans<transactions; trans++){
    int i = (int) (RandomNumberGenerator(gen) * (double) Nagents);
    int j = (int) (RandomNumberGenerator(gen) * (double) Nagents);
    double eps = (double) (RandomNumberGenerator(gen));

    double mi = agents(i);
    double mj = agents(j);

    double dm = eps*mj - (1.0-eps)*mi;

    agents(i) += dm;
    agents(j) -= dm;

    // find the variance for each transaction, and check if equilibrium is reached
    var_m = var(agents);    // why variance???
    if (trans%1000 == 0){
      exp_m = var_m/trans;

      // check for equilibrium state
      if ((fabs(prev_exp_m - exp_m)/fabs(prev_exp_m)) > 0.01){
        prev_exp_m = exp_m;
      }

      else{
        cout << "Equilibrium state is reached, at number of transaction " << trans << endl;
        break;
      }

      var_m = 0;

    }
  }
}

void StockMarked::initialize(int Nagents, double m0)
{
  vec agents = zeros(Nagents);
  for (int i=0; i<Nagents; i++){
    agents(i) = m0;
  }
}

void StockMarked::WriteToFile(int Nagents, vec mean_agents, string filename)
{
  ofstream mfile;
  mfile.open(filename, ios::app | ios::out);
  for (int i=0; i<Nagents; i++){
    mfile << setw(15) << setprecision(8) << mean_agents(i) << "\n";
  }
  mfile.close();
}

vec StockMarked::Agents()
{
  return agents;
}
