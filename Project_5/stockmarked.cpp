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

//void StockMarked::Model(int Nagents, int transactions, double m0, vec agents)
vec StockMarked::Model(int Nagents, int transactions, double m0, vec agents)
{
  // Set up uniform distribution from 0 to 1
  random_device rd;
  mt19937_64 gen(rd());
  uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);

  // initialize vector with the agents with startup capital of m0
  agents.fill(m0);

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
    var_m = var(agents);
    if (trans%10000 == 0){
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

  return agents;
}

void StockMarked::Simulation(int Nagents, int runs, int transactions, double m0, string filename)
{
  vec mean_agents = zeros(Nagents);
  clock_t start, stop;
  start = clock();
  for (int i=0; i<runs; i++){
    // Metropolis/Monte carlo stuff:
    vec agents = zeros(Nagents);

    agents = Model(Nagents, transactions, m0, agents);

    mean_agents += (agents);

  }

  mean_agents = mean_agents/(runs);       // Find mean value of the agents

  // write to file:
  WriteToFile(Nagents, mean_agents, filename);

  stop = clock();
  double time_used = (double)(stop - start)/(CLOCKS_PER_SEC );
  cout << "Time used: t = " << time_used << " s." << endl;

}


void StockMarked::WriteToFile(int Nagents, vec mean_agents, string filename)
{
  ofstream mfile;
  filename.append(".txt");
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
