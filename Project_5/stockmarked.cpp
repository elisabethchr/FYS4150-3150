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
vec StockMarked::Model(int Nagents, int transactions, double m0, vec agents, double lmbd, double alpha)
{
  // Set up uniform distribution from 0 to 1
  random_device rd;
  mt19937_64 gen(rd());
  uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);

  // initialize vector with the agents with startup capital of m0
  agents.fill(m0);

  double beta = 1.0/m0;
  double dm = 0.01;
  double var_m, exp_m, prev_exp_m;
  prev_exp_m = 1e8; // some big number.

  // run transactions
  //int j;
  for (int trans=1; trans<transactions; trans++){

    int i = (int) (RandomNumberGenerator(gen) * (double) Nagents);
    int j = (int) (RandomNumberGenerator(gen) * (double) Nagents);
    double eps = (double) (RandomNumberGenerator(gen));

    double mi = agents(i);
    double mj = agents(j);
    double randomnumber = (double) (RandomNumberGenerator(gen));

    if (trans < 2){
      cout <<trans<< ": "<< i << " " << j << " - "<< mi << " " << mj << endl;
    }

    double dm = (eps*mj - (1.0 - eps)*mi)*(1.0 - lmbd);

    // Check if the agents are of the same transaction interest
    //if (pow(fabs(dm), -alpha) < randomnumber){

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
      } // end if test 1
    //}
  }
  return agents;
}

void StockMarked::Simulation(int Nagents, int runs, int transactions, double m0, string filename, double lmbd, double alpha)
{
  vec mean_agents = zeros(Nagents);
  clock_t start, stop;
  start = clock();
  for (int i=0; i<runs; i++){
    // Metropolis/Monte carlo stuff:
    cout << "Run no. " << i<< ": ";
    vec agents = zeros(Nagents);

    agents = Model(Nagents, transactions, m0, agents, lmbd, alpha);

    mean_agents += sort(agents);      // sort or not??

  }

  mean_agents = mean_agents/(runs);       // Find mean value of the agents
  //vec Gibbs = beta*exp(-beta*mean_agents);
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
