// Economy project:

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>

#include "stockmarked.h"

using namespace std;

int main(int argc, char* argv[])
{
  string filename;
  int exponent = 7; double m0;
  if (argc <= 1){
    cout << "Bad usage: \n";
    cout << "type filename for output, start money" << endl;
    exit(1);
  }
  else{
    filename = argv[1];
    m0 = atof(argv[2]);                      // Initial money

  }

  int Nagents = 500;
  int runs = 1;             // Experiments                    // Number of agents
  int transactions = pow(10, exponent);       // Number of MC cycles
  vec mean_agents = zeros(Nagents);

  cout << "======================================================" << endl;
  StockMarked SM;

  SM.Simulation(Nagents, runs, transactions, m0, filename);
  /*
  clock_t start, stop;
  start = clock();
  for (int i=1; i<runs; i++){
    // Metropolis/Monte carlo stuff:
    vec outvalues = zeros(Nagents);

    SM.Model(Nagents, transactions, filename, m0, outvalues);

    mean_agents += sort(outvalues);

  }
  mean_agents.print("Mean agents: ");
  mean_agents = mean_agents/runs;

  // write to file:
  //SM.WriteToFile(Nagents, mean_agents, filename);

  stop = clock();
  double time_used = (double)(stop - start)/(CLOCKS_PER_SEC );
  cout << "Time used: t = " << time_used << " s." << endl;
  */
  cout << "======================================================" << endl;


  // End main function
  return 0;
}
