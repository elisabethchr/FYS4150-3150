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
  int exponent = 7; double m0; double lmbd; double alpha;
  if (argc <= 1){
    cout << "Bad usage: \n";
    cout << "type filename for output, start money, saving parameter (0.25,0.5,0.9), alpha (0.5,1.0,1.5,2.0)" << endl;
    exit(1);
  }
  else{
    filename = argv[1];
    m0 = atof(argv[2]);                      // Initial money
    lmbd = atof(argv[3]);
    alpha = atof(argv[4]);
  }

  int Nagents = 500;                    // Number of agents
  int runs = 10;                       // Experiments
  int transactions = 1e7;               // Number of MC cycles
  vec mean_agents = zeros(Nagents);


  cout << "======================================================" << endl;
  StockMarked SM;

  SM.Simulation(Nagents, runs, transactions, m0, filename, lmbd, alpha);

  cout << "======================================================" << endl;


  // End main function
  return 0;
}
