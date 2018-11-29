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
  int exponent = 7; double m0; double lmbd;
  if (argc <= 1){
    cout << "Bad usage: \n";
    cout << "type filename for output, start money and saving parameter (0.25,0.5,0.9)" << endl;
    exit(1);
  }
  else{
    filename = argv[1];
    m0 = atof(argv[2]);                      // Initial money
    lmbd = atof(argv[3]);

  }

  int Nagents = 500;                    // Number of agents
  int runs = 1e3;                       // Experiments
  int transactions = 1e7;               // Number of MC cycles
  vec mean_agents = zeros(Nagents);


  cout << "======================================================" << endl;
  StockMarked SM;

  SM.Simulation(Nagents, runs, transactions, m0, filename, lmbd);

  cout << "======================================================" << endl;


  // End main function
  return 0;
}
