#ifndef STOCKMARKED_H
#define STOCKMARKED_H

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
#include <mpi.h>

using namespace std;
using namespace arma;

class StockMarked
{
  public:
    void Model(int Nagents, int transactions, string filename, double m0, vec outvalues);
    void initialize(int Nagents, double m0);
    void WriteToFile(int Nagents, vec mean_agents, string filename);
    vec Agents();

  private:
    vec agents;

};
#endif // STOCKMARKED_H
