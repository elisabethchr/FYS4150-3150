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
    //void Model(int Nagents, int transactions, double m0, vec agents);
    vec Model(int Nagents, int transactions, double m0, vec agents, double lmbd);
    void Simulation(int Nagents, int runs, int transactions, double m0, string filename, double lmbd);
    void WriteToFile(int Nagents, vec mean_agents, string filename);

    vec Agents();

  private:
    vec agents;
    vec outvalues;

};
#endif // STOCKMARKED_H
