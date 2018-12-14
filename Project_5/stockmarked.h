#ifndef STOCKMARKED_H
#define STOCKMARKED_H

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
//#include <mpi.h>

using namespace std;
using namespace arma;

class StockMarked
{
  public:
    //void Model(int Nagents, int transactions, double m0, vec agents);
    vec Model(int Nagents, int transactions, double m0, vec agents, double lmbd, double alpha, double gamma);
    void Simulation(int Nagents, int runs, int transactions, double m0, string filename, double lmbd, double alpha, double gamma);
    void WriteToFile(int Nagents, vec mean_agents, string filename);

    vec Agents();
    double Probability();

  private:
    vec agents;
    vec outvalues;
    vec prob;
    double p_ij;

};
#endif // STOCKMARKED_H
