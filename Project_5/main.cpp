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
    double m0, lmbd, alpha, gamma;
    if (argc <= 1){
        cout << "Bad usage: \n";
        cout << "type filename for output, inital capital, saving parameter [0,1), wealth parameter [0, 1), and transactions parmeter [0, 1)" << endl;
        exit(1);
    }
    else{
        filename = argv[1];
        m0 = atof(argv[2]);                      // Initial money
        lmbd = atof(argv[3]);                   // Saving parameter
        alpha = atof(argv[4]);                  // Wealth parameter
        gamma = atof(argv[5]);                  // Number of transactions parameter
    }

    int Nagents = 1000;                    // Number of agents
    int runs = 1000;                       // Number of simulations
    int transactions = 1e7;               // Number of MC cycles
    vec mean_agents = zeros(Nagents);

    //string filename = "data/Task_d/d_stockmarket";
    string NumberOfAgents = to_string(Nagents);
    string lambda = to_string(lmbd);
    string alph = to_string(alpha);
    string gam = to_string(gamma);
    filename.append("_NAgents_");
    filename.append(NumberOfAgents);
    filename.append("_lmbd_");
    filename.append(lambda);
    filename.append("_a_");
    filename.append(alph);
    filename.append("_g_");
    filename.append(gam);
    cout << filename << endl;

    cout << "======================================================" << endl;
    StockMarked SM;
//Simulation(int Nagents, int runs, int transactions, double m0, string filename, double lmbd, double alpha, double gamma)
    SM.Simulation(Nagents, runs, transactions, m0, filename, lmbd, alpha, gamma);

    cout << "======================================================" << endl;

    // End main function
    return 0;
}
