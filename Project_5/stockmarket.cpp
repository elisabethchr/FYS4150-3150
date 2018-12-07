#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
//#include <mpi.h>

#include "stockmarked.h"

using namespace std;
using namespace arma;


vec StockMarked::Model(int Nagents, int transactions, double m0, vec agents, double lmbd, double alpha, double gamma)
{

    // Compute the exchange model for a trading process between two partners with an initial capital m0


    // Set up uniform distribution from 0 to 1
    unsigned __int64 i;
    i = __rdtsc();
    mt19937_64 gen(i);
    uniform_real_distribution<double> RandomNumberGenerator(0,1);

    // Initialize all agents with a startup capital of m0
    agents.fill(m0);

    // Stating parameters involved in the computations
    double var_m, exp_m, prev_exp_m;
    prev_exp_m = 1e10;
    double dm, mj, mi;
    double p_ij, c_ij, max_c_ij;
    double eps, r;
    int count = 0;

    // Number of transactions matrix
    mat c_matrix = zeros(Nagents, Nagents);

    // Run transactions
    for (int trans=1; trans<transactions; trans++)
    {
        // Pick a pair of agents at random
        int i = (int) (RandomNumberGenerator(gen) * (double) Nagents);
        int j = (int) (RandomNumberGenerator(gen) * (double) Nagents);
        eps = RandomNumberGenerator(gen);
        r = RandomNumberGenerator(gen);

        // If same agent being drawn, draw again
        if(i==j)
        {
            int i = (int) (RandomNumberGenerator(gen) * (double) Nagents);
            int j = (int) (RandomNumberGenerator(gen) * (double) Nagents);
            count ++;
        }

        /*if(trans == 1){
            cout << "eps = " << eps << "\n";
            cout << "rdtsc() = " << i << "\n";
            cout << "rd() = " << rd() << "\n";
            cout << "i: " << i << " j: " << j << "\n";
        }*/

        // Extracting wealth of each agent
        mi = agents(i);
        mj = agents(j);

        // The current transaction
        c_ij = c_matrix(i, j);

        // Calculate probabilities of interaction:
        if(mi != mj)
        {
            p_ij = pow(fabs(mi - mj), -alpha)*pow((c_ij + 1), gamma);
        }
        else
        {
            p_ij = 1;
        }

        // Compute the transaction if:
        // - the probabiltiy is greater than some random number r
        // - the two agents is not the same agent
        if (i!=j && p_ij>r)
        {
            // Update number of transactions made for each agent
            c_matrix(i, j) += 1;
            c_matrix(j, i) += 1;

            /*if(c_ij + 1 > max_c_ij)
            {
                max_c_ij = c_ij + 1;
            }*/

            // Update the wealth of each agent
            dm = (eps*mj - (1.0 - eps)*mi)*(1.0 - lmbd);
            agents(i) += dm;
            agents(j) -= dm;
        }

        // Find the average variance for each transaction
        var_m += var(agents);
        if (trans%10000 == 0)
        {
            exp_m = var_m/((double) trans);

            // Check if:
            // -change in average variance > 0.05% --> minimize the change by setting previous value equal to the current average variance
            // -change in average variance < 0.05% --> reached equilibrium state and break transactions loop
            if ((fabs(prev_exp_m - exp_m)/fabs(prev_exp_m)) > 0.005)
            {
                prev_exp_m = exp_m;
            }
            else
            {
                cout << "Equilibrium state is reached, at number of transaction " << trans << endl;
                break;
            }

            var_m = 0;
        }
    }
    //cout << count << endl;
    return agents;
}

void StockMarked::Simulation(int Nagents, int runs, int transactions, double m0, string filename, double lmbd, double alpha)
{

    // Run the simulations and write data to file

    vec mean_agents = zeros(Nagents);
    clock_t start, stop;
    start = clock();

    // Run the simulations for the transactions;
    for (int run=0; run<runs; run++){
        cout << "Run no. " << run << ": ";
        vec in_agents = zeros(Nagents);

        agents = Model(Nagents, transactions, m0, in_agents, lmbd, alpha);
        mean_agents += sort(agents);
    }
    // Find mean value of the agents
    mean_agents = mean_agents/(runs);


    // Write to file:
    //WriteToFile(Nagents, mean_agents, filename);

    stop = clock();
    double time_used = (double)(stop - start)/(CLOCKS_PER_SEC );
    cout << "Time used: t = " << time_used << " s." << endl;

}


void StockMarked::WriteToFile(int Nagents, vec mean_agents, string filename)
{

    // Writing data to file

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

double StockMarked::Probability()
{
    return p_ij;
}
