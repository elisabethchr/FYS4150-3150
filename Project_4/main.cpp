#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>

#include "Metropolis.h"
#include "System.h"

using namespace std;
using namespace arma;

int main(int argc, char* argv[])
{
    string filename;
    int n_spin, MCs;      // number of spins and number of Monte carlo simulations+
    double Temp0, Temp_final, Temp_step;

    //  if (argc <= 1){
    //    cout << "Bad usage: type filename for output" << endl;
    //    exit(1);
    //  }
    //  else{
    //    filename = argv[1];
    //  }

    filename = "Ising";
    n_spin = 2; MCs = 1000000; Temp0 = 1.; Temp_final = 1.; Temp_step = 0.05;
    int nTemp = (Temp_final-Temp0)/Temp_step;

    Metropolis metrop;
    System sys;

    //    for (double T = Temp0; T<=Temp_final; T+=Temp_step){
    vec ExpValues = zeros(5);
    //metrop.metropolis(n_spin, MCs, T, ExpValues, nTemp, filename);
    metrop.metropolis(n_spin, MCs, 1, ExpValues, 1, filename);

    //sys.writefile(n_spin, MCs, T, ExpValues, filename);
    //    }

    //ofile.close();
    return 0;
}

