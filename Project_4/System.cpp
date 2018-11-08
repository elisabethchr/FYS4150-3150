#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
#include <string>

#include "System.h"
#include "Metropolis.h"

using namespace  std;
using namespace arma;

// output file
ofstream ofile;
void System::writefile(int n_spin, int MCs, double Temp, vec expValues, int cycle, int nTemp, string filename)
{
//expValues = (E, E^2, M, M^2, |M|)
    double norm = 1.0/((double) cycle);
    double averageE = expValues(0)*norm;
    double averageE2 = expValues(1)*norm;
    double averageM = expValues(2)*norm;
    double averageM2 = expValues(3)*norm;
    double averageMabs = expValues(4)*norm;
    //std::cout << "<E^2> = " << averageE2 << " E*E = " << averageE*averageE << std::endl;
    //double varE = (averageE2 - averageE*averageE)/n_spin/n_spin;
    //double varM = (averageM2 - averageMabs*averageMabs)/n_spin/n_spin;
    double Cv = (averageE2 - averageE*averageE)/(Temp*Temp*n_spin*n_spin);// /n_spin/n_spin;
    double chi = (averageM2 - averageM*averageM)/(Temp*n_spin*n_spin); // /n_spin/n_spin;

    string output = filename;
    string arg = to_string(n_spin);
    string NumberOfMCs = to_string(MCs);
    string NumberOfTemps = to_string(nTemp);
    output.append("_nSpin_");
    output.append(arg);
    output.append("_nTemp_");
    output.append(NumberOfTemps);
    output.append("_MC_");
    output.append(NumberOfMCs);
    output.append("_.txt");
    //ofile.open(output, ios::trunc);
    ofile.open(output, ios::app | ios::out);

    ofile << setiosflags(ios::showpoint | ios::uppercase);
//    ofile << setw(15) << "T" << setw(15) << "<E>" << setw(15) << "Cv" << setw(15) << "<M>" << setw(15) << "Chi" << setw(15) << "|<M>|";
    ofile << setw(20) << setprecision(8) << Temp;       //print temperatures
    ofile << setw(20) << setprecision(8) << averageE/n_spin/n_spin;  //print mean energies
    ofile << setw(20) << setprecision(8) << Cv/Temp/Temp;   //print specific heat capacities (energy variance)
    ofile << setw(20) << setprecision(8) << averageM/n_spin/n_spin; //print average magnetic moment
    ofile << setw(20) << setprecision(8) << chi/Temp;   //print susceptibility (magnetic variance)
    ofile << setw(20) << setprecision(8) << averageMabs/n_spin/n_spin;  //print absolute value of magnetic moment
    ofile << setw(20) << setprecision(8) << cycle << "\n";    //print Monte Carlo cycle
    ofile.close();
}


void System::initialize(int n_spin, double Temp, double& E1, double& M1)
{
    /*
 * Initialize lattice matrices containing spins
 */
    E = E1; M = M1;
    spin_matrix.zeros(n_spin,n_spin);
    for (int y=0; y<n_spin; y++){
        for (int x=0; x<n_spin; x++){
            spin_matrix(y,x) = 1.0;
            M += (double) spin_matrix(y,x);
        }
    }
    for (int y=0; y<n_spin; y++){
        for (int x=0; x<n_spin; x++){
            E -= (double) spin_matrix(y,x)*(spin_matrix(periodic(y,n_spin,-1), x) + spin_matrix(y,periodic(x,n_spin,-1)));
        }
    }
}

int System::periodic(int i, int n_spin, int add)
{
    return (i+n_spin+add) % (n_spin);       //integer remainder
}

arma::mat System::Lattice(){
    return spin_matrix;
}

double System::Energy(){
    return E;
}

double System::MagneticMoment(){
    return M;
}
