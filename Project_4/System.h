#ifndef SYSTEM_H
#define SYSTEM_H

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
#include <string>

using namespace arma;

class System
{
public:
  int periodic(int i, int n_spin, int add);
  void writefile(int n_spin, int mcs, double Temp, vec average, std::string filename, int cycles, int N_accepted);
  void writefileMPI(int n_spin, int nExperiments, double Temp, vec ExpValues, std::string filename);
  void initialize(int n_spin, double Temp, double Ein, double Min, int choise);//, double &E, double &M);
  void ComputeAverage(int mcs, int n_spin, vec &ExpValues, double Temp);

  mat Lattice();
  double MagneticMoment();
  double Energy();
  double meanMag();
  double meanMag2();
  double meanEnergy();
  double meanMagAbs();
  double meanEnergy2();
  double HeatCapacity();
  double Suseptibility();



private:
  mat spin_matrix;
  double E;
  double M;

  double averageE;
  double averageE2;
  double averageM;
  double averageM2;
  double averageMabs;
  double varE;
  double varM;
  double Cv;
  double chi;
};
#endif
