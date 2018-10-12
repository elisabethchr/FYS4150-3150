#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "Method_Earth_sun.hpp"


using namespace std;
using namespace arma;


int main (int argc, char* argv[]){
  string filename1; string filename2;
  int year; int N; int dim;
  if (argc <= 4){
    cout << "Not enough arguments, need output filename for Euler and verlet, N years and dimensions" << endl;
    exit(1);
  }
  else{
    filename1 = argv[1];
    filename2 = argv[2];
    year = atoi(argv[3]);
    dim = 3; //atoi(argv[4]);

  }
  //if (year >= 4){
  N = 365*year;// + (int)0.26*4;
  cout << N << endl;
  //}
  //else {N = 365*year;}
  string nyear = to_string(year);
  cout << "Eulers method:" << endl;
  string dataEuler = filename1;

  dataEuler.append(nyear);
  dataEuler.append("yr.txt");
  //mat velE = zeros(dim,N); mat posE = zeros(dim,N);mat accE = zeros(dim,N);
  Euler(N, dim, dataEuler);
  cout << "---------" << endl;


  cout << "Velocity verlet method:" << endl;
  string dataVerlet = filename2;
  dataVerlet.append(nyear);
  dataVerlet.append("yr.txt");
  //mat velV = zeros(dim, N); mat posV = zeros(dim, N);mat accV = zeros(dim,N);
  VelocityVerlet(N, dim, dataVerlet);


  return 0;
}
