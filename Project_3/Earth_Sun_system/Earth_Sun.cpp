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
  int year; int N; int dim = 3;
  string PrYear;
  if (argc <= 4){
    cout << "Not enough arguments, need output filename for Euler and verlet, N years and dimensions" << endl;
    exit(1);
  }
  else{
    filename1 = argv[1];
    filename2 = argv[2];
    year = atoi(argv[3]);
    PrYear = (argv[4]);

  }

  N = 365*year+2;
  cout << N << endl;

  double stepsPrYear = 365.25;
  double epsilon = 1e-5;

  cout << "Eulers method:" << endl;
  /*
  string dataEuler1 = filename1;

  //dataEuler1.append("yr.txt");
  Euler(N, dim, dataEuler1, 0.1);
  cout << "---------" << endl;
  */
  string dataEuler2 = filename1;

  dataEuler2.append(".txt");
  Euler(N, dim, dataEuler2, epsilon, 1);
  cout << "---------" << endl;
  /*
  string dataEuler3 = filename1;

  //dataEuler3.append("yr.txt");
  Euler(N, dim, dataEuler3, 10);
  cout << "---------" << endl;
*/
  // Verlet
  cout << "Velocity verlet method:" << endl;
  /*
  string dataVerlet1 = filename2;

  //dataVerlet1.append("yr.txt");
  VelocityVerlet(N, dim, dataVerlet1, 0.1);
  cout << "---------" << endl;
*/
  string dataVerlet2 = filename2;

  dataVerlet2.append(".txt");
  VelocityVerlet(N, dim, dataVerlet2, epsilon, 1);
  cout << "---------" << endl;
/*
  string dataVerlet3 = filename2;

  //dataVerlet3.append("yr.txt");
  VelocityVerlet(N, dim, dataVerlet3, 10);
  cout << "---------" << endl;
*/


  // Escape Velocity:
  string Vesc = filename2;
  Vesc.append("Vesc");
  //EscapeVelocity(N, dim, Vesc);


  // Alterative gravitational Force:
  string beta = filename2;
  beta.append("beta");
  AlterativeForce(N, dim);

  return 0;
}
