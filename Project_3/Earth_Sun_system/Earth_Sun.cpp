//main function for Earth_Sun system --> this we can push to the SolarSystem when we have a generalized code for all planets
#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "Method_Earth_sun.hpp"
#include "forwardeuler.h"
#include "velocityverlet.h"
#include "solarsystem.h"

using namespace std;
using namespace arma;


int main (int argc, char* argv[]){
  string filename1; string filename2;
  vec filenames;
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
//    vector<string> filenames;
 //   int i, n;
//    char* filenames;  //array of pointers to strings, i.e. can store multiple elements, char filenames[argc] cannot (older compilers won't like this line)
 //   filenames = (char*) malloc(i+1);
    //char**command = malloc(argc, sizeof(char*)); //older compilers would likt this
/*
    for (int i=0; i<2; i++)
    {
        filenames(i) = strcpy(filename, argv[1]);
    }
*/
//    filenames.push_back(filename1); filenames.push_back(filename2);
  }
  cout << filenames << endl;
  N = 365*year+2;
  cout << N << endl;

  double stepsPrYear = 365.25;
  double epsilon = 1e-5;
  double dt = 0.001;
  mat pos_euler;
  mat pos_verlet;

//  SolarSystem solarsystem;
// CelestialObject &sun = solarsystem.createCelestialBody( vec3(0,0,0), vec3(0,0,0), 1.0 );

  class ForwardEuler integrate_euler;     //need object of the class (integrate_euler) to use the member function Integrate within the class
  pos_euler = integrate_euler.Integrate(N, dim, filename1, epsilon, dt);
  cout << pos_euler << endl;

  class VelocityVerlet integrate_verlet;  //need object of the class (integrate_verlet) to use the member function Integrate within the class
  //pos_verlet =
  integrate_verlet.Integrate(N, dim, filename1, epsilon, dt);



//for (int i=0, i<filenames.size());
//  for(int i<0; i<filenames.size(), i++)
 // {
/*
  FILE* myfile = fopen(filename1.c_str(), "w");
  double x; double y;  double z;
  for(int i=0; i<N-1; i++){
      for (int j=0; j < dim ; j++){
          if(j==0){
          x = A(j,i);
          }
          else if(j==1){
          y = A(j, i);
          }
          else if(j==2){
          z = A(j, i);
          }
  fprintf(myfile, "%.8f     %.8f     %.8f\n", x,y,z);
}
}
fclose(myfile);
//}
*/

//Euler
  /*
  string dataEuler1 = filename1;
  //dataEuler1.append("yr.txt");
  Euler(N, dim, dataEuler1, 0.1);
  cout << "---------" << endl;
  */
//  string dataEuler2 = filename1;


/*
  dataEuler2.append(".txt");
  Euler(N, dim, dataEuler2, epsilon, 1);
  cout << "---------" << endl;
*/
  /*
  string dataEuler3 = filename1;
  //dataEuler3.append("yr.txt");
  Euler(N, dim, dataEuler3, 10);
  cout << "---------" << endl;
*/



  // Verlet
//  cout << "Velocity verlet method:" << endl;
  /*
  string dataVerlet1 = filename2;
  //dataVerlet1.append("yr.txt");
  VelocityVerlet(N, dim, dataVerlet1, 0.1);
  cout << "---------" << endl;
*/

/*
  string dataVerlet2 = filename2;

  dataVerlet2.append(".txt");
  VelocityVerlet(N, dim, dataVerlet2, epsilon, 1);
  cout << "---------" << endl;
*/
/*
  string dataVerlet3 = filename2;
  //dataVerlet3.append("yr.txt");
  VelocityVerlet(N, dim, dataVerlet3, 10);
  cout << "---------" << endl;
*/

/*
  // Escape Velocity:
  string Vesc = filename2;
  Vesc.append("Vesc");
  //EscapeVelocity(N, dim, Vesc);


  // Alterative gravitational Force:
  string beta = filename2;
  beta.append("beta");
  AlterativeForce(N, dim);
*/
//  return 0;
//}
}
