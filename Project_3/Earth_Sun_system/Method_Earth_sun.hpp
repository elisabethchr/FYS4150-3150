// Header file:

#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "WriteToFile.hpp"

using namespace std;
using namespace arma;
//ofstream ofile;     // create file for output

mat Initialvelocity(mat vel){

  vel(0,0) = -3.988919278934853E-03;      // AU/day
  vel(1,0) = 1.674541445029773E-02;       // AU/day
  vel(2,0) = -1.119098958990659E-06;      // AU/day

  cout <<"Initial velosity in x, y, z direction:" << endl;
  cout << vel(0,0) << " " << vel(1,0)<<" " << vel(2,0) << endl;
  return vel;
}
mat InitialPosition(mat pos){

  pos(0,0) = 9.763619062330592E-01;       // AU
  pos(1,0) = 2.225327099640603E-01;       // AU
  pos(2,0) = -8.864996762397152E-05;      // AU
  cout <<"Initial position in x, y, z direction:" << endl;
  cout << pos(0,0) << " " << pos(1,0)<<" " << pos(2,0) << endl;
  return pos;
}

double force (double pos, double r){
  /*
  Compute the Gravitational force, and retun F/M_planet
  */
  double G = 6.67e-11;      // N m^2 kg^2 Gravitational constant
  double Msun = 1.99e+30;   // kg
  double pi = acos(-1); double AU = 1.0; double yr = 1.0;
  double GtimesMsun = 4*pi*pi*Msun; //AU^3/yr^2
  double F = GtimesMsun * pos/(r*r*r);

  return F;
}

void Euler(int N, int dim, mat acc, mat vel, mat pos, string filename){
  /*
  Compute the position of the planet using forward Euler method.
  */
  double h = 1.0/((double) N);
  //mat vel = zeros(dim, N); mat pos = zeros(dim, N);
  vel = Initialvelocity(vel);
  pos = InitialPosition(pos);
  cout << pos(0,0) << endl;

  vec t = zeros(N);
  for (int i; i<N-1; i++){
    double r = sqrt(pow(pos(0, i),2) + pow(pos(1, i),2) + pow(pos(2, i),2));
    
    for (int j=0; j<dim; j++){

      acc(j,i) = -force(pos(j,i), r);
      vel(j, i+1) = vel(j, i) + h*acc(j,i);//force(pos(j,i), r);
      pos(j, i+1) = pos(j, i) + h*vel(j,i);

    }
    t(i+1) = t(i) + h;
    cout << vel(0,i+1) << " " << vel(1,i+1) << " " << r << endl;
  }
  cout << t(N-1) << " "<< h <<" " << N<< endl;
  WriteFile(filename, N, dim, pos);
}


void VelocityVerlet(int N, int dim,mat acc, mat vel, mat pos, string filename){
  /*
  Compute the position of the planets using velocity verlet method.
  */

  double h = 1.0/(N + 1.0);
  double hh_half = h*h/2.0;
  double h_half = h/2.0;

  vel = Initialvelocity(vel);
  pos = InitialPosition(pos);

  vec t = zeros(N);
  for (int i=0; i<N-1; i++){
    double r = sqrt(pos(0, i)*pos(0, i) + pos(1, i)*pos(1, i) + pos(2, i)*pos(2, i));

    for (int j=0; j<dim; j++){
      //acc(j,i) = -force(pos(j,i), r)
      pos(j, i+1) = pos(j, i) + h*vel(j, i) + hh_half*acc(j,i);
      //pos(1, i+1) = pos(1, i) + h*vel(1, i) + hh_half*force(pos(1, i), r);
      //pos(2, i+1) = pos(2, i) + h*vel(2, i) + hh_half*force(pos(2, i), r);
    }
    double rnew = sqrt(pow(pos(0,i+1),2) + pow(pos(1,i+1),2) + pow(pos(2,i+1),2));
    //cout << rnew << " " << i << endl;
    for (int j=0; j<dim; j++){
      acc(j,i+1) = -force(pos(j,i+1), r);
      vel(j, i+1) = vel(j,i) + h_half*(acc(j,i+1) - acc(j,i));
      //vel(1, i+1) = vel(1,i) + h_half*(force(pos(1,i+1), rnew) + force(pos(1,i), r));
      //vel(2, i+1) = vel(2,i) + h_half*(force(pos(2,i+1), rnew) + force(pos(2,i), r));

    }
    t(i+1) = t(i) + h;
    //cout << pos(0,i+1) << " " << pos(1,i+1) << " " << r << endl;
  }
  cout << t(N-1) << endl;
  //pos.save("testpos.txt", arma_ascii);
  WriteFile(filename, N, dim, pos);
}
