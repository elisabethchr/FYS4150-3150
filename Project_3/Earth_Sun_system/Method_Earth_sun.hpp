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

mat Initialvelocity(mat vel, int N){

  vel(0,0) = -5.428888690270241E-03;      // AU/day
  vel(1,0) = 1.636353485946535E-02;       // AU/day
  vel(2,0) = -4.491683144318728E-07;      // AU/day

  cout <<"Initial velosity in x, y, z direction:" << endl;
  cout << vel(0,0) << " " << vel(1,0)<<" " << vel(2,0) << endl;
  return vel;
}
mat InitialPosition(mat pos, int N){

  pos(0,0) = 9.528047055398201E-01;       // AU
  pos(1,0) = 3.053612869840809E-01;       // AU
  pos(2,0) = -9.272902073041313E-05;      // AU

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
  double pi = acos(-1);
  double GtimesMsun = 4*pi*pi; //AU^3/yr^2
  double F = -GtimesMsun * pos/((double)r*(double)r*(double)r);

  return F;
}

void Euler(int N, int dim, string filename){
  /*
  Compute the position of the planet using forward Euler method.
  */
  double h = 1.0/((double) N);
  mat vel = zeros(dim, N); mat pos = zeros(dim, N); mat acc = zeros(dim, N);
  vel = Initialvelocity(vel, N);
  pos = InitialPosition(pos, N);

  vec t = zeros(N);

  for (int i=0; i<N-1; i++){
    double r = sqrt(pow(pos(0, i),2) + pow(pos(1, i),2) + pow(pos(2, i),2));

    t(i+1) = t(i) + h;
    for (int j=0; j<dim; j++){

      acc(j,i) = force(pos(j,i), r);
      vel(j, i+1) = vel(j, i) + h*acc(j,i);
      pos(j, i+1) = pos(j, i) + h*vel(j,i);

    }
    cout << acc(0,i) << " " << acc(1,i) << " " << acc(2,i) << endl;
  }
  cout << t(N-1) << " " << h <<" " << N << endl;
  WriteFile(filename, N, dim, pos);
}


void VelocityVerlet(int N, int dim, string filename){
  /*
  Compute the position of the planets using velocity verlet method.
  */

  double h = 1.0/((double) N);
  double hh_half = h*h/2.0;
  double h_half = h/2.0;
  mat acc = zeros(dim, N); mat vel = zeros(dim, N); mat pos = zeros(dim, N);
  vel = Initialvelocity(vel, N);
  pos = InitialPosition(pos, N);

  vec t = zeros(N);
  for (int i=0; i<N-1; i++){
    double r = sqrt(pos(0, i)*pos(0, i) + pos(1, i)*pos(1, i) + pos(2, i)*pos(2, i));

    for (int j=0; j<dim; j++){
      acc(j,i) = force(pos(j,i), r);
      pos(j, i+1) = pos(j, i) + h*vel(j, i) + hh_half*acc(j,i);

    }
    double rnew = sqrt(pow(pos(0,i+1),2) + pow(pos(1,i+1),2) + pow(pos(2,i+1),2));
    for (int j=0; j<dim; j++){
      acc(j,i+1) = force(pos(j,i+1), rnew);
      vel(j, i+1) = vel(j,i) + h_half*(acc(j,i+1) + acc(j,i));

    }
    t(i+1) = t(i) + h;
    //cout << pos(0,i+1) << " " << pos(1,i+1) << " " << r << endl;
  }
  cout << t(N-1) << endl;
  //pos.save("testpos.txt", arma_ascii);
  WriteFile(filename, N, dim, pos);
}
