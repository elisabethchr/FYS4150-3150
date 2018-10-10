// Header file:

#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <math.h>

using namespace std;
using namespace arma;
ofstream ofile;     // create file for output

void Initialconditions(int N, mat vel, mat pos){//vel, vec vel, vec pos, vec pos, vec pos) {
  //vec vel = zeros(3,N); vec pos = zeros(3,N);
  vel(0,0) = -3.988919278934853E-03;      // AU/day
  vel(1,0) = 1.674541445029773E-02;       // AU/day
  vel(2,0) = -1.119098958990659E-06;      // AU/day
  pos(0,0) = 9.763619062330592E-01;       // AU
  pos(1,0) = 2.225327099640603E-01;       // AU
  pos(2,0) = -8.864996762397152E-05;      // AU
}


double force (double pos, double r){
  /*
  Compute the Gravitational force, and retun F/M_planet
  */
  double G = 6.67e-11;      // N m^2 kg^2 Gravitational constant
  double Msun = 1.99e+30;   // kg
  double F = G*Msun * pos/((double) pow(r, 3));

  return F;
}

void Euler(int N, mat vel, mat pos){
  /*
  Compute the position of the planet using forward Euler method.
  */
  double h = 1.0/(N + 1.0);
  //mat vel = zeros(3,N); mat pos = zeros(3,N);
  Initialconditions(N, vel, pos);//, vel(0,0),vel(1,0),vel(2,0), pos(0,0),pos(1,0),pos(2,0));

  for (int i; i<N; i++){
    double r = sqrt(pos(0, i)*pos(0, i) + pos(1, i)*pos(1, i) + pos(2, i)*pos(2, i));

    double ax = force(pos(0, i), r);
    double ay = force(pos(1, i), r);
    double az = force(pos(2, i), r);

    vel(0, i+1) = vel(0, i) + h*ax;
    vel(1, i+1) = vel(1, i) + h*ay;
    vel(2, i+1) = vel(2, i) + h*az;

    pos(0, i+1) = pos(0, i) + h*vel(0,i);
    pos(1, i+1) = pos(1, i) + h*vel(1,i);
    pos(2, i+1) = pos(2, i) + h*vel(2,i);

  }
}


void VelocityVerlet(int N, mat vel, mat pos){
  /*
  Compute the position of the planets using velocity verlet method.
  */

  double h = 1.0/(N + 1.0);
  double hh_half = h*h/2.0;
  double h_half = h/2.0;
  //mat vel = zeros(3,N); mat pos = zeros(3,N);
  Initialconditions(N, vel, pos);

  for (int i=0; i<N; i++){
    double r = sqrt(pos(0, i)*pos(0, i) + pos(1, i)*pos(1, i) + pos(2, i)*pos(2, i));

    pos(0, i+1) = pos(0, i) + h*vel(0, i) + hh_half*force(pos(0, i), r);
    pos(1, i+1) = pos(1, i) + h*vel(1, i) + hh_half*force(pos(1, i), r);
    pos(2, i+1) = pos(2, i) + h*vel(2, i) + hh_half*force(pos(2, i), r);

    double rnew = sqrt(pow(pos(0, i+1),2) + pow(pos(1, i+1),2) + pow(pos(2, i),2));

    vel(0, i+1) = vel(0,i) + h_half*(force(pos(0,i+1), rnew) + force(pos(0,i), r));
    vel(1, i+1) = vel(1,i) + h_half*(force(pos(1,i+1), rnew) + force(pos(1,i), r));
    vel(2, i+1) = vel(2,i) + h_half*(force(pos(2,i+1), rnew) + force(pos(2,i), r));

  }
}

void WriteFile(string filename, int N, vec vel, vec pos){

  ofile.open(filename);
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  for(int i=0; i<N; i++){
    ofile << setw(15) << setprecision(8) << vel(0,i) << "  ";
    ofile << setw(15) << setprecision(8) << vel(1,i) << "  ";
    ofile << setw(15) << setprecision(8) << vel(2,i) << "  ";
    ofile << setw(15) << setprecision(8) << pos(0,i) << "  ";
    ofile << setw(15) << setprecision(8) << pos(1,i) << "  ";
    ofile << setw(15) << setprecision(8) << pos(2,i) << "  ";
    ofile << "\n";
  }
  ofile.close();
}
