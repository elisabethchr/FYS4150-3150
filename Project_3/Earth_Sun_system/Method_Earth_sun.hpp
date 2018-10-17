// Header file - not object oriented code:

#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "WriteToFile.hpp"
//#include "UnitTesting.hpp"

using namespace std;
using namespace arma;
//ofstream ofile;     // create file for output

mat Initialvelocity(mat vel, int N){
  // JPL values are AU/day so multiply with 365.25
  vel(0,0) = -5.428888690270241E-03*365.25;      // AU/yr
  vel(1,0) = 1.636353485946535E-02*365.25;       // AU/yr
  vel(2,0) = -4.491683144318728E-07*365.25;      // AU/yr

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

double step (double stepsPrYear){
  return 1.0/stepsPrYear;
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



void Euler(int N, int dim, string filename, double eps, double dt){
  /*
  Compute the position of the planet using forward Euler method.
  */
  int r2;
  double h = step(365.25)*(dt); //1.0/((double) 365.25);
  mat vel = zeros(dim, N); mat pos = zeros(dim, N); mat acc = zeros(dim, N);
  vel = Initialvelocity(vel, N);
  pos = InitialPosition(pos, N);

  vec t = zeros(N);
  double r0 = sqrt(pow(pos(0,0),2) + pow(pos(1,0),2) + pow(pos(2,0),2));
  double v20 = vel(0,0)*vel(0,0) + vel(1,0)*vel(1,0) + vel(1,0)*vel(1,0);

  clock_t start, finish;
  start = clock();  // start timing

  // Integrating loop:
  for (int i=0; i<N-1; i++){
    double r = sqrt(pow(pos(0, i),2) + pow(pos(1, i),2) + pow(pos(2, i),2));
    double v2 = vel(0,i)*vel(0,i) + vel(1,i)*vel(1,i) + vel(1,i)*vel(1,i);
    t(i+1) = t(i) + h;
    for (int j=0; j<dim; j++){

      acc(j,i) = force(pos(j,i), r);
      vel(j, i+1) = vel(j, i) + h*acc(j,i);
      pos(j, i+1) = pos(j, i) + h*vel(j,i);

    }

    //cout << pos(0,i) << " " << pos(1,i) << " " << pos(2,i) << endl;
  }
  finish = clock();   // end timing
  double time_used = (double)(finish - start)/(CLOCKS_PER_SEC );
  cout << setprecision(10) << "Time used: " << time_used << " s at " << N/365 <<" yr" << endl;

  cout << t(N-1) << " " << h <<" " << N << endl;
  // Final position and velocity
  double rN = sqrt(pow(pos(0, N-1),2) + pow(pos(1, N-1),2) + pow(pos(2, N-1),2));
  double v2N = vel(0,N-1)*vel(0,N-1) + vel(1,N-1)*vel(1,N-1) + vel(1,N-1)*vel(1,N-1);
  cout << "--> UnitTesting -->" << endl;
/*
  PositionTest(pos(0,0), pos(0,366), eps*100);
  EnergyTest(v20, v2N, r0, rN, eps);
  AngularMomentumTest(pos, vel, eps, N);
  WriteFile(filename, N, dim, pos);
*/
}


void VelocityVerlet(int N, int dim, string filename, double eps, double dt){
  /*
  Compute the position of the planets using velocity verlet method.
  */
  double r2;
  double h = step(365.25)* dt;//1.0/(365.25);
  double hh_half = h*h/2.0;
  double h_half = h/2.0;
  mat acc = zeros(dim, N); mat vel = zeros(dim, N); mat pos = zeros(dim, N);
  vel = Initialvelocity(vel, N);
  pos = InitialPosition(pos, N);

  vec t = zeros(N);

  double r0 = sqrt(pow(pos(0,0),2) + pow(pos(1,0),2) + pow(pos(2,0),2));
  double v20 = vel(0,0)*vel(0,0) + vel(1,0)*vel(1,0) + vel(1,0)*vel(1,0);

  clock_t start, finish;
  start = clock();  // start timing

  for (int i=0; i<N-1; i++){
    double r = sqrt(pos(0, i)*pos(0, i) + pos(1, i)*pos(1, i) + pos(2, i)*pos(2, i));
    acc(0,i) = acc(1,i) = acc(2,i) = 0.0;
    for (int j=0; j<dim; j++){
      acc(j,i) = force(pos(j,i), r);
      pos(j, i+1) = pos(j, i) + h*vel(j, i) + hh_half*acc(j,i);


    }
    double rnew = sqrt(pow(pos(0,i+1),2) + pow(pos(1,i+1),2) + pow(pos(2,i+1),2));
    for (int j=0; j<dim; j++){
      acc(j,i+1) = force(pos(j,i+1), rnew);
      vel(j, i+1) = vel(j,i) + h_half*(acc(j,i+1) + acc(j,i));

    }
    //cout << pos(0,i+1) << " " << pos(1,i+1) << " " << pos(2, i+1) << endl;
  }
  finish = clock();   // end timing
  double time_used = (double)(finish - start)/(CLOCKS_PER_SEC );
  cout << setprecision(10) << "Time used: " << time_used << " s at " << N/365 <<" yr" << endl;

  cout << t(N-1) << " " << r2<< endl;

  double rN = sqrt(pow(pos(0, N-1),2) + pow(pos(1, N-1),2) + pow(pos(2, N-1),2));
  double v2N = vel(0,N-1)*vel(0,N-1) + vel(1,N-1)*vel(1,N-1) + vel(1,N-1)*vel(1,N-1);
  cout << "--> UnitTesting -->" << endl;
/*
  PositionTest(pos(0,0), pos(0,366), eps*100);
  EnergyTest(v20, v2N, r0, rN, eps);
  AngularMomentumTest(pos, vel, eps, N);
  //pos.save("testpos.txt", arma_ascii);
  WriteFile(filename, N, dim, pos);
*/
}


double Gforce(double pos, double r, double beta){
  double pi = acos(-1);
  double GtimesMsun = 4*pi*pi; //AU^3/yr^2
  double F = -GtimesMsun * pos/((double) pow(r, beta));

  return F;
}

void EscapeVelocity(int N, int dim, string filename){
  mat acc = zeros(dim, N); mat vel = zeros(dim, N); mat pos = zeros(dim, N);
  pos(0,0) = 1.0; pos(1,0) = 0.0; pos(2,0) = 0.0;
  vel(0,0) = 0.0; vel(2,0) = 0.0;

  double h = step(365.25);//1.0/(365.25);
  double hh_half = h*h/2.0;
  double h_half = h/2.0;
  double V_escape;
  double pi = acos(-1);
  double GMsun = 4*pi*pi;
  double r0 = sqrt(pos(0,0)*pos(0,0) + pos(1,0)*pos(1,0)+ pos(2,0)*pos(2,0));
  V_escape = sqrt(2*GMsun/r0);
  cout << "Analytic escape velocity, V_e = " << V_escape << "AU/yr" << endl;

  for (int vy0=5; vy0 <13; vy0++){
    vel(1,0) = vy0;

    for (int i=0; i<N-1; i++){
      double r = sqrt(pos(0, i)*pos(0, i) + pos(1, i)*pos(1, i) + pos(2, i)*pos(2, i));
      acc(0,i) = acc(1,i) = acc(2,i) = 0.0;
      for (int j=0; j<dim; j++){
        acc(j,i) = force(pos(j,i), r);
        pos(j, i+1) = pos(j, i) + h*vel(j, i) + hh_half*acc(j,i);

      }
      double rnew = sqrt(pow(pos(0,i+1),2) + pow(pos(1,i+1),2) + pow(pos(2,i+1),2));
      for (int j=0; j<dim; j++){
        acc(j,i+1) = force(pos(j,i+1), rnew);
        vel(j, i+1) = vel(j,i) + h_half*(acc(j,i+1) + acc(j,i));

      }
      //cout << pos(0,i+1) << " " << pos(1,i+1) << " " << pos(2, i+1) << endl;
    }
    string filename = "Vescape";
    string arg = to_string(vy0);
    filename.append(arg);
    filename.append(".txt");
    WriteFile(filename, N, dim, pos);
    filename.clear();
  }
}


void AlterativeForce(int N, int dim){
  double h = step(365.25);//1.0/(365.25);
  double hh_half = h*h/2.0;
  double h_half = h/2.0;
  mat acc = zeros(dim, N); mat vel = zeros(dim, N); mat pos = zeros(dim, N);

  vel = Initialvelocity(vel, N);
  pos = InitialPosition(pos, N);

  vec beta = zeros(5);
  beta(0) = 2.0; beta(4) = 3.0; beta(1) = 2.25; beta(2) = 2.5; beta(3) = 2.75;

  for (int b=0; b<5; b++){
    for (int i=0; i<N-1; i++){
      double r = sqrt(pos(0, i)*pos(0, i) + pos(1, i)*pos(1, i) + pos(2, i)*pos(2, i));
      acc(0,i) = acc(1,i) = acc(2,i) = 0.0;
      for (int j=0; j<dim; j++){
        acc(j,i) = Gforce(pos(j,i), r, beta(b));
        pos(j, i+1) = pos(j, i) + h*vel(j, i) + hh_half*acc(j,i);

      }
      double rnew = sqrt(pow(pos(0,i+1),2) + pow(pos(1,i+1),2) + pow(pos(2,i+1),2));
      for (int j=0; j<dim; j++){
        acc(j,i+1) = Gforce(pos(j,i+1), rnew, beta(b));
        vel(j, i+1) = vel(j,i) + h_half*(acc(j,i+1) + acc(j,i));

      }
      //cout << pos(0,i+1) << " " << pos(1,i+1) << " " << pos(2, i+1) << endl;
    }
    string filename = "Beta";
    string arg = to_string(b);
    filename.append(arg);
    filename.append(".txt");
    WriteFile(filename, N, dim, pos);
    filename.clear();
  }



}
