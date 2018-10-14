#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace arma;


void EnergyTest (double v20, double v2N, double r0, double rN, double eps){
  cout << "Testing energy conservation after one year" << endl;
  double pi = acos(-1);
  double GMsun = 4*pi*pi;
  double mass = 6e24/1.99e30;
  double K_init = 0.5 *mass*v20;
  double U_init = -GMsun*mass/r0;
  double K_final = 0.5 *mass*v2N;
  double U_final = -GMsun*mass/rN;

  double relerrorKin = fabs((K_init - K_final));
  double relerrorPot = fabs((U_init - U_final));

  if (relerrorKin < eps) {
    cout << "   PASSED: Kinetic Energy is conserved, abs(K0-Kfinal) = "<< relerrorKin << endl;
  }
  else{
    cout << "   NOT PASSED: Kinetic Energy is not conserved, abs(K0-Kfinal) = " << relerrorKin << endl;
  }
  if (relerrorPot < eps) {
    cout << "   PASSED: Potetial Energy is conserved, abs(U0-Ufinal) = "<< relerrorPot << endl;
  }
  else{
    cout << "   NOT PASSED: Potential Energy is not conserved, abs(U0-Ufinal) = " << relerrorPot << endl;
  }
}

void PositionTest(double r1, double r2, double eps){
  cout << "Test for same x position after 1 year" << endl;

  if (fabs(r1 - r2) > eps){
    cout << "   NOT PASSED: To big difference in radius after a year " << fabs(r1 - r2)<< endl;
  }
  else{
    cout << "   PASSED: Circular orbit " << fabs(r1 - r2)<< endl;
  }
}


void AngularMomentumTest(mat pos, mat vel, double eps, int N){
  cout << "Test for conservation of angualar momentum, r x p after one year" << endl;

  vec p0 = zeros(3); vec p1 = zeros(3); double diff=0;

  p0[0] = pos(1,0)*vel(2,0)-pos(2,0)*vel(1,0);
  p0[1] = pos(0,0)*vel(2,0)-pos(2,0)*vel(0,0);
  p0[2] = pos(0,0)*vel(1,0)-pos(1,0)*vel(0,0);

  p1[0] = pos(1,N-1)*vel(2,N-1)-pos(2,N-1)*vel(1,N-1);
  p1[1] = pos(0,N-1)*vel(2,N-1)-pos(2,N-1)*vel(0,N-1);
  p1[2] = pos(0,N-1)*vel(1,N-1)-pos(1,N-1)*vel(0,N-1);

  for (int i=0; i<3; i++){
    diff = p1(i)-p0(i);
  }
  if (fabs(diff) < eps){
    cout << "   PASSED: Angular momentum is conserved" << endl;
  }
  else{
    cout << "   NOT PASSED: Angular momentum is not conserved. |p1-p0| = "<< fabs(diff) << " > eps" << endl;
  }


}
