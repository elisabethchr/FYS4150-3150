//clas for Initialize --> this class is not really necessary, we can include the functions in it within main()
#ifndef INITIALIZE_H
#define INITIALIZE_H
#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <math.h>

using namespace std;
using namespace arma;

class Initialize
{
//private:
public:
    double x0, y0, z0;
    double vx0, vy0, vz0;
//    mat pos, vel, acc;

//public:
   static mat Initialvelocity(mat vel, int dim, int N);
   static mat Initialposition(mat pos, int dim, int N);
};

#endif // INITIALIZE_H
