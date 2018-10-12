#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <math.h>

using namespace std;
using namespace arma;


int WriteFile(string filename, int N, int dim, mat A){
  ofstream ofile;
  ofile.open(filename);
  //int dim = 3;
  
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  for(int i=0; i<N-1; i++){
    for (int j=0; j < dim ; j++){
      ofile << setw(15) << setprecision(8) << A(j,i) << "  ";
    }
    ofile << "\n";
  }
  ofile.close();
  return 0;
}
