#include "readfile_test.h"
#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>

arma::mat Readfile_test::Readfile_(std::string file)
{
  arma::mat A = arma::zeros(10,3);    // Matrix with the x,y,z elements of each planet in the Solar System
  //A.print();
  std::ifstream infile;
  std::string name; arma::vec vx = arma::zeros(10); arma::vec vy=arma::zeros(10); arma::vec vz = arma::zeros(10);
  infile.open(file);//("Initialposition.txt");

  while (!infile.fail())
  {
    for (int i=0; i<10; i++)
    {
      infile >> name >> vx(i) >> vy(i)>>vz(i);
      //cout << name << " "<<vx(i) <<" "<<vy(i)<<" "<<vz(i) << endl;
      A(i,0) = vx(i);
      A(i,1) = vy(i);
      A(i,2) = vz(i);
      }
    }

    infile.close();
  A.print("output from Readfile_test");
  std::cout << "Readfile_test okay" << std::endl;
  std::cout << "size of matrix n_cols:" << sizeof(A.n_cols);
  return A;
}
