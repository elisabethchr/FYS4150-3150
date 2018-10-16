#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>

#include "readfile.hpp"
using namespace std;
using namespace arma;

//int main()
mat Readfile(string file)
{
  mat A = zeros(10,3);    // Matrix with the x,y,z elements of each planet in the Solar System
  //A.print();
  ifstream infile;
  string name; vec vx = zeros(10); vec vy=zeros(10); vec vz = zeros(10);
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
  //A.print("poop");
  return A;
}
