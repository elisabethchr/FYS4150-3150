#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <cstdio>

using namespace std;
using namespace arma;

/*
void writeToFileTest (){
cout << "I'm writing to file! " << endl;
    string filename = "hello.txt";

    ofstream t_file(filename);
    t_file.open(filename);
    t_file << "Hello" << endl;
    t_file.close();

    ofstream taxfile;
    taxfile.open("tax.txt", ios::out | ios::trunc);
    taxfile << "tax" << endl;
    taxfile.close();
*/
/*
    FILE* myfile = fopen(filename.c_str(), "w");
    for(int i=0; i<N-1; i++){
        for (int j=0; j < dim ; j++){
            ofile << setw(15) << setprecision(8) << A(j,i) << "  ";
    fprintf(myfile, "hello\n");
    fclose(myfile);
}
*/
int WriteFile(string filename, int N, int dim, mat A){
    cout << "Hello" << endl;
    FILE* myfile = fopen(filename.c_str(), "w");
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
return 0;


/*
  vec a = linspace<vec>(0, 5, 6);
  ofstream ofile;
  ofile.open("test.txt");
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  for(int i=0; i<N-1; i++){
//  cout << "hello" <<endl;
    ofile << setw(15) << setprecision(8) << a(i) << "  ";
    }
  ofile.close();
*/
/*  ofstream ofile;
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
  */
    //string filename = "hello.txt";
}
