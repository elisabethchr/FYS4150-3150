/*
Compute the 2D ising model with parallelization
*/
#include <mpi.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>

#include "Metropolis.h"
#include "System.h"

using namespace std;
using namespace arma;


int main(int argc, char* argv[])
{
  string filename;
  int n_spin, MCs, myRank, numprocs;      // number of spins and number of Monte carlo simulations
  double Temp0, Temp_final, Temp_step;

  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &myRank);
  if (argc <= 1){
    cout << "Bad usage: \n";
    cout << "type filename for output, number of cycles, size of lattice, start temperature and final temperature" << endl;
    exit(1);
  }
  else{
    filename = argv[1];
    MCs = atoi(argv[2]);
    n_spin = atoi(argv[3]);
    Temp0 = atof(argv[4]);
    Temp_final = atof(argv[5]);
  }
  //n_spin = 2; //MCs = 100;
  //Temp0 = 1.0; Temp_final = 1.1;
  Temp_step = 0.01;
  double choise = 0; // 0 = orderd, 1= random
  int noIntervalls = MCs/numprocs;
  int myLoopBegin = myRank*noIntervalls + 1;
  int myLoopEnd = (myRank+1)*noIntervalls;
  if ((myRank==numprocs-1) && (myLoopEnd < MCs)) myLoopEnd = MCs;
  cout << myLoopBegin << " " << myLoopEnd<< " " << noIntervalls << endl;

  // Broadcast to all nodes common variable:
  MPI_Bcast (&n_spin, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&Temp0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (&Temp_final, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (&Temp_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  cout << "======================================================" << endl;
  cout << "The Ising model using Number of MC cycles = " << MCs << endl;
  cout << "With parallelization and rank = "<< myRank << endl;
  Metropolis metroplis;
  System sys;


  double TimeStart, TimeEnd, TotTime;
  TimeStart = MPI_Wtime();

  for (double T = Temp0; T<=Temp_final; T+=Temp_step){
    vec ExpValues = zeros(5);
    //vec totExpValues = zeros(5):
    metroplis.metropolisMPI(n_spin, myLoopBegin, myLoopEnd, MCs, T, ExpValues, filename, choise, myRank);

  }
  TimeEnd = MPI_Wtime();
  TotTime = TimeEnd - TimeStart;
  if (myRank==0){
    cout << "Time used: T = " << TotTime << " s on number of processors: "<< numprocs << endl;
  }
  cout << "======================================================" << endl;
  MPI_Finalize();

  return 0;
}
