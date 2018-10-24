#include "writetofile.h"


int WriteToFile::WritetoFileMatrix(std::string filename, arma::mat A, int N){
    int dim = sizeof(A.n_rows);
    FILE* myfile = fopen(filename.c_str(), "w");
    double x; double y;  double z;
    for(int i=0; i<N-1; i++){
        x = A(0,i);
        y = A(1,i);
        z = A(2,i);

        fprintf(myfile, "%.8f     %.8f     %.8f\n", x,y,z);
    }
    fclose(myfile);
    return 0;
}

int WriteToFile::WritetoFile(std::string filename, double x, double y, double z){
    std::ofstream myfile;
    myfile.open(filename.c_str(), std::fstream::app);
    myfile << std::setw(15) << std::setprecision(8) << x << "  " << y << "  " << z << "\n";
    myfile.close();
    return 0;
}

int WriteToFile::WritetoFile_Energy_AngMom(std::string filename, double E, double time){
    std::ofstream myfile;
    myfile.open(filename.c_str(), std::fstream::app);
    myfile << std::setw(15) << std::setprecision(8) << time << "  " << E << "\n";
    myfile.close();
    return 0;

}
