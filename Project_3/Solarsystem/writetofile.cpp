#include "writetofile.h"

//WriteToFile::WriteToFile(std::string filename, arma::mat A)
//{
int WriteToFile::WritetoFileMatrix(std::string filename, arma::mat A, int N){

    std::cout << "N = " << N << std::endl;
    int dim = sizeof(A.n_rows);
    std::cout << "dim = " << dim << std::endl;
    FILE* myfile = fopen(filename.c_str(), "w");
    double x; double y;  double z;
    for(int i=0; i<N-1; i++){
        x = A(0,i);
        y = A(1, i);
        z = A(2, i);

        fprintf(myfile, "%.8f     %.8f     %.8f\n", x,y,z);
    }
    fclose(myfile);
    return 0;
}

int WriteToFile::WritetoFile(std::string filename, double x, double y, double z){
    std::ofstream myfile;
    myfile.open(filename.c_str(), std::fstream::app);
    //FILE* myfile = fopen(filename.c_str(), std::fstream::app);//"w");
    myfile << std::setw(15) << std::setprecision(8) << x << "  " << y << "  " << z << "\n";
    //fprintf(myfile, "%.8f     %.8f     %.8f\n", x,y,z);
    //fclose(myfile);
    myfile.close();
    return 0;
}
