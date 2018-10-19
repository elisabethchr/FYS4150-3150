#include "writetofile.h"

//WriteToFile::WriteToFile(std::string filename, arma::mat A)
//{
int WriteToFile::WritetoFile(std::string filename, arma::mat A){
    int N = sizeof(A.n_cols);
    std::cout << "N = " << N << std::endl;
    int dim = sizeof(A.n_rows);
    std::cout << "dim = " << dim << std::endl;
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

}
