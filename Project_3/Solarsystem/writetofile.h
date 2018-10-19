#ifndef WRITETOFILE_H
#define WRITETOFILE_H

#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "vec3.h"
#include "solarsystem.h"

class WriteToFile
{
public:
    int WritetoFile(std::string filename, arma::mat A);
};

#endif // WRITETOFILE_H
