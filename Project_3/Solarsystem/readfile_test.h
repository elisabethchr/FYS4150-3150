#ifndef READFILE_TEST_H
#define READFILE_TEST_H

#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>

class Readfile_test
{
public:
    arma::mat Readfile_(std::string file);
};

#endif // READFILE_TEST_H
