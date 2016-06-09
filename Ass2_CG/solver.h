
#ifndef SOLVER_H
#define SOLVER_H

#include<iostream>
using std::cout;using std::cin;using std::endl;

#include<vector>
using std::vector;

#include "matrix_generator.h"


class Solver{

public:
    Solver();
    Matrix_Generator mg;
    void solverProcess(void);
    void printVector(vector<int> data);



};

#endif
