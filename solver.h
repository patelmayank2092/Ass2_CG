
#ifndef SOLVER_H
#define SOLVER_H

#include "matrix_generator.h"
#include "type.h"

#include<iostream>
using std::cout;using std::cin;using std::endl;

#include<vector>
using std::vector;

#include <cmath>

#include<string>
using std::string;

#include<map>
using std::map;

#include<fstream>
using std::ifstream;using std::ofstream;

#include<cmath>
#include<math.h>
#include<stdlib.h>

#include <iomanip>
using  std::setprecision;

#include"Vertex.h"
#include"initialisation.h"

#include"vector.h"

class Solver{


public:
    real numPoint;
    Initialisation i;
    vector<map<intg,real>> mass;
    vector<map<intg,real>> stiff;

    Vector u;
    Vector Au;
    Vector f;
    Vector r;
    Vector d;
    Vector Adk;

    Matrix_Generator mg1;
    Solver(Matrix_Generator &mg,intg &ver,Initialisation &init);
    void doProcess(long double &eph);
    void display(Vector &dis);
    void leftTerm(void);
    void rightTerm(void);
    void computeADK(void);
    void initU(void);
    void inversePower(long double &eph);
    void storeFile(Vector u1,real ev1);
};




#endif
