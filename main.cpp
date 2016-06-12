
#include <iostream>
using std::cout;using std::cin;using std::endl;

#include <vector>
using std::vector;

#include "type.h"
#include "initialisation.h"
#include "Vertex.h"

#include"matrix_generator.h"

#include "solver.h"

//#include"Colsamm.h"
//using namespace ::_COLSAMM_;


//typedef double real;


int main(int argc, char *argv[]){

    real delta = std::stod(argv[1]);
    long double eph= std::stold(argv[2]);
    string filename=argv[3];

    Initialisation init;

    Matrix_Generator M;
    init=M.generate(filename,delta);

   intg ver=init.numPoint;
   cout<<ver<<endl;

    Solver S(M,ver,init);
    S.inversePower(eph);

    //M.inversePower(eph);

 return 0;
}
