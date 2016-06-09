
#include <iostream>
using std::cout;using std::cin;using std::endl;

#include <vector>
using std::vector;

#include "type.h"
#include "initialisation.h"
#include "Vertex.h"

#include"matrix_generator.h"



//#include"Colsamm.h"
//using namespace ::_COLSAMM_;


//typedef double real;


int main(int argc, char *argv[]){



    real delta = std::stod(argv[1]);
    long double eph= std::stold(argv[2]);

    Matrix_Generator M;
    M.generate();
    M.inversePower(eph);

 return 0;
}
