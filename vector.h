
#ifndef VECTOR_H
#define VECTOR_H

#include"type.h"
#include <vector>
using std::vector;

#include<cmath>

#define num 1039

class Vector{

private:


public:
  vector<real> data;
       Vector();
       Vector operator-(Vector &vec);
       real operator*(Vector &vec);
       Vector operator*(real val);
       Vector operator+(Vector vec);
};




#endif
