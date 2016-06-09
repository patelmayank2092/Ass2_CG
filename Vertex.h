
/*

Here each vertex has two coordinate x,y.that is intialize by
intialize constructor.

*/



#ifndef VERTEX_H
#define VERTEX_H

#include "type.h"

#include<iostream>
using std::cout;using std::cin;using std::endl;

//typedef double real;

struct vertex
{

    real x=0.0;
    real y=0.0;
    int verNum=0;

    /*vertex(real xh_= 0.0, real yh_=0.0)
    {
        x= xh_;
        y= yh_;
    }*/

};


struct face
{
    size_t v1,v2,v3;
};


#endif //VERTEX_H
