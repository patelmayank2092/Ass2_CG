
#ifndef INITIALISATION_H
#define INITIALISATION_H

#include "type.h"

#include<iostream>
using std::cout;using std::cin;using std::endl;

#include<fstream>
using std::ifstream;using std::ofstream;

#include<string>
using std::string;

#include <vector>
using std::vector;

#include "Vertex.h"

#include <cmath>

//typedef double real;


class Initialisation{

  private:
    vector<real> corners;

  public:
    vector<vertex> V;
    vector<face> F;
    vector<real> K;
    vertex v;
    face f;
    intg numPoint;


    Initialisation();
    void process(string &filename,real &delta);
    void readFile(string &filename);
    void k_waveNumber(real &delta);
};

#endif  // INITIALISATION_H
