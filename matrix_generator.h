#ifndef MATRIX_GENERATOR_H
#define MATRIX_GENERATOR_H

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

class Matrix_Generator
{
   public:
    Matrix_Generator();

    vector<intg> rowM_;
    vector<intg> colM_;
    vector<real> valueM_;

    vector<intg> rowA_;
    vector<intg> colA_;
    vector<real> valueA_;

    vector<vector<double>> Laplacian_Stencil ;
    vector<vector<double>> Helmholtz_Stencil ;

    vector<vector<double>> my_local_stiffness_matrix;
    vector<vector<double>> my_local_mass_matrix;

    vector<map<int,double>> global_stiffness_matrix;
    vector<map<int,double>> global_mass_matrix;

    vector<double> corners;

    Initialisation i;
    void local_stiffness_matrix(face f);

    void print_stiffness_matrix();

     void print_mass_matrix();

    void generate();

    void fillMassVector(void);

    void fillStiffnessVector(void);

    void printVector(void);

 //-------------------------------------- CG ------------------------

  vector<real>& CGprocess(long double &eph);

    vector<real> load_vec;         //load vec
    vector<real> init_vec;
    vector<real> res_vec;
    vector<real> stiff_vec;
    vector<real> dir_vec;
    vector<real> Adk;
    vector<real> xk_new;
    vector<real> rk_new;
     vector<real> dir_new;
     vector<real> u_app;
     vector<real> u_h;


    void setLoad(void);
    void setStiffness(void);
    void residualandDirection(void);
    real normResidual(void);
    real setRo(void);
    real stepSize(real &ro);
    void init_new(real &tk);
    void res_new(real &tk);
    real setBeta(real &ro);
    void setNewdirection(real &beta);
    void CGprint(vector<real> &data);
    void refreshVector(void);

//------------------------------------- inverse power iteration ----------

void inversePower(long double &eph);
vector<real>& U_appNormalize(vector<real>&u_app);
real norm(vector<real> &normVec);
real computeEV(vector<real> &u_h);
void storeFile(vector<real> &u_h,real &ev);

};


#endif // MATRIX_GENERATOR_H

