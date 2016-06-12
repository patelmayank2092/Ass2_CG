#include"matrix_generator.h"

//#define numPoint 1039


#include"Colsamm.h"
using namespace ::_COLSAMM_;

Matrix_Generator::Matrix_Generator(){
;
}


double k_func(double x, double y)
{
    double delta=0.01;
    real a=(100 + delta);
    real k;
    k=  a*exp(-50*((x*x) + (y*y)))-100;
    return k;
}

void Matrix_Generator::local_stiffness_matrix(face f)
{
    static int face_num =0;
    ELEMENTS::Triangle my_element;
    corners.resize(6,0);

    corners[0] = i.V[f.v1].x; corners[1] = i.V[f.v1].y;
    corners[2] = i.V[f.v2].x; corners[3] = i.V[f.v2].y;
    corners[4] = i.V[f.v3].x; corners[5] = i.V[f.v3].y;


    /*cout<< "Vertices " << endl;
    cout << face_num << '\t'<< f.v1 << '\t' << f.v2 << '\t' << f.v3 << endl;

    cout<< "vertix co-ordinates" <<endl;
    for(int i=0 ; i<5 ;i+=2)
        cout<< corners[i] << '\t' << corners[i+1] << endl;*/

   my_element(corners);


   Laplacian_Stencil = my_element.integrate(grad(v_()) * grad(w_()));

   Helmholtz_Stencil = my_element.integrate(func<double>(k_func) *v_() * w_());

   my_local_stiffness_matrix.resize(3);

   for(auto i = my_local_stiffness_matrix.begin(); i < my_local_stiffness_matrix.end(); ++i)
       i->resize(3);

   for(size_t i = 0; i < 3; ++i)
   {
        for(size_t j = 0; j < 3; ++j)
        {
        my_local_stiffness_matrix[i][j] = Laplacian_Stencil[i][j] - Helmholtz_Stencil[i][j];
        }
   }

   /*cout<< "Local stiffness Matrix" << endl;
   for(auto i = my_local_stiffness_matrix.begin(); i < my_local_stiffness_matrix.end(); ++i)
   {
       for(auto j = i->begin(); j< i->end() ; ++j)
           cout<< *j << '\t';
   cout<<endl;
   }*/

   my_local_mass_matrix = my_element.integrate(v_() * w_());

   face_num++;


   global_stiffness_matrix.resize(numPointMG);

   global_stiffness_matrix[f.v1][f.v1] += my_local_stiffness_matrix[0][0];
   global_stiffness_matrix[f.v2][f.v2] += my_local_stiffness_matrix[1][1];
   global_stiffness_matrix[f.v3][f.v3] += my_local_stiffness_matrix[2][2];

   global_stiffness_matrix[f.v1][f.v2] += my_local_stiffness_matrix[0][1];
   global_stiffness_matrix[f.v2][f.v1] += my_local_stiffness_matrix[1][0];
   global_stiffness_matrix[f.v1][f.v3] += my_local_stiffness_matrix[0][2];
   global_stiffness_matrix[f.v3][f.v1] += my_local_stiffness_matrix[2][0];
   global_stiffness_matrix[f.v2][f.v3] += my_local_stiffness_matrix[1][2];
   global_stiffness_matrix[f.v3][f.v2] += my_local_stiffness_matrix[2][1];

   /*
   cout<< "\nGlobal stiffness matrix" << endl;
   cout<< "row " << f.v1 << " columns " << endl;
   cout << f.v1 << '\t' << '\t' << f.v2 << '\t'<< '\t' << f.v3 << endl;
   cout << global_stiffness_matrix[f.v1][f.v1] << '\t' <<  global_stiffness_matrix[f.v1][f.v2] << '\t' << global_stiffness_matrix[f.v1][f.v3] <<endl;
   */

   global_mass_matrix.resize(numPointMG);

   global_mass_matrix[f.v1][f.v1] +=  my_local_mass_matrix[0][0];
   global_mass_matrix[f.v2][f.v2] +=  my_local_mass_matrix[1][1];
   global_mass_matrix[f.v3][f.v3] +=  my_local_mass_matrix[2][2];

   global_mass_matrix[f.v1][f.v2] +=  my_local_mass_matrix[0][1];
   global_mass_matrix[f.v2][f.v1] +=  my_local_mass_matrix[1][0];
   global_mass_matrix[f.v1][f.v3] +=  my_local_mass_matrix[0][2];
   global_mass_matrix[f.v3][f.v1] +=  my_local_mass_matrix[2][0];
   global_mass_matrix[f.v2][f.v3] +=  my_local_mass_matrix[1][2];
   global_mass_matrix[f.v3][f.v2] +=  my_local_mass_matrix[2][1];

   /*cout<< "\nGlobal mass matrix" << endl;
   cout<< "row " << f.v1 << " columns " << endl;
   cout << f.v1 << '\t' << '\t' << f.v2 << '\t'<< '\t' << f.v3 << endl;
   cout << global_mass_matrix[f.v1][f.v1] << '\t' <<  global_mass_matrix[f.v1][f.v2] << '\t' << global_mass_matrix[f.v1][f.v3] <<endl;
   */
}

void Matrix_Generator::print_stiffness_matrix()
{

      ofstream file_A("A.txt");
      if(file_A.is_open())
      {
      for(size_t i=0;i<numPointMG;i++)
        for(size_t j=0;j<numPointMG;j++)
        {
        if(global_stiffness_matrix[i][j]!=0)
            file_A<<i<<"\t"<< j<<"\t"<<global_stiffness_matrix[i][j]<<"\n";
        }
      }
    else {cout<<"unable to open file"<<endl;}

file_A.close();
}

void Matrix_Generator::print_mass_matrix()
{

      ofstream file_M("M.txt");
      if(file_M.is_open())
      {
      for(size_t i=0;i<numPointMG;i++)
        for(size_t j=0;j<numPointMG;j++)
        {
        if(global_mass_matrix[i][j]!=0)
        file_M<<i<<"\t"<< j<<"\t"<<global_mass_matrix[i][j]<<"\n";
        }
      }
    else {cout<<"unable to open file"<<endl;}

file_M.close();
}


Initialisation& Matrix_Generator::generate(string &filename, real &delta)
{
    i.process(filename,delta);
    //int count =0;
    numPointMG=i.numPoint;
    cout<<numPointMG<<endl;

    for(auto iter=i.F.begin();iter<i.F.end();++iter)
        {

            //if(count<2)
           // {
            local_stiffness_matrix(*iter); //++count;
           // }
           // else
             //   break;
        }
    print_stiffness_matrix();
    print_mass_matrix();

return i;
}





