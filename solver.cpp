
#include "solver.h"


Solver::Solver(Matrix_Generator &mg,intg &ver,Initialisation &init){

    mg1=mg;
    i=init;
    numPoint=ver;
   // Au.resize(numPoint,0);
   // f.resize(numPoint,0);
   // r.resize(numPoint,0);
   // d.resize(numPoint,0);
    stiff.resize(numPoint);
    mass.resize(numPoint);


 //generate stiff matrix with nonzero
   for(intg i=0;i<numPoint;i++){
    std::map<intg,real>::iterator it_start=mg.global_stiffness_matrix[i].begin();
    std::map<intg,real>::iterator it_end=mg.global_stiffness_matrix[i].end();

    for (auto it=it_start; it!=it_end; ++it){
          if(it->second!=0) stiff[i][it->first]=it->second;
        }
     }

 //generate mass matrix with nonzero
   for(intg i=0;i<numPoint;i++){
    std::map<intg,real>::iterator it_start=mg.global_mass_matrix[i].begin();
    std::map<intg,real>::iterator it_end=mg.global_mass_matrix[i].end();

    for (auto it=it_start; it!=it_end; ++it){
          if(it->second!=0); mass[i][it->first]=it->second;
        }
     }
}

void Solver::display(Vector &dis){

    for (intg it=1; it<5; ++it)
        std::cout << dis.data[it] << "\n";

 std::cout<<std::endl;

 /*

    std::map<intg,real>::iterator it_start=mg.global_stiffness_matrix[0].begin();
    std::map<intg,real>::iterator it_end=mg.global_stiffness_matrix[0].end();

    for (auto it=it_start; it!=it_end; ++it)
        std::cout << it->first << " => " << it->second << '\n';

*/
}


void Solver::doProcess(long double &eph){

    real mul=-1.0,resNorm=0.0,ro=0.0,tk=0.0,beta=0.0;
    intg i=1;


    this->leftTerm();
    this->rightTerm();

    r = Au - f;

    d = r*(-1);

    resNorm = sqrt(r * r);
   // cout<<resNorm<<endl;

     while(resNorm>eph)
    {
       // cout<<"iteration is "<<i<<endl;
        this->computeADK();
        ro = d * Adk;

        tk = (r * d)/(ro*(-1));

        u =u + (d*tk);

        r =r + (Adk*tk);

        beta = (r*Adk)/ro;

        d = (r*(-1)) + (d*beta);

        resNorm = sqrt(r * r);
       // cout<<resNorm<<endl;

        i++;
    }


}


void Solver::inversePower(long double &eph){

    real ev=2.0,ev_old=1.0,norm=0.0;
    int i=1;
    this->initU();
    //cout<<"Outside"<<endl;
    //while(i<3)
       while(fabs((ev-ev_old)/ev_old)>(1e-10))
       {       //cout<<"inside loop"<<endl;
               ev_old = ev;
               //ev = 0.0;

               this->doProcess(eph);

               norm  = sqrt(u*u);
               u = u*(1.0/norm);

               this->leftTerm();
               this->rightTerm();

               ev = (u * Au)/(u *f);

               cout<<"////////////////////////////////////////////"<<endl;
               //cout<<setprecision(10) <<"\n";
               cout<<"    Number of iteration :   "<<setprecision(10)<<i<<"\n";
               cout<<"    New Eigen value :       "<<setprecision(10)<<ev<<"\n";
               cout<<"    Old Eigen value :       "<<setprecision(10)<<ev_old<<"\n";
               //cout<<"    Relative error :        "<<setprecision(10)<<(((ev-ev_old))/ev_old)<<"\n";
               cout<<"    Fabs Relative error :   "<<setprecision(10)<<fabs(((ev-ev_old))/ev_old)<<"\n";
               cout<<"////////////////////////////////////////////\n \n"<<endl;

       i++;
       }

       storeFile(u,ev);

   }


void Solver::leftTerm(void){

  real fval=0.0;

  for(intg i=0;i<numPoint;i++){
    std::map<intg,real>::iterator it_start=stiff[i].begin();
    std::map<intg,real>::iterator it_end=stiff[i].end();

    for (auto it=it_start; it!=it_end; ++it){
          fval += u.data[it->first]*it->second;
    }

      Au.data[i] = fval;
      fval=0.0;
  }
}



void Solver::rightTerm(void){

  real fval=0.0;

    for(intg i=0;i<numPoint;i++){
      std::map<intg,real>::iterator it_start=mass[i].begin();
      std::map<intg,real>::iterator it_end=mass[i].end();

      for (auto it=it_start; it!=it_end; ++it){
            fval += u.data[it->first]*it->second;
      }

        f.data[i] = fval;
        fval=0.0;
    }
}



void Solver::computeADK(void){

  real fval=0.0;

    for(intg i=0;i<numPoint;i++){
      std::map<intg,real>::iterator it_start=stiff[i].begin();
      std::map<intg,real>::iterator it_end=stiff[i].end();

      for (auto it=it_start; it!=it_end; ++it){
            fval += d.data[it->first]*it->second;
      }

        Adk.data[i] = fval;
        fval=0.0;
    }
}




void Solver::initU(void){

    for(int i=0;i<numPoint;i++){
        u.data[i] = u.data[i] + (1.0/sqrt(numPoint));
    }
}





void Solver::storeFile(Vector u1,real ev1){

    ofstream em("eigenmode.txt");
    ofstream eigenv("eigenValue.txt");

    if(em.is_open())
    {
      for(intg j=0;j<numPoint;j++)
      {
          em<<i.V[j].x<<"\t"<<i.V[j].y<<"\t"<<setprecision(10)<<u1.data[j]<<"\n";
            //em<<i.V[j].x<<"\t"<<i.V[j].y<<"\n";
      }
   }
  else {cout<<"unable to open file"<<endl;}

em.close();
///---------------------------------- EV -------------
    if(eigenv.is_open()){

        eigenv<<setprecision(10)<<ev1;
    }
    else throw "unable to open eigenValue.txt";
eigenv.close();

cout<<"storeFile Done"<<endl;

}
