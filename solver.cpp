
#include "solver.h"

Solver::Solver(){
    ;
}

void Solver::solverProcess(void){


    printVector(mg.rowA_);

    for(auto it=mg.rowA_.begin();it!=mg.rowA_.end();it++){
        cout<<*it<<endl;
    }

    cout<<"here"<<endl;

}
void Solver::printVector(vector<int> data){

    for(auto it=data.begin();it!=data.end();it++){
        cout<<"data is..."<<*it<<endl;
    }


}
