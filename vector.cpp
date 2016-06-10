
#include "vector.h"


Vector::Vector(){
    data.resize(num,0.0);
}

Vector Vector::operator -(Vector &vec){

    Vector v;
    for(intg i=0;i<num;i++){
        v.data[i] = this->data[i] - vec.data[i];
    }

return v;
}

real Vector::operator*(Vector &vec){

    real fval=0.0;

    for(intg i=0;i<num;i++){
        fval += this->data[i] * vec.data[i];
    }

return fval;
}


Vector Vector::operator*(real val){

    Vector v;

    for(intg i=0;i<num;i++){
      v.data[i]  = this->data[i] * val;
    }

return v;
}

Vector Vector::operator +(Vector vec){

    Vector v;

    for(intg i=0;i<num;i++){
      v.data[i]=this->data[i]  +  vec.data[i];
    }
return v;
}
