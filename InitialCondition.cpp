#include "InitialCondition.h"


double InitialCondition::DoEval(double x) const {

    if (a) {return ff(x);} else {return 0;}  

};

InitialCondition::InitialCondition(){
    a = false;
  };

InitialCondition::InitialCondition(TF1 f){
  ff = f;
  a = true;
};
