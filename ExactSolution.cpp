#include "ExactSolution.h"

double ExactSolution::DoEval(const double *x) const{
    if (a) {return ff(x[0],x[1]);} else return 0;
  };

ExactSolution::ExactSolution(){
    a=false;
  };
ExactSolution::ExactSolution(TF2 f){
    a = true;
    ff = f;
  }; 
