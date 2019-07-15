#include "InitialCondition.h"


double InitialCondition::DoEval(double x) const {
    double a = 50;
    double b = 5;
    if (a_tf1) {return ff(x);}
     else {
       if (deriv0 == "phi0" ) // initial condition for phi(t==0)
       {
         return  cos(1.0*(x-a))*exp(-(x-a)*(x-a)/(2*b));
       } else if(deriv0 == "dphi0") //initial conditions for d/dt phi(t==0)
       {
          return  1.*sin(1.0*(x-a))*exp(-(x-a)*(x-a)/(2*b)) + (x-a)/b* cos(1.0*(x-a))*exp(-(x-a)*(x-a)/(2*b));
       } else
       {
         //cout << "try 0 for phi(t=0,x) or 1 for d/dt phi(t=0,x) " << endl;
         return -5555555;}
           }

};

InitialCondition::InitialCondition(){
    a_tf1 = false;
  };

  InitialCondition::InitialCondition(string deriv){
      a_tf1 = false;
      deriv0 = deriv;
    };

InitialCondition::InitialCondition(TF1 f){
  ff = f;
  a_tf1 = true;
};
