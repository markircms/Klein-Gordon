#ifndef _HEADER_FILE_NAME_H
#define _HEADER_FILE_NAME_H


#include <string>
#include <stdlib.h>
#include <iostream> 
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "TF1.h"
#include "TF2.h"
#include "Math/IFunction.h"


#endif


class ExactSolution:public ROOT::Math::IBaseFunctionMultiDim{
public:
  double DoEval(const double *x) const;/*{
    if (a) {return ff(x[0],x[1]);} else return 0;
  };*/
  
  unsigned int NDim() const
   {
      return 2;
   };
  ROOT::Math::IBaseFunctionMultiDim* Clone() const{
    return new ExactSolution();
  };
  
  ExactSolution();/*{
    a=false;
  };*/
  ExactSolution(TF2 f);/*{
    a = true;
    ff = f;
  };  */
private:
  TF2 ff = TF2("ff", "0",0,100);
  bool a = false;
};
