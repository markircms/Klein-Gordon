#ifndef _HEADER_FILE_NAME_H
#define _HEADER_FILE_NAME_H

#include <string>
#include <stdlib.h>
#include <iostream> 
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "TH1F.h"
#include "TH2F.h"
#include "TDirectory.h"
#include "Math/IFunction.h"

#endif

class InitialCondition:public ROOT::Math::IBaseFunctionOneDim{
public:
  double DoEval(double x) const;
  
  ROOT::Math::IBaseFunctionOneDim* Clone() const{
    return new InitialCondition();
  };
  InitialCondition();
  InitialCondition(TF1 f);
private:
  bool a = false;
  TF1 ff = TF1("ff", "0",0,100);
};
