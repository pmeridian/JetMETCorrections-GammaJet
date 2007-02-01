#ifndef GammaJetHelper_h
#define GammaJetHelper_h
#include "JetMETCorrections/Utilities/interface/GetParameters.h"
#include <map>
#include <string>
#include <vector>

class GammaJetHelper{

 public:
  GammaJetHelper(){}
  ~GammaJetHelper(){}
  GammaJetHelper(JetParameters theParam)
  {
    theCalibrationType = getCalibrationType(theParam); 
  };
  
  std::string getCalibrationType(JetParameters);
  std::string getCalibrationType(){return theCalibrationType;}

 private:

  std::string theCalibrationType;
  
};
#endif
