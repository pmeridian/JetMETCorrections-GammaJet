#ifndef JetCalibratorGammaJet_h
#define JetCalibratorGammaJet_h

#include <map>
#include <string>
#include <vector>
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/JetObjects/interface/CaloJetCollection.h"
#include "DataFormats/JetObjects/interface/CaloJet.h"
#include "CLHEP/Vector/LorentzVector.h"


class ParametrizationGammaJet;
class CaloJet;

///
/// jet energy corrections from GammaJet calibration
///

class JetCalibratorGammaJet
{
public:  

  JetCalibratorGammaJet() : parametrization(),
                         theCalibrationType() {};
  virtual ~JetCalibratorGammaJet();
  CaloJet applyCorrection (const CaloJet& fJet);
  void setParameters(std::string );
   
private:
  
  std::map<double,ParametrizationGammaJet *> parametrization;

  std::string theCalibrationType;
};

#endif
