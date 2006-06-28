#ifndef JetCalibratorGammaJet_h
#define JetCalibratorGammaJet_h

#include <map>
#include <string>
#include <vector>
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
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
  reco::CaloJet applyCorrection (const reco::CaloJet& fJet);
  void setParameters(std::string );
   
private:
  
  std::map<double,ParametrizationGammaJet *> parametrization;

  std::string theCalibrationType;
};

#endif
