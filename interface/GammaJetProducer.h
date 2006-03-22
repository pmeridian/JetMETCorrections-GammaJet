#ifndef JetProducers_GammaJet_h
#define JetProducers_GammaJet_h

/* Template producer to correct jet
    F.Ratnikov (UMd)
    Mar 2, 2006
*/

#include "FWCore/Framework/interface/EDProducer.h"
#include "DataFormats/Common/interface/EDProduct.h"

#include "JetMETCorrections/GammaJet/interface/JetCalibratorGammaJet.h"

namespace edm {
  class ParameterSet;
  class Event;
  class EventSetup;
}

namespace cms
{
  class GammaJet : public edm::EDProducer
  {
  public:

    // The following is not yet used, but will be the primary
    // constructor when the parameter set system is available.
    //
    explicit GammaJet (const edm::ParameterSet& ps);

    virtual ~GammaJet () {}

    virtual void produce(edm::Event& e, const edm::EventSetup& c);

  private:
    JetCalibratorGammaJet mAlgorithm;
    std::string mInput;
    std::string mTag;
  };
}


#endif
