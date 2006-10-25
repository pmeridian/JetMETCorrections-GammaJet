#include "PluginManager/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "JetMETCorrections/GammaJet/interface/GammaJetProducer.h"
using cms::GammaJet;
DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(GammaJet);
