#include "PluginManager/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "JetMETCorrections/GammaJet/interface/GammaJetProducer.h"
#include "JetMETCorrections/GammaJet/interface/GammaJetAnalysis.h"
using cms::GammaJet;
using cms::GammaJetAnalysis;
DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(GammaJet);
DEFINE_ANOTHER_FWK_MODULE(GammaJetAnalysis);
