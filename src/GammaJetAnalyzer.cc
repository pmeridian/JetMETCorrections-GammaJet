// -*- C++ -*-
//
// Package:    GammaJetAnalyzer
// Class:      GammaJetAnalyzer
// 
/**\class GammaJetAnalyzer GammaJetAnalyzer.cc MyAnalysis/GammaJetAnalyzer/src/GammaJetAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Daniele del Re
//         Created:  Thu Sep 13 16:00:15 CEST 2007
// $Id: GammaJetAnalyzer.cc,v 1.71 2013/06/12 17:42:53 meridian Exp $
//
//

//
// constructors and destructor
//

#include "JetMETCorrections/GammaJet/interface/GammaJetAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
//#include "DataFormats/VertexReco/interface/VertexCollection.h"
//#include "DataFormats/JetReco/interface/CaloJetCollection.h"
//#include "DataFormats/JetReco/interface/JPTJetCollection.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/JPTJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/ClusterShape.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h" 
#include "RecoCaloTools/Selectors/interface/CaloConeSelector.h"
#include "RecoCaloTools/MetaCollections/interface/CaloRecHitMetaCollections.h"

//#include "RecoEcal/EgammaCoreTools/interface/ClusterShapeAlgo.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

#include "MyAnalysis/IsolationTools/interface/SuperClusterHitsEcalIsolation.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
#include "TrackingTools/GsfTools/interface/GsfPropagatorAdapter.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"

// HLT trigger
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include <FWCore/Common/interface/TriggerNames.h>
#include <DataFormats/Common/interface/TriggerResults.h>
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "PFIsolation/SuperClusterFootprintRemoval/interface/SuperClusterFootprintRemoval.h"

#include "DataFormats/PatCandidates/interface/Photon.h"

#include "TLorentzVector.h"
#include "TRegexp.h"
#include "TString.h"
#include "TVector3.h"

//For HggVertexAnalysis
#include "Analysis/VertexAnalysis/interface/PhotonInfo.h"
#include "JetMETCorrections/GammaJet/interface/GlobalVertexInfo.h"

#include <set>
#include <algorithm>

#include "DataFormats/BTauReco/interface/JetTag.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/PhotonTkIsolation.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

//For BeamHALO
#include "DataFormats/METReco/interface/BeamHaloSummary.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"

//EnergyCorrections
#include "RecoEgamma/EgammaTools/interface/EGEnergyCorrector.h"
//#include "RecoEgamma/EgammaTools/interface/PhotonFix.h"
// #include "HiggsAnalysis/HiggsToGammaGamma/interface/EGEnergyCorrector.h"
// #include "HiggsAnalysis/HiggsToGammaGamma/interface/PhotonFix.h"

#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "RecoEgamma/EgammaTools/interface/ConversionInfo.h"

#include "RecoCaloTools/Navigation/interface/CaloNavigator.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"
#include "RecoCaloTools/Navigation/interface/EcalPreshowerNavigator.h"

#include "DataFormats/Math/interface/deltaR.h"

using namespace edm;
using namespace reco;

//values for Jet Vertex PU ID
#define DEF_GOODVTX_NDOF 4
#define DEF_GOODVTX_Z 24

// 2011 or 2012
#define IS2012 1

// Signed version of delta_phi
inline float GammaJetAnalyzer::delta_phi(float phi1, float phi2) {

  float dphi = fabs(phi1 - phi2);
  float sgn = (phi1 >= phi2 ? +1. : -1.);
  return sgn * (dphi <= TMath::Pi() ? dphi : TMath::TwoPi() - dphi);
}

// Mirror-symmetric delta_eta
inline float GammaJetAnalyzer::delta_eta(float eta1, float eta2) {

  return (eta2 >= 0 ? eta1 - eta2 : eta2 - eta1);
}

inline double GammaJetAnalyzer::oplus(double a, double b) {
  return sqrt(a*a + b*b);
}

// in HF: emf -> R = (Eshort - Elong) / (Ehort + ELong)
// emf = Eshort / (Eshort + Elong)
inline double GammaJetAnalyzer::fixEMF(double emf, double eta) {
  return (fabs(eta)<=3.0 ? emf : 0.5*(1+emf));
}

//used for E2E9 value calculations
inline float GammaJetAnalyzer::recHitE( const  DetId id,  const EcalRecHitCollection &recHits )
{
  if ( id == DetId(0) ) {
    return 0;
  } else {
    EcalRecHitCollection::const_iterator it = recHits.find( id );
    if ( it != recHits.end() ) return (*it).energy();
  }
  return 0;
}


inline float GammaJetAnalyzer::recHitE( const DetId id, const EcalRecHitCollection & recHits,int di, int dj )
{
  // in the barrel:   di = dEta   dj = dPhi
  // in the endcap:   di = dX     dj = dY

  DetId nid;
  if( id.subdetId() == EcalBarrel) nid = EBDetId::offsetBy( id, di, dj );
  else if( id.subdetId() == EcalEndcap) nid = EEDetId::offsetBy( id, di, dj );

  return ( nid == DetId(0) ? 0 : recHitE( nid, recHits ) );
}



inline float GammaJetAnalyzer::recHitApproxEt(  const DetId id,  const EcalRecHitCollection &recHits )
{
  // for the time being works only for the barrel
  if ( id.subdetId() == EcalBarrel ) {
    return recHitE( id, recHits ) / cosh( EBDetId::approxEta( id ) );
  }
  return 0;
}


float GammaJetAnalyzer::GetE2OverE9( const DetId id, const EcalRecHitCollection & recHits)
{ ///////////start calculating e2/e9
  ////http://cmslxr.fnal.gov/lxr/source/RecoLocalCalo/EcalRecAlgos/src/EcalSeverityLevelAlgo.cc#240
  // compute e2overe9
  //   | | | |
  //   +-+-+-+
  //   | |1|2|
  //   +-+-+-+
  //   | | | |
  //   1 - input hit,  2 - highest energy hit in a 3x3 around 1
  //   rechit 1 must have E_t > recHitEtThreshold
  //   rechit 2 must have E_t > recHitEtThreshold2
  //   function returns value of E2/E9 centered around 1 (E2=energy of hits 1+2) if energy of 1>2
  //   if energy of 2>1 and KillSecondHit is set to true, function returns value of E2/E9 centered around 2
  //   *provided* that 1 is the highest energy hit in a 3x3 centered around 2, otherwise, function returns 0
  
float recHitEtThreshold = 10.0;
float recHitEtThreshold2 = 1.0;
bool avoidIeta85=false;
bool KillSecondHit=true;

if ( id.subdetId() == EcalBarrel ) {
  EBDetId ebId( id );
  // avoid recHits at |eta|=85 where one side of the neighbours is missing
  if ( abs(ebId.ieta())==85 && avoidIeta85) return 0;
  // select recHits with Et above recHitEtThreshold
  float e1 = recHitE( id, recHits );
  float ete1=recHitApproxEt( id, recHits );
  // check that rechit E_t is above threshold
  if (ete1 < std::min(recHitEtThreshold,recHitEtThreshold2) ) return 0;
  if (ete1 < recHitEtThreshold && !KillSecondHit ) return 0;

  float e2=-1;
  float ete2=0;
  float s9 = 0;
  // coordinates of 2nd hit relative to central hit
  int e2eta=0;
  int e2phi=0;

  // LOOP OVER 3x3 ARRAY CENTERED AROUND HIT 1
  for ( int deta = -1; deta <= +1; ++deta ) {
    for ( int dphi = -1; dphi <= +1; ++dphi ) {
      // compute 3x3 energy
      float etmp=recHitE( id, recHits, deta, dphi );
      s9 += etmp;
      EBDetId idtmp=EBDetId::offsetBy(id,deta,dphi);
      float eapproxet=recHitApproxEt( idtmp, recHits );
      // remember 2nd highest energy deposit (above threshold) in 3x3 array
      if (etmp>e2 && eapproxet>recHitEtThreshold2 && !(deta==0 && dphi==0)) {
	e2=etmp;
	ete2=eapproxet;
	e2eta=deta;
	e2phi=dphi;
      }
    }
  }

  if ( e1 == 0 )  return 0;
  // return 0 if 2nd hit is below threshold
  if ( e2 == -1 ) return 0;
  // compute e2/e9 centered around 1st hit
  float e2nd=e1+e2;
  float e2e9=0;

  if (s9!=0) e2e9=e2nd/s9;
  // if central hit has higher energy than 2nd hit
  //  return e2/e9 if 1st hit is above E_t threshold
  if (e1 > e2 && ete1>recHitEtThreshold) return e2e9;
  // if second hit has higher energy than 1st hit
  if ( e2 > e1 ) {
    // return 0 if user does not want to flag 2nd hit, or
    // hits are below E_t thresholds - note here we
    // now assume the 2nd hit to be the leading hit.

    if (!KillSecondHit || ete2<recHitEtThreshold || ete1<recHitEtThreshold2) {
      return 0;
    }
    else {
      // LOOP OVER 3x3 ARRAY CENTERED AROUND HIT 2
      float s92nd=0;
      float e2nd_prime=0;
      int e2prime_eta=0;
      int e2prime_phi=0;

      EBDetId secondid=EBDetId::offsetBy(id,e2eta,e2phi);

      for ( int deta = -1; deta <= +1; ++deta ) {
	for ( int dphi = -1; dphi <= +1; ++dphi ) {

	  // compute 3x3 energy
	  float etmp=recHitE( secondid, recHits, deta, dphi );
	  s92nd += etmp;

	  if (etmp>e2nd_prime && !(deta==0 && dphi==0)) {
	    e2nd_prime=etmp;
	    e2prime_eta=deta;
	    e2prime_phi=dphi;
	  }
	}
      }
      // if highest energy hit around E2 is not the same as the input hit, return 0;
      if (!(e2prime_eta==-e2eta && e2prime_phi==-e2phi))
	{
	  return 0;
	}
      // compute E2/E9 around second hit
      float e2e9_2=0;
      if (s92nd!=0) e2e9_2=e2nd/s92nd;
      //   return the value of E2/E9 calculated around 2nd hit
      return e2e9_2;
    }
  }
 } else if ( id.subdetId() == EcalEndcap ) {
  // only used for EB at the moment
  return 0;
 }
return 0;
}

std::vector<float> GammaJetAnalyzer::getESHits(double X, double Y, double Z, std::map<DetId, EcalRecHit> rechits_map, const CaloGeometry& geometry, CaloSubdetectorTopology *topology_p, int row) {
  std::vector<float> esHits;

  const GlobalPoint point(X,Y,Z);

  const CaloSubdetectorGeometry *geometry_p ;
  geometry_p = geometry.getSubdetectorGeometry (DetId::Ecal,EcalPreshower) ;

  DetId esId1 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_p))->getClosestCellInPlane(point, 1);
  DetId esId2 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_p))->getClosestCellInPlane(point, 2);
  ESDetId esDetId1 = (esId1 == DetId(0)) ? ESDetId(0) : ESDetId(esId1);
  ESDetId esDetId2 = (esId2 == DetId(0)) ? ESDetId(0) : ESDetId(esId2);

  std::map<DetId, EcalRecHit>::iterator it;
  ESDetId next;
  ESDetId strip1;
  ESDetId strip2;

  strip1 = esDetId1;
  strip2 = esDetId2;
    
  EcalPreshowerNavigator theESNav1(strip1, topology_p);
  theESNav1.setHome(strip1);
    
  EcalPreshowerNavigator theESNav2(strip2, topology_p);
  theESNav2.setHome(strip2);

  if (row == 1) {
    if (strip1 != ESDetId(0)) strip1 = theESNav1.north();
    if (strip2 != ESDetId(0)) strip2 = theESNav2.east();
  } else if (row == -1) {
    if (strip1 != ESDetId(0)) strip1 = theESNav1.south();
    if (strip2 != ESDetId(0)) strip2 = theESNav2.west();
  }

  // Plane 1 
  if (strip1 == ESDetId(0)) {
    for (int i=0; i<31; ++i) esHits.push_back(0);
  } else {

    it = rechits_map.find(strip1);
    if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
    else esHits.push_back(0);
    //cout<<"center : "<<strip1<<" "<<it->second.energy()<<endl;      

    // east road 
    for (int i=0; i<15; ++i) {
      next = theESNav1.east();
      if (next != ESDetId(0)) {
        it = rechits_map.find(next);
        if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
        else esHits.push_back(0);
        //cout<<"east "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;
      } else {
        for (int j=i; j<15; j++) esHits.push_back(0);
        break;
        //cout<<"east "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }

    // west road 
    theESNav1.setHome(strip1);
    theESNav1.home();
    for (int i=0; i<15; ++i) {
      next = theESNav1.west();
      if (next != ESDetId(0)) {
        it = rechits_map.find(next);
        if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
        else esHits.push_back(0);
        //cout<<"west "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;
      } else {
        for (int j=i; j<15; j++) esHits.push_back(0);
        break;
        //cout<<"west "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }
  }

  if (strip2 == ESDetId(0)) {
    for (int i=0; i<31; ++i) esHits.push_back(0);
  } else {

    it = rechits_map.find(strip2);
    if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
    else esHits.push_back(0);
    //cout<<"center : "<<strip2<<" "<<it->second.energy()<<endl;      

    // north road 
    for (int i=0; i<15; ++i) {
      next = theESNav2.north();
      if (next != ESDetId(0)) {
        it = rechits_map.find(next);
        if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
        else esHits.push_back(0);
        //cout<<"north "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;  
      } else {
        for (int j=i; j<15; j++) esHits.push_back(0);
        break;
        //cout<<"north "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }

    // south road 
    theESNav2.setHome(strip2);
    theESNav2.home();
    for (int i=0; i<15; ++i) {
      next = theESNav2.south();
      if (next != ESDetId(0)) {
        it = rechits_map.find(next);
        if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
        else esHits.push_back(0);
        //cout<<"south "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;
      } else {
        for (int j=i; j<15; j++) esHits.push_back(0);
        break;
        //cout<<"south "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }
  }

  return esHits;
}

std::vector<float> GammaJetAnalyzer::getESShape(std::vector<float> ESHits0, float* esXshape, float* esYshape )
{
  std::vector<float> esShape;
    
  const int nBIN = 21;
  float esRH_F[nBIN];
  float esRH_R[nBIN];
  for (int idx=0; idx<nBIN; idx++) {
    esRH_F[idx] = 0.;
    esRH_R[idx] = 0.;
  }

  for(int ibin=0; ibin<((nBIN+1)/2); ibin++) {
    if (ibin==0) {
      esRH_F[(nBIN-1)/2] = ESHits0[ibin];
      esRH_R[(nBIN-1)/2] = ESHits0[ibin+31];
    } else {
      esRH_F[(nBIN-1)/2+ibin] = ESHits0[ibin];
      esRH_F[(nBIN-1)/2-ibin] = ESHits0[ibin+15];
      esRH_R[(nBIN-1)/2+ibin] = ESHits0[ibin+31];
      esRH_R[(nBIN-1)/2-ibin] = ESHits0[ibin+31+15];
    }
  } 

  // ---- Effective Energy Deposit Width ---- //
  double EffWidthSigmaXX = 0.;
  double EffWidthSigmaYY = 0.;
  double totalEnergyXX   = 0.;
  double totalEnergyYY   = 0.;
  double EffStatsXX      = 0.;
  double EffStatsYY      = 0.;

  std::copy(esRH_F,esRH_F+21,esXshape);
  std::copy(esRH_R,esRH_R+21,esYshape);
  for (int id_X=0; id_X<21; id_X++) {
    totalEnergyXX  += esRH_F[id_X];
    EffStatsXX     += esRH_F[id_X]*(id_X-10)*(id_X-10);
    totalEnergyYY  += esRH_R[id_X];
    EffStatsYY     += esRH_R[id_X]*(id_X-10)*(id_X-10);
  }
  EffWidthSigmaXX  = (totalEnergyXX>0.)  ? sqrt(fabs(EffStatsXX  / totalEnergyXX))   : 0.;
  EffWidthSigmaYY  = (totalEnergyYY>0.)  ? sqrt(fabs(EffStatsYY  / totalEnergyYY))   : 0.;

  esShape.push_back(EffWidthSigmaXX);
  esShape.push_back(EffWidthSigmaYY);
    
  return esShape;
}

bool GammaJetAnalyzer::PhotonMITPreSelection( int photon_index, bool electronVeto) {

  int photon_category = PhotonCategory(photon_index);

  float mitCuts_hoe[4]                 = {0.082,0.075,0.075,0.075};                                        
  float mitCuts_sieie[4]               = {0.014,0.014,0.034,0.034};                                        
  float mitCuts_ecaliso[4]             = {50,4,50,4};                                                      
  float mitCuts_hcaliso[4]             = {50,4,50,4};                                                      
  float mitCuts_trkiso[4]              = {50,4,50,4};                                                      
  //float mitCuts_hcalecal[4]            = {3,3,3,3};                                                        
  //float mitCuts_abstrkiso[4]           = {2.8,2.8,2.8,2.8};                                                
  //float mitCuts_trkiso_hollow03[4]     = {4,4,4,4};                                                       
  //float mitCuts_drtotk_25_99[4]	= {0.26,0.029,0.0062,0.0055};
  float mitCuts_pfiso[4]               = {4,4,4,4};
  
  float val_hoe        = pid_HoverE[photon_index];
  float val_sieie      = pid_etawid[photon_index];                                                          
  float val_ecaliso = pid_jurECAL03[photon_index] - 0.012*ptPhot[photon_index];                              
  float val_hcaliso = pid_twrHCAL03[photon_index] - 0.005*ptPhot[photon_index]; 
  float val_trkiso  = pid_hlwTrack03[photon_index] - 0.002*ptPhot[photon_index]; 
  
  //float val_hcalecal   = (pho_ecalsumetconedr03[photon_index]+pho_hcalsumetconedr03[photon_index]-rho_algo1*rhofac);                                             
  //float val_abstrkiso  = (*pho_tkiso_recvtx_030_002_0000_10_01)[photon_index][vertex_index];                
  //float val_trkiso_hollow03 = pho_trksumpthollowconedr03[photon_index];                                    
  //float val_drtotk_25_99 = pho_drtotk_25_99[photon_index];
  int val_pho_isconv = !hasMatchedPromptElePhot[photon_index];
//   float val_pfiso02 = pid_pfIsoCharged02ForCiC[photon_index][vertex_index];
  
  if (val_hoe             >= mitCuts_hoe[photon_category]         ) return false;                                           
  if (val_sieie           >= mitCuts_sieie[photon_category]       ) return false;
  if (val_ecaliso         >= mitCuts_ecaliso[photon_category]     ) return false;
  if (val_hcaliso         >= mitCuts_hcaliso[photon_category]     ) return false;                                           
  if (val_trkiso          >= mitCuts_trkiso[photon_category]      ) return false;
  //if (val_hcalecal        >= mitCuts_hcalecal[photon_category]    ) return false;
  //if (val_abstrkiso       >= mitCuts_abstrkiso[photon_category]   ) return false;                   
  // if (val_drtotk_25_99    <  mitCuts_drtotk_25_99[photon_category]   ) return false; // Electron Rejection based on CiC for now
  if ((!val_pho_isconv && electronVeto) ) return false; // Electron Rejection based Conversion Safe Veto
  //if (val_trkiso_hollow03 >= mitCuts_trkiso_hollow03[photon_category]) return false;                                        
//   if (val_pfiso02 >= mitCuts_pfiso[photon_category]) return false;            
  
  return true;
}


GammaJetAnalyzer::GammaJetAnalyzer(const edm::ParameterSet& iConfig)
{
  h1_hbherh_detid = new TH1F("hbherh_detid", "", 6, -0.5, 5.5);
    h1_etaPhot = new TH1F("etaPhot", "", 50, -5.5, 5.5);
    h2_n_vs_eta = new TH2D("n_vs_eta", "", 10, 0., 2.5, 1000, 0., 1000.);
  _debug = iConfig.getParameter<bool>("debug");
  outFileName= iConfig.getUntrackedParameter<std::string>("outFileName","output.root");
  // regressionWeights= iConfig.getUntrackedParameter<std::string>("regressionWeights","http://cern.ch/meridian/regweights/gbrv2ph.root");
  regressionWeights= iConfig.getUntrackedParameter<std::string>("regressionWeights","http://cern.ch/meridian/regweights/gbrv3ph_52x.root"); 
  puSummaryInfo_ = iConfig.getParameter<edm::InputTag>("PUSummaryInfoCollection");
  MCTruthCollection_ = iConfig.getUntrackedParameter<edm::InputTag>("MCTruthCollection");
  triggerTag_ = iConfig.getUntrackedParameter<edm::InputTag>("TriggerTag");
  trackTags_ = iConfig.getUntrackedParameter<edm::InputTag>("tracks");
  Vertexsrc_ = iConfig.getUntrackedParameter<edm::InputTag>("vertices");
  Photonsrc_ = iConfig.getUntrackedParameter<edm::InputTag>("Photonsrc");
  Conversionsrc_ = iConfig.getUntrackedParameter<edm::InputTag>("Conversionsrc");
  Electronsrc_ = iConfig.getUntrackedParameter<edm::InputTag>("Electronsrc");
  Jetsrckt4_ = iConfig.getUntrackedParameter<edm::InputTag>("jetskt4");
  Jetsrckt6_ = iConfig.getUntrackedParameter<edm::InputTag>("jetskt6");
  Jetsrcakt5_ = iConfig.getUntrackedParameter<edm::InputTag>("jetsakt5");
  Jetsrcakt7_ = iConfig.getUntrackedParameter<edm::InputTag>("jetsakt7");
  JetJPTsrcak5_ = iConfig.getUntrackedParameter<edm::InputTag>("jetsjptak5");
  JetPFsrckt4_ = iConfig.getUntrackedParameter<edm::InputTag>("jetspfkt4");
  JetPFsrckt6_ = iConfig.getUntrackedParameter<edm::InputTag>("jetspfkt6");
  JetPFsrcakt5_ = iConfig.getUntrackedParameter<edm::InputTag>("jetspfakt5");
  JetPFsrcakt7_ = iConfig.getUntrackedParameter<edm::InputTag>("jetspfakt7");
  JetGensrckt4_ = iConfig.getUntrackedParameter<edm::InputTag>("jetsgenkt4");
  JetGensrckt6_ = iConfig.getUntrackedParameter<edm::InputTag>("jetsgenkt6");
  JetGensrcakt5_ = iConfig.getUntrackedParameter<edm::InputTag>("jetsgenakt5");
  JetGensrcakt7_ = iConfig.getUntrackedParameter<edm::InputTag>("jetsgenakt7");
  METsrc_ = iConfig.getUntrackedParameter<edm::InputTag>("met");
  METGensrc_ = iConfig.getUntrackedParameter<edm::InputTag>("genMet");
  HBhitsrc_ = iConfig.getUntrackedParameter<edm::InputTag>("hbhits");
  recoCollection_    = iConfig.getParameter<std::string>("recoCollection");
  recoProducer_      = iConfig.getParameter<std::string>("recoProducer");
  JetCorrector_akt5_ = iConfig.getParameter<std::string>("JetCorrectionService_akt5");
  JetCorrector_akt7_ = iConfig.getParameter<std::string>("JetCorrectionService_akt7");
  JetCorrector_jptak5_ = iConfig.getParameter<std::string>("JetCorrectionService_jptak5");
  JetCorrector_pfakt5_ = iConfig.getParameter<std::string>("JetCorrectionService_pfakt5");
  JetCorrector_pfakt7_ = iConfig.getParameter<std::string>("JetCorrectionService_pfakt7");
  genjetptthr_ = iConfig.getParameter<double>("genjetptthr");
  calojetptthr_ = iConfig.getParameter<double>("calojetptthr");
  pfjetptthr_ = iConfig.getParameter<double>("pfjetptthr");
  jptjetptthr_ = iConfig.getParameter<double>("jptjetptthr");
  genjetnmin_ = iConfig.getParameter<int>("genjetnmin");
  pfjetnmin_ = iConfig.getParameter<int>("pfjetnmin");
  jptjetnmin_ = iConfig.getParameter<int>("jptjetnmin");
  Xsec_ = iConfig.getParameter<double>("Xsec");
  jetID_ = new reco::helper::JetIDHelper(iConfig.getParameter<ParameterSet>("JetIDParams"));
  dumpBeamHaloInformations_=iConfig.getUntrackedParameter<bool>("dumpBeamHaloInformations",false);
  dumpAKT5Jets_=iConfig.getUntrackedParameter<bool>("dumpAKT5Jets",false);
  dumpAKT7Jets_=iConfig.getUntrackedParameter<bool>("dumpAKT7Jets",false);
  
  dumpJPTAKT5Jets_=iConfig.getUntrackedParameter<bool>("dumpJPTAKT5Jets",false);
  dumpJPTAKT7Jets_=iConfig.getUntrackedParameter<bool>("dumpJPTAKT7Jets",false);
  
  dumpPFAKT5Jets_=iConfig.getUntrackedParameter<bool>("dumpPFAKT5Jets",true);
  dumpPFAKT7Jets_=iConfig.getUntrackedParameter<bool>("dumpPFAKTyJets",true);
  
  dumpKT4Jets_=iConfig.getUntrackedParameter<bool>("dumpKT4Jets",false);
  dumpKT6Jets_=iConfig.getUntrackedParameter<bool>("dumpKT4Jets",false);

  aHLTNames = new std::vector<std::string>;
  aHLTResults = new std::vector<bool>;

  jetMVAAlgos = iConfig.getParameter<std::vector<edm::ParameterSet> >("puJetIDAlgos");
  jetId_algos.resize(jetMVAAlgos.size());

  for(unsigned int imva=0; imva<jetMVAAlgos.size(); imva++){
    jetId_algos[imva] = new PileupJetIdAlgo((jetMVAAlgos.at(imva)));
  }

  //Vertex analysis inizialization
  vtxPar.fixTkIndex=0;

  vtxPar.rescaleTkPtByError=0;
  vtxPar.trackCountThr=0;
  vtxPar.highPurityOnly=0;

  vtxPar.maxD0Signif=999.;
  vtxPar.maxDzSignif=999.;

  vtxPar.removeTracksInCone=1;
  vtxPar.coneSize=0.05;

  vtxPar.useAllConversions=1;

  vtxPar.sigma1Pix=0.016;
  vtxPar.sigma1Tib=0.572;
  vtxPar.sigma1Tob=4.142;
  vtxPar.sigma1PixFwd=0.082;
  vtxPar.sigma1Tid=0.321;
  vtxPar.sigma1Tec=3.109;

  vtxPar.sigma2Pix=0.035;
  vtxPar.sigma2Tib=0.331;
  vtxPar.sigma2Tob=1.564;
  vtxPar.sigma2PixFwd=0.189;
  vtxPar.sigma2Tid=0.418;
  vtxPar.sigma2Tec=0.815;

  vtxPar.singlelegsigma2Pix=0.036;
  vtxPar.singlelegsigma2Tib=0.456;
  vtxPar.singlelegsigma2Tob=0.362;
  vtxPar.singlelegsigma2PixFwd=0.130;
  vtxPar.singlelegsigma2Tid=0.465;
  vtxPar.singlelegsigma2Tec=1.018;

  vtxPar.vtxProbFormula="1.-0.49*(x+1)*(y>0.)";


  bool mvaVertexSelection=1;
  bool addConversionToMva=1;
  std::string tmvaPerVtxWeights="TMVAClassification_BDTCat_conversions_tmva_407.weights.xml";
  tmvaPerVtxMethod="BDTCat";
  //  std::string tmvaPerEvtWeights="TMVAClassification_evtBDTG_conversions_tmva_407.weights.xml";
  std::string tmvaPerEvtWeights="TMVAClassification_BDTvtxprob2012.weights.xml";
  //tmvaPerEvtMethod="evtBDTG";
  tmvaPerEvtMethod="BDTvtxprob2012";

  vtxAna=  new HggVertexAnalyzer(vtxPar);
  vtxAnaFromConv= new HggVertexFromConversions(vtxPar);
  rankVariables.push_back("ptbal"), rankVariables.push_back("ptasym"), rankVariables.push_back("logsumpt2");
  if( tmvaPerVtxWeights != ""  ) {
    tmvaPerVtxVariables_.push_back("ptbal"), tmvaPerVtxVariables_.push_back("ptasym"), tmvaPerVtxVariables_.push_back("logsumpt2");
    if( addConversionToMva ) {
      tmvaPerVtxVariables_.push_back("limPullToConv");
      tmvaPerVtxVariables_.push_back("nConv");
    }
    tmvaPerVtxReader_ = new TMVA::Reader( "!Color:!Silent" );
    HggVertexAnalyzer::bookVariables( *tmvaPerVtxReader_, tmvaPerVtxVariables_ );
    tmvaPerVtxReader_->BookMVA( tmvaPerVtxMethod, tmvaPerVtxWeights );
  } else {
    tmvaPerVtxReader_ = 0;
  }
  if( tmvaPerEvtWeights != "" ) {
    tmvaPerEvtReader_ = new TMVA::Reader( "!Color:!Silent" );
    HggVertexAnalyzer::bookPerEventVariables( *tmvaPerEvtReader_ );
    tmvaPerEvtReader_->BookMVA( tmvaPerEvtMethod, tmvaPerEvtWeights );
  } else {
    tmvaPerEvtReader_ = 0;
  }
  assert( !mvaVertexSelection || tmvaPerVtxReader_ != 0 );

  ecorr_=new EGEnergyCorrector();
  //  PhotonFix::initialiseParameters(iConfig);
  cicPhotonId = new CiCPhotonID(iConfig);


  // chiara --------------------------
  // electron IDs MVAs                                                                                               
  if (IS2012) {
    std::vector<std::string> myManualCatWeigths;
    myManualCatWeigths.push_back(edm::FileInPath("JetMETCorrections/GammaJet/data/Electrons_BDTG_NonTrigV0_Cat1.weights.xml").fullPath());
    myManualCatWeigths.push_back(edm::FileInPath("JetMETCorrections/GammaJet/data/Electrons_BDTG_NonTrigV0_Cat2.weights.xml").fullPath());
    myManualCatWeigths.push_back(edm::FileInPath("JetMETCorrections/GammaJet/data/Electrons_BDTG_NonTrigV0_Cat3.weights.xml").fullPath());
    myManualCatWeigths.push_back(edm::FileInPath("JetMETCorrections/GammaJet/data/Electrons_BDTG_NonTrigV0_Cat4.weights.xml").fullPath());
    myManualCatWeigths.push_back(edm::FileInPath("JetMETCorrections/GammaJet/data/Electrons_BDTG_NonTrigV0_Cat5.weights.xml").fullPath());
    myManualCatWeigths.push_back(edm::FileInPath("JetMETCorrections/GammaJet/data/Electrons_BDTG_NonTrigV0_Cat6.weights.xml").fullPath());
    
    Bool_t manualCat = true;
    myMVANonTrig = new EGammaMvaEleEstimator();
    myMVANonTrig->initialize("BDT",
			     EGammaMvaEleEstimator::kNonTrig,
			     manualCat,
			     myManualCatWeigths);
    
    std::vector<std::string> myManualCatWeigthsTrig;
    myManualCatWeigthsTrig.push_back(edm::FileInPath("JetMETCorrections/GammaJet/data/Electrons_BDTG_TrigV0_Cat1.weights.xml").fullPath());
    myManualCatWeigthsTrig.push_back(edm::FileInPath("JetMETCorrections/GammaJet/data/Electrons_BDTG_TrigV0_Cat2.weights.xml").fullPath());
    myManualCatWeigthsTrig.push_back(edm::FileInPath("JetMETCorrections/GammaJet/data/Electrons_BDTG_TrigV0_Cat3.weights.xml").fullPath());
    myManualCatWeigthsTrig.push_back(edm::FileInPath("JetMETCorrections/GammaJet/data/Electrons_BDTG_TrigV0_Cat4.weights.xml").fullPath());
    myManualCatWeigthsTrig.push_back(edm::FileInPath("JetMETCorrections/GammaJet/data/Electrons_BDTG_TrigV0_Cat5.weights.xml").fullPath());
    myManualCatWeigthsTrig.push_back(edm::FileInPath("JetMETCorrections/GammaJet/data/Electrons_BDTG_TrigV0_Cat6.weights.xml").fullPath());
    myMVATrig = new EGammaMvaEleEstimator();
    myMVATrig->initialize("BDT",
			EGammaMvaEleEstimator::kTrig,
			  manualCat,
			  myManualCatWeigthsTrig);
  }
  // chiara --------------------
}
  
  
GammaJetAnalyzer::~GammaJetAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   //TFile* file_prova = TFile::Open("prova.root", "recreate");
   //file_prova->cd();
   //h1_hbherh_detid->Write();
   //h1_etaPhot->Write();
   //h2_n_vs_eta->Write();
   //file_prova->Close();
  delete cicPhotonId;
  delete ecorr_;

  // chiara
  if (IS2012) {
    delete eIsoFromPFCandsValueMap_;
    delete myMVANonTrig;
    delete myMVATrig;
  }
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
GammaJetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   nMC = nPhot = nElePhot = nJet_akt5 = nJet_kt4 = nJet_kt6 = nJet_pfkt4 = nJet_pfakt5 = nJet_pfkt6 = nJetGen_kt4 = nJetGen_akt5 = nJetGen_kt6 = nEle = 0;
   nJet_jptak5  = 0;
   nJet_pfakt7 = nJet_akt7 = nJetGen_akt7 = 0;
   nvertex = 0;


   //Cleaning double index variables
   for (int iii=0;iii<100;++iii)
     {
       for (int iij=0;iij<100;++iij)
	 {
	   beta_pfakt5[iii][iij]=-1;
	   betaStar_pfakt5[iii][iij]=-1;
	 }
     }

   for (int iii=0;iii<40;++iii)
     {
       for (int iij=0;iij<100;++iij)
	 {
	   pid_hlwTrackForCiC[iii][iij]=-1;
	   pid_hlwTrack03ForCiC[iii][iij]=-1;
	 }
     }

   aHLTNames->clear();
   aHLTResults->clear();

   using reco::TrackCollection;
  
   //   // get generated pt hat
   //   Handle<double> genEventScale; 
   //   iEvent.getByLabel("genEventScale", genEventScale); 

   isMC = !iEvent.isRealData(); // separate MC processing
   store = iEvent.eventAuxiliary().storeNumber(); // study stability across a store
   lbn = iEvent.luminosityBlock(); // sum LBN lumi for normalization
   bx = iEvent.bunchCrossing(); // study effect of out-of-time pile-up
   orbit = iEvent.orbitNumber(); // study beam heating with time (longer bunches)
   run = iEvent.id().run(); // unique ID - part 1
   event = iEvent.id().event(); // unique ID - part 2

   if (!ecorr_->IsInitialized()) ecorr_->Initialize(iSetup,regressionWeights.c_str());  
   //   PhotonFix::initialiseGeometry(iSetup);

   // rho from fast jet
   edm::Handle<double> rhoH;
   //iEvent.getByLabel(edm::InputTag("kt6PFJets","rho","Iso"),rhoH); 
   if( iEvent.getByLabel(edm::InputTag("kt6PFJetsForIso","rho"),rhoH) )
     rho = *rhoH;
   else 
     rho = 0;
   
    edm::Handle<double> rhocaloH;
    //iEvent.getByLabel(edm::InputTag("kt6PFJets","rho","Iso"),rhoH); 
    if( iEvent.getByLabel(edm::InputTag("kt6CaloJetsForIso","rho"),rhocaloH) )
      rhoCalo = *rhocaloH;
    else 
      rhoCalo = 0;

   edm::Handle<double> rhoAllH;
   //iEvent.getByLabel(edm::InputTag("kt6PFJets","rho","Iso"),rhoH); 
   if( iEvent.getByLabel(edm::InputTag("kt6PFJets","rho"),rhoAllH) )
     rhoAllJets = *rhoAllH;
   else 
     rhoAllJets = 0;

   Handle<GenEventInfoProduct> hEventInfo;
   if( isMC ) iEvent.getByLabel("generator", hEventInfo);

   // ------ MC INFORMATION:

   // get MC info from GenParticleCandidates 
   Handle<GenParticleCollection> genParticles;
   if( isMC ) iEvent.getByLabel("genParticles", genParticles);

  // get GEANT sim tracks and vertices (includes conversions)
   Handle<SimTrackContainer> simTracks_h;
   const SimTrackContainer* simTracks;
   if( isMC ) iEvent.getByLabel("g4SimHits", simTracks_h);
   simTracks = (simTracks_h.isValid()) ? simTracks_h.product() : 0;

   Handle<SimVertexContainer> simVert_h;
   const SimVertexContainer* simVertices;
   if( isMC ) iEvent.getByLabel("g4SimHits", simVert_h);
   simVertices = (simVert_h.isValid()) ? simVert_h.product() : 0;

   // ------ RECO INFORMATION:

   // get tracks
   Handle<TrackCollection> tracks;
   iEvent.getByLabel(trackTags_,tracks);
   ConversionFinder convFinder;
   Double_t bfield = 0;
   edm::ESHandle<MagneticField> magneticField;
   iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
   const  MagneticField *mField = magneticField.product();
   bfield = mField->inTesla(GlobalPoint(0.,0.,0.)).z();

   Handle<TrackCollection> tracksHP;
   iEvent.getByLabel("highPurityTracks",tracksHP);

   // get primary vertices
   Handle<VertexCollection> VertexHandle;
   //Handle<vector<Vertex> > VertexHandle;
   iEvent.getByLabel("offlinePrimaryVerticesWithBS", VertexHandle);
   //iEvent.getByLabel(Vertexsrc_, VertexHandle);


   // For jet ID mva -------
   // get primary vertices without BS correction
   Handle<VertexCollection> VertexHandleJetId;
   iEvent.getByLabel("offlinePrimaryVertices", VertexHandleJetId);
   const reco::VertexCollection vertexCollection = *(VertexHandleJetId.product()); 
   const reco::Vertex* selectedVtx  = &(*vertexCollection.begin());;

   // mva computation
   PileupJetIdAlgo jetMVACalculator(*(jetMVAAlgos.begin()));
   
   reco::BeamSpot vertexBeamSpot;
   edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
   iEvent.getByLabel("offlineBeamSpot",recoBeamSpotHandle);
   vertexBeamSpot = *recoBeamSpotHandle;
    
   edm::Handle<reco::ConversionCollection> hConversions;
   iEvent.getByLabel(Conversionsrc_, hConversions);

   // get photons
   Handle<PhotonCollection>  PhotonHandle;
   iEvent.getByLabel(Photonsrc_, PhotonHandle);

   // get electrons
   Handle<GsfElectronCollection>  ElectronHandle;
   iEvent.getByLabel(Electronsrc_, ElectronHandle);

   // get PFCandidates
   Handle<PFCandidateCollection>  PFCandidates;
   iEvent.getByLabel("particleFlow", PFCandidates);

   // get calo jet collection
   Handle<CaloJetCollection> jetsakt5;
   if (dumpAKT5Jets_)
     iEvent.getByLabel(Jetsrcakt5_, jetsakt5);
   Handle<CaloJetCollection> jetsakt7;
   if (dumpAKT7Jets_)
     iEvent.getByLabel(Jetsrcakt7_, jetsakt7);
   Handle<CaloJetCollection> jetskt4;
   if (dumpKT4Jets_)
     iEvent.getByLabel(Jetsrckt4_, jetskt4);
   Handle<CaloJetCollection> jetskt6;
   if (dumpKT6Jets_)
     iEvent.getByLabel(Jetsrckt6_, jetskt6);

   // get JPT collection
   Handle<JPTJetCollection> jptjetsak5;
   if (dumpJPTAKT5Jets_)
     iEvent.getByLabel(JetJPTsrcak5_, jptjetsak5);

   // get PF jets collection
   Handle<PFJetCollection> pfjetskt4;
   if (dumpKT4Jets_)
   iEvent.getByLabel(JetPFsrckt4_, pfjetskt4);
   Handle<PFJetCollection> pfjetsakt5;
   if (dumpPFAKT5Jets_)
     iEvent.getByLabel(JetPFsrcakt5_, pfjetsakt5);
   Handle<PFJetCollection> pfjetsakt7;
   if (dumpPFAKT7Jets_)
     iEvent.getByLabel(JetPFsrcakt7_, pfjetsakt7);
   Handle<PFJetCollection> pfjetskt6;
   if (dumpKT6Jets_)
     iEvent.getByLabel(JetPFsrckt6_, pfjetskt6);

   //get jet correctors
   const JetCorrector* corrector_akt5 = 0;
   const JetCorrector* corrector_jptak5 = 0;
   const JetCorrector* corrector_pfakt5 = 0;

   if (dumpAKT5Jets_)
     corrector_akt5 = JetCorrector::getJetCorrector (JetCorrector_akt5_, iSetup);
   if (dumpJPTAKT5Jets_)
     corrector_jptak5 = JetCorrector::getJetCorrector (JetCorrector_jptak5_, iSetup);
   if (dumpPFAKT5Jets_)
     corrector_pfakt5 = JetCorrector::getJetCorrector (JetCorrector_pfakt5_, iSetup);
 
   // get gen jet collection
   Handle<GenJetCollection> jetsgenkt4;
   if( isMC && dumpKT4Jets_ ) iEvent.getByLabel(JetGensrckt4_, jetsgenkt4);
   Handle<GenJetCollection> jetsgenkt6;
   if( isMC && dumpKT6Jets_ ) iEvent.getByLabel(JetGensrckt6_, jetsgenkt6);
   Handle<GenJetCollection> jetsgenakt5;
   if( isMC && (dumpAKT5Jets_ || dumpJPTAKT5Jets_ || dumpPFAKT5Jets_) ) iEvent.getByLabel(JetGensrcakt5_, jetsgenakt5);
   Handle<GenJetCollection> jetsgenakt7;
   if( isMC && (dumpAKT7Jets_ || dumpJPTAKT7Jets_ || dumpPFAKT7Jets_) ) iEvent.getByLabel(JetGensrcakt7_, jetsgenakt7);

   // get caloMET
   Handle<CaloMETCollection> calomethandle;
   iEvent.getByLabel(METsrc_, calomethandle);

   // get tcMET
   Handle< View<MET> > tcmethandle; 
   iEvent.getByLabel("tcMet", tcmethandle);

   // get pfMET
   Handle<PFMETCollection> pfmethandle;
   iEvent.getByLabel("pfMet", pfmethandle);

   // get gen MET
   Handle<GenMETCollection> genmethandle;
   if( isMC ) iEvent.getByLabel(METGensrc_, genmethandle);

   Handle<GenMETCollection> genmethandle2;
   if( isMC ) iEvent.getByLabel("genMetCalo", genmethandle2);

   // new for 52X
   // Handle<CaloMETCollection>  muJESCorrMEThandle;
   // iEvent.getByLabel("metMuonJESCorAK5", muJESCorrMEThandle);

   // new for 52X
   Handle<PFMETCollection>  pfMetType1Handle;
   // iEvent.getByLabel("metJESCorPFAK5", pfMetType1Handle);
   iEvent.getByLabel("pfType1CorrectedMet", pfMetType1Handle);

   Handle<CaloMETCollection>  caloMEThandleNoHF;
   iEvent.getByLabel("metNoHF", caloMEThandleNoHF);

   Handle<CaloMETCollection>  muCorrMEThandle;
   iEvent.getByLabel("corMetGlobalMuons", muCorrMEThandle);

   // get HCAL info
   Handle<HBHERecHitCollection> hbhe;
   //   iEvent.getByLabel(HBhitsrc_, hbhe);
   iEvent.getByLabel("reducedHcalRecHits","hbhereco",hbhe);

   const HBHERecHitMetaCollection mhbhe(*hbhe);


   // get ECAL reco hits
   Handle<EBRecHitCollection> ecalhitseb;
   const EBRecHitCollection* rhitseb=0;
//    iEvent.getByLabel(recoProducer_, recoCollection_, ecalhitseb);
   iEvent.getByLabel("reducedEcalRecHitsEB", ecalhitseb);
   //const EcalRecHitMetaCollection mecalhits(*ecalhits);    
   rhitseb = ecalhitseb.product(); // get a ptr to the product

   Handle<EERecHitCollection> ecalhitsee;
   const EERecHitCollection* rhitsee=0;
//    iEvent.getByLabel(recoProducer_, "EcalRecHitsEE", ecalhitsee);
   iEvent.getByLabel("reducedEcalRecHitsEE", ecalhitsee);
   rhitsee = ecalhitsee.product(); // get a ptr to the product

   Handle<ESRecHitCollection> ecalhitses;
   const ESRecHitCollection* rhitses=0;
//    iEvent.getByLabel(recoProducer_, "EcalRecHitsES", ecalhitses);
   iEvent.getByLabel("reducedEcalRecHitsES", ecalhitses);
   rhitses = ecalhitses.product(); // get a ptr to the product

   std::map<DetId,EcalRecHit> rechits_map_;
   rechits_map_.clear();

   if (ecalhitses.isValid()) {
     EcalRecHitCollection::const_iterator it;
     for (it = ecalhitses->begin(); it != ecalhitses->end(); ++it) {
       // remove bad ES rechits
       if (it->recoFlag()==1 || it->recoFlag()==14 || (it->recoFlag()<=10 && it->recoFlag()>=5)) continue;
       //Make the map of DetID, EcalRecHit pairs
       rechits_map_.insert(std::make_pair(it->id(), *it));
     }
   }

      edm::Handle<reco::JetTagCollection> combinedSecondaryVertexBJetTags,
                                          combinedSecondaryVertexMVABJetTags,
                                          jetBProbabilityBJetTags,
                                          jetProbabilityBJetTags,
                                          simpleSecondaryVertexHighEffBJetTags,
                                          simpleSecondaryVertexHighPurBJetTags,
                                          softMuonBJetTags,
                                          softMuonByIP3dBJetTags,
                                          softMuonByPtBJetTags,
                                          softElectronBJetTags,
                                          softElectronByIP3dBJetTags,
                                          softElectronByPtBJetTags,
                                          trackCountingHighPurBJetTags,
                                          trackCountingHighEffBJetTags;

      iEvent.getByLabel("combinedSecondaryVertexBJetTags", combinedSecondaryVertexBJetTags);
      iEvent.getByLabel("combinedSecondaryVertexMVABJetTags", combinedSecondaryVertexMVABJetTags);
      iEvent.getByLabel("jetBProbabilityBJetTags", jetBProbabilityBJetTags);
      iEvent.getByLabel("jetProbabilityBJetTags", jetProbabilityBJetTags);
      iEvent.getByLabel("simpleSecondaryVertexHighEffBJetTags", simpleSecondaryVertexHighEffBJetTags);
      iEvent.getByLabel("simpleSecondaryVertexHighPurBJetTags", simpleSecondaryVertexHighPurBJetTags);
      iEvent.getByLabel("softMuonBJetTags", softMuonBJetTags);
      iEvent.getByLabel("softMuonByIP3dBJetTags", softMuonByIP3dBJetTags);
      iEvent.getByLabel("softMuonByPtBJetTags", softMuonByPtBJetTags);
//      iEvent.getByLabel("softElectronBJetTags", softElectronBJetTags);
      iEvent.getByLabel("softElectronByIP3dBJetTags", softElectronByIP3dBJetTags);
      iEvent.getByLabel("softElectronByPtBJetTags", softElectronByPtBJetTags);
      iEvent.getByLabel("trackCountingHighPurBJetTags", trackCountingHighPurBJetTags);
      iEvent.getByLabel("trackCountingHighEffBJetTags", trackCountingHighEffBJetTags);

      edm::Handle<std::vector<PileupSummaryInfo> > pileupHandle;
      if( isMC ) iEvent.getByLabel(puSummaryInfo_, pileupHandle);
      pu_n = 0;
      pu_true_n = 0;

      if( pileupHandle.isValid() ) 
	{
	  std::vector<PileupSummaryInfo>::const_iterator PVI;
	  
	  for(PVI = pileupHandle->begin(); PVI != pileupHandle->end(); ++PVI) 
	    {
	      
	      if(PVI->getBunchCrossing() != 0) 
		continue;
	      
	      pu_n = PVI->getPU_NumInteractions();
	      pu_true_n = PVI->getTrueNumInteractions();
	      //pu_bunchcrossing = PVI->getBunchCrossing();
	      int sv = PVI->getPU_zpositions().size() < 50 ? PVI->getPU_zpositions().size() : 50;
	      for (int iPU=0;iPU<sv;++iPU)
		{
		  pu_zpos[iPU]=PVI->getPU_zpositions()[iPU];
		  pu_sumpt_lowpt[iPU]=PVI->getPU_sumpT_lowpT()[iPU];
		  pu_sumpt_highpt[iPU]=PVI->getPU_sumpT_highpT()[iPU];
		  pu_ntrks_lowpt[iPU]=PVI->getPU_ntrks_lowpT()[iPU];
		  pu_ntrks_highpt[iPU]=PVI->getPU_ntrks_highpT()[iPU];
		}
	    }
	}

  // get geometry
   edm::ESHandle<CaloGeometry> geoHandle;
   //   iSetup.get<IdealGeometryRecord>().get(geoHandle);
   iSetup.get<CaloGeometryRecord>().get(geoHandle);
   const CaloGeometry* geometry = geoHandle.product();

   // get topology
// const CaloSubdetectorTopology *topology_p;
// edm::ESHandle<CaloTopology> topoHandle;
// iSetup.get<CaloTopologyRecord>().get(topoHandle);
// topology_p = topoHandle->getSubdetectorTopology(DetId::Ecal, EcalBarrel);

   edm::ESHandle<CaloTopology> pTopology;
   iSetup.get<CaloTopologyRecord>().get(pTopology);
   const CaloTopology *topology = pTopology.product();

   //   edm::ESHandle<TrackerGeometry> trackerHandle_;
   //edm::ESHandle<MagneticField> theMagField;
   //iSetup.get<IdealMagneticFieldRecord>().get(theMagField);
   //   iSetup.get<TrackerDigiGeometryRecord>().get(trackerHandle_);
  
   //ClusterShapeAlgo algo;

   // chiara -------------------------------------------------------------------
   // for electrons PF isolation and MVA ID
   if (IS2012) {
     eIsoFromPFCandsValueMap_ = new isoContainer(9);
     iEvent.getByLabel( "myElectronPFIsoChHad03",(*eIsoFromPFCandsValueMap_)[0] );
     iEvent.getByLabel( "myElectronPFIsoNHad03",(*eIsoFromPFCandsValueMap_)[1] );
     iEvent.getByLabel( "myElectronPFIsoPhoton03",(*eIsoFromPFCandsValueMap_)[2] );
     iEvent.getByLabel( "myElectronPFIsoChHad04",(*eIsoFromPFCandsValueMap_)[3] );
     iEvent.getByLabel( "myElectronPFIsoNHad04",(*eIsoFromPFCandsValueMap_)[4] );
     iEvent.getByLabel( "myElectronPFIsoPhoton04",(*eIsoFromPFCandsValueMap_)[5] );
     iEvent.getByLabel( "myElectronPFIsoChHad05",(*eIsoFromPFCandsValueMap_)[6] );
     iEvent.getByLabel( "myElectronPFIsoNHad05",(*eIsoFromPFCandsValueMap_)[7] );
     iEvent.getByLabel( "myElectronPFIsoPhoton05",(*eIsoFromPFCandsValueMap_)[8] );
   }

   // transient track builder needed for ele ID MVA
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilder_);  
   const TransientTrackBuilder thebuilder = *(trackBuilder_.product());
   // chiara -------------------------------------------------------------------


//---------------------------HLT Trigger ---------------------------------------------------------------------------------------------
// You Can See HLT Name list ->  " JetMETCorrections/GammaJet/test/HLTList.txt " file 

   hltPass = false;



   hltNamesLen = 0;

   edm::Handle<edm::TriggerResults> hltTriggerResultHandle;
   iEvent.getByLabel(triggerTag_, hltTriggerResultHandle);

    if (!hltTriggerResultHandle.isValid())
    {
      //     std::cout << "invalid handle for HLT TriggerResults" << std::endl;
    }
    else
    {

     edm::TriggerNames HLTNames;
//   HLTNames.init(*hltTriggerResultHandle);
     HLTNames = iEvent.triggerNames(*hltTriggerResultHandle);
     std::string tempnames;
     hltCount = hltTriggerResultHandle->size();
     //std::cout << "hltTriggerResult->size(): " << hltCount << std::endl;
     std::vector<TRegexp> reg;
     reg.push_back(TRegexp(TString(".*Photon.*")));
     reg.push_back(TRegexp(TString(".*Ele.*")));
     reg.push_back(TRegexp(TString(".*MET.*")));
     reg.push_back(TRegexp(TString(".*cosmic.*")));
     reg.push_back(TRegexp(TString(".*Cosmic.*")));
     reg.push_back(TRegexp(TString(".*halo.*")));
     reg.push_back(TRegexp(TString(".*Halo.*")));

     for (int i = 0 ; i != hltCount; ++i) {

        TString hltName_tstr(HLTNames.triggerName(i));
        std::string hltName_str(HLTNames.triggerName(i));
	
	for (unsigned int ireg=0;ireg<reg.size();ireg++)
	  {
	    if ( hltName_tstr.Contains(reg[ireg]) ) 
	      {
		aHLTNames->push_back(hltName_str);
		aHLTResults->push_back(hltTriggerResultHandle->accept(i));
	      }
	  }
     } // for i

     hltNamesLen = tempnames.length();

    } // HLT isValid


    //added Trigger Objects
    edm::InputTag trigEventTag("hltTriggerSummaryAOD","","HLT"); 
    edm::Handle<trigger::TriggerEvent> trigEvent; 
    iEvent.getByLabel(trigEventTag,trigEvent);

    std::vector<std::string> temp_names;
    temp_names.clear();


    //temp_names.push_back("hltEG18CaloId10Iso50TrackIsoDoubleLastFilterUnseeded");
    temp_names.push_back("hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter");
    temp_names.push_back("hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter");
    temp_names.push_back("hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4PMMassFilter");
    temp_names.push_back("hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4TrackIsoFilter");
  
    //new electron triggers:
    temp_names.push_back("hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter");
    temp_names.push_back("hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter");
    temp_names.push_back("hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter");
    temp_names.push_back("hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsoFilter");

    temp_names.push_back("hltEG20CaloIdVLHEFilter");
    temp_names.push_back("hltEG30CaloIdVLHEFilter");
    temp_names.push_back("hltPhoton50CaloIdVLHEFilter");
    temp_names.push_back("hltPhoton75CaloIdVLHEFilter");
    temp_names.push_back("hltPhoton90CaloIdVLHEFilter");

    std::vector<unsigned int>* filter_pass = new std::vector<unsigned int>;
    std::vector<std::string>* filter_names_HLT1 = new std::vector<std::string>;
//     filter_pass->clear();
//     filter_names_HLT1->clear();
    std::vector<std::string>::iterator filter_it;

    //new electron trigger
    std::vector<trigger::TriggerObject> ElectronRefs0;
    std::vector<trigger::TriggerObject> ElectronRefs1;
    std::vector<trigger::TriggerObject> ElectronRefs2;
    std::vector<trigger::TriggerObject> ElectronRefs3;
    std::vector<trigger::TriggerObject> ElectronRefs4;
    std::vector<trigger::TriggerObject> ElectronRefs5;
    std::vector<trigger::TriggerObject> ElectronRefs6;
    std::vector<trigger::TriggerObject> ElectronRefs7;

    std::vector<trigger::TriggerObject> PhotonRefs0;
    std::vector<trigger::TriggerObject> PhotonRefs1;
    std::vector<trigger::TriggerObject> PhotonRefs2;
    std::vector<trigger::TriggerObject> PhotonRefs3;
    std::vector<trigger::TriggerObject> PhotonRefs4;


    //std::vector<trigger::TriggerObject> ElectronRefs00;

    if ( trigEvent.isValid() )
      {
	for (filter_it = temp_names.begin(); filter_it != temp_names.end(); ++filter_it){
	  filter_names_HLT1->push_back((std::string)(*filter_it));
	  const trigger::TriggerObjectCollection & triggerObjects = trigEvent -> getObjects();
	  trigger::size_type filter1_idx = trigEvent -> filterIndex (edm::InputTag(*filter_it,"","HLT") ) ;   
	  trigger::size_type n_filters    = trigEvent -> sizeFilters();
	  if ( filter1_idx < n_filters ) {
	    const trigger::Keys & triggerKeys ( trigEvent -> filterKeys ( filter1_idx ) );
	    const int nkeys = triggerKeys.size();
	    filter_pass->push_back(nkeys);
	    for (int ikey = 0; ikey < nkeys; ++ikey ) {
	      
	    if (*filter_it == "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter") ElectronRefs0.push_back(triggerObjects[ triggerKeys [ikey] ]);
	    else if (*filter_it == "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter") ElectronRefs1.push_back(triggerObjects[ triggerKeys [ikey]]);
	    else if (*filter_it == "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter") ElectronRefs2.push_back(triggerObjects[ triggerKeys [ ikey ] ]);
	    else if (*filter_it == "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter") ElectronRefs3.push_back(triggerObjects[ triggerKeys [ ikey ] ]);
	    else if (*filter_it == "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter") ElectronRefs4.push_back(triggerObjects[ triggerKeys [ ikey ] ]);
	    else if (*filter_it == "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsoFilter") ElectronRefs5.push_back(triggerObjects[ triggerKeys [ ikey ] ]);
	    else if (*filter_it == "hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4PMMassFilter") ElectronRefs6.push_back(triggerObjects[ triggerKeys [ikey] ]);
	    else if (*filter_it == "hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4TrackIsoFilter") ElectronRefs7.push_back(triggerObjects[ triggerKeys [ikey]]);
	    else if (*filter_it == "hltEG20CaloIdVLHEFilter") PhotonRefs0.push_back(triggerObjects[ triggerKeys [ikey]]);
	    else if (*filter_it == "hltEG30CaloIdVLHEFilter") PhotonRefs1.push_back(triggerObjects[ triggerKeys [ikey]]);
	    else if (*filter_it == "hltPhoton50CaloIdVLHEFilter") PhotonRefs2.push_back(triggerObjects[ triggerKeys [ikey]]);
	    else if (*filter_it == "hltPhoton75CaloIdVLHEFilter") PhotonRefs3.push_back(triggerObjects[ triggerKeys [ikey]]);
	    else if (*filter_it == "hltPhoton90CaloIdVLHEFilter") PhotonRefs4.push_back(triggerObjects[ triggerKeys [ikey]]);
	    }
	  }
	  else 
	    filter_pass->push_back(0);
	}
      }


    //final mass filter with L1seeded and SC
    ElectronRefs0_n= 0;
    for (unsigned int i=0; i<ElectronRefs0.size(); i++) {
      if (ElectronRefs0_n >= 8) break;
      trigger::TriggerObject pho=ElectronRefs0[i];
      ElectronRefs0_eta[i] = pho.eta();
      ElectronRefs0_phi[i] = pho.phi();
      ElectronRefs0_et[i] = pho.et();
      ElectronRefs0_n++;
    }

    //L1seeded
    ElectronRefs1_n= 0;
    for (unsigned int i=0; i<ElectronRefs1.size(); i++) {
      if (ElectronRefs1_n >= 8) break;
      trigger::TriggerObject pho=ElectronRefs1[i];
      ElectronRefs1_eta[i] = pho.eta();
      ElectronRefs1_phi[i] = pho.phi();
      ElectronRefs1_et[i] = pho.et();
      ElectronRefs1_n++;
    }


    //Ele17_ele8
    //final mass filter 
    ElectronRefs2_n= 0;
    for (unsigned int i=0; i<ElectronRefs2.size(); i++) {
      if (ElectronRefs2_n >= 8) break;
      trigger::TriggerObject pho=ElectronRefs2[i];
      ElectronRefs2_eta[i] = pho.eta();
      ElectronRefs2_phi[i] = pho.phi();
      ElectronRefs2_et[i] = pho.et();
      ElectronRefs2_n++;
    }

    //L1seeded
    ElectronRefs3_n= 0;
    for (unsigned int i=0; i<ElectronRefs3.size(); i++) {
      if (ElectronRefs3_n >= 8) break;
      trigger::TriggerObject pho=ElectronRefs3[i];
      ElectronRefs3_eta[i] = pho.eta();
      ElectronRefs3_phi[i] = pho.phi();
      ElectronRefs3_et[i] = pho.et();
      ElectronRefs3_n++;
    }
    //final mass filter 
    ElectronRefs4_n= 0;
    for (unsigned int i=0; i<ElectronRefs4.size(); i++) {
      if (ElectronRefs4_n >= 8) break;
      trigger::TriggerObject pho=ElectronRefs4[i];
      ElectronRefs4_eta[i] = pho.eta();
      ElectronRefs4_phi[i] = pho.phi();
      ElectronRefs4_et[i] = pho.et();
      ElectronRefs4_n++;
    }

    //L1seeded
    ElectronRefs5_n= 0;
    for (unsigned int i=0; i<ElectronRefs5.size(); i++) {
      if (ElectronRefs5_n >= 8) break;
      trigger::TriggerObject pho=ElectronRefs5[i];
      ElectronRefs5_eta[i] = pho.eta();
      ElectronRefs5_phi[i] = pho.phi();
      ElectronRefs5_et[i] = pho.et();
      ElectronRefs5_n++;
    }

    //L1seeded
    ElectronRefs6_n= 0;
    for (unsigned int i=0; i<ElectronRefs6.size(); i++) {
      if (ElectronRefs6_n >= 8) break;
      trigger::TriggerObject pho=ElectronRefs6[i];
      ElectronRefs6_eta[i] = pho.eta();
      ElectronRefs6_phi[i] = pho.phi();
      ElectronRefs6_et[i] = pho.et();
      ElectronRefs6_n++;
    }

    //L1seeded
    ElectronRefs7_n= 0;
    for (unsigned int i=0; i<ElectronRefs7.size(); i++) {
      if (ElectronRefs7_n >= 8) break;
      trigger::TriggerObject pho=ElectronRefs7[i];
      ElectronRefs7_eta[i] = pho.eta();
      ElectronRefs7_phi[i] = pho.phi();
      ElectronRefs7_et[i] = pho.et();
      ElectronRefs7_n++;
    }

    //Photon20
    PhotonRefs0_n= 0;
    for (unsigned int i=0; i<PhotonRefs0.size(); i++) {
      if (PhotonRefs0_n >= 8) break;
      trigger::TriggerObject pho=PhotonRefs0[i];
      PhotonRefs0_eta[i] = pho.eta();
      PhotonRefs0_phi[i] = pho.phi();
      PhotonRefs0_et[i] = pho.et();
      PhotonRefs0_n++;
    }

    //Photon30
    PhotonRefs1_n= 0;
    for (unsigned int i=0; i<PhotonRefs1.size(); i++) {
      if (PhotonRefs1_n >= 8) break;
      trigger::TriggerObject pho=PhotonRefs1[i];
      PhotonRefs1_eta[i] = pho.eta();
      PhotonRefs1_phi[i] = pho.phi();
      PhotonRefs1_et[i] = pho.et();
      PhotonRefs1_n++;
    }
    //Photon50
    PhotonRefs2_n= 0;
    for (unsigned int i=0; i<PhotonRefs2.size(); i++) {
      if (PhotonRefs2_n >= 8) break;
      trigger::TriggerObject pho=PhotonRefs2[i];
      PhotonRefs2_eta[i] = pho.eta();
      PhotonRefs2_phi[i] = pho.phi();
      PhotonRefs2_et[i] = pho.et();
      PhotonRefs2_n++;
    }

    //Photon75
    PhotonRefs3_n= 0;
    for (unsigned int i=0; i<PhotonRefs3.size(); i++) {
      if (PhotonRefs3_n >= 8) break;
      trigger::TriggerObject pho=PhotonRefs3[i];
      PhotonRefs3_eta[i] = pho.eta();
      PhotonRefs3_phi[i] = pho.phi();
      PhotonRefs3_et[i] = pho.et();
      PhotonRefs3_n++;
    }
    //Photon90
    PhotonRefs4_n= 0;
    for (unsigned int i=0; i<PhotonRefs4.size(); i++) {
      if (PhotonRefs4_n >= 8) break;
      trigger::TriggerObject pho=PhotonRefs4[i];
      PhotonRefs4_eta[i] = pho.eta();
      PhotonRefs4_phi[i] = pho.phi();
      PhotonRefs4_et[i] = pho.et();
      PhotonRefs4_n++;
    }



    delete filter_pass;
    delete filter_names_HLT1;

    if (!isMC)
      {


	edm::Handle<bool> ecalLaserFilter;
	iEvent.getByLabel("ecalLaserCorrFilter",ecalLaserFilter);
// 	if (ecalLaserFilter.isValid())
 	passEcalLaserFilter =  (*ecalLaserFilter);
	
	edm::Handle<bool> HBHENoiseFilter;
	iEvent.getByLabel("HBHENoiseFilterResultProducer","HBHENoiseFilterResult",HBHENoiseFilter);
	// 	if (HBHENoiseFilter.isValid())
	passHBHENoiseFilter =  (*HBHENoiseFilter);
	
// 	edm::Handle<bool> CSCTightHaloFilter;
// 	iEvent.getByLabel("CSCTightHaloFilter",CSCTightHaloFilter);
// // 	if (CSCTightHaloFilter.isValid())
// 	  passCSCTightHaloFilter =  (*CSCTightHaloFilter);
	
	edm::Handle<bool> hcalLaserEventFilter;
	iEvent.getByLabel("hcalLaserEventFilter",hcalLaserEventFilter);
// 	if (hcalLaserEventFilter.isValid())
	  passhcalLaserEventFilter =  (*hcalLaserEventFilter);

	edm::Handle<bool> EcalDeadCellTriggerPrimitiveFilter;
	iEvent.getByLabel("EcalDeadCellTriggerPrimitiveFilter",EcalDeadCellTriggerPrimitiveFilter);
// 	if (EcalDeadCellTriggerPrimitiveFilter.isValid())
	  passEcalDeadCellTriggerPrimitiveFilter =  (*EcalDeadCellTriggerPrimitiveFilter);

// 	edm::Handle<bool> trackingFailureFilter;
// 	iEvent.getByLabel("trackingFailureFilter",trackingFailureFilter);
// // 	if (trackingFailureFilter.isValid())
// 	  passtrackingFailureFilter =  (*trackingFailureFilter);

	edm::Handle<bool> eeBadScFilter;
	iEvent.getByLabel("eeBadScFilter",eeBadScFilter);
// 	if (eeBadScFilter.isValid())
	  passeeBadScFilter =  (*eeBadScFilter);

      }

//--------------------------------------------------------------------------------------------------------------------------------------

   // Loop over MC truth
   genpt = 0.;
   
   if( isMC ) {

   
     //   genpt = *genEventScale;   
     if (hEventInfo->binningValues().size() > 0)
       genpt = hEventInfo->binningValues()[0];

     genProcessId=hEventInfo->signalProcessID();
     genQScale=hEventInfo->qScale();
     // Figure out the true vertex from the partons:
     // Each event has exactly eight partons (status==3);
     // 0-1: incoming protons (pdgId==2212),
     // 2-3: ISR partons (|pdgId|<=21, udsg only, no c/b),
     // 4-5: incoming partons (|pdgId|<=21, udscg, no b?),
     // 6-7: outgoing partons (|pdgId|<=22, udscg, photons)
     // Every parton 2-7 has the same production vertex (0 for 0-1)
     //assert(genParticles->at(2).status()==3);
     //assert(fabs(genParticles->at(2).pdgId())<100); //<25-><100 to include Z,W
     vxMC = genParticles->at(2).vx();
     vyMC = genParticles->at(2).vy();
     vzMC = genParticles->at(2).vz();

     // Momentum conservation in px (py, pt?): 
     // outgoing + ISR - incoming ~ 0 (couple of GeV)
     // Is the small discrepancy due to the virtual mass of the partons?
     set<int> mothers;
     std::map<const GenParticle*, int> mapMC;

     //Only for debugging purposes
     std::map<const GenParticle*, int> indexMapMC;
     int ipar=0;

     for (GenParticleCollection::const_iterator p = genParticles->begin();
          p != genParticles->end(); ++p) {

       if (_debug)  
	 {
	   // indexMapMC.insert(std::make_pair<const GenParticle*,int>(&(*p),ipar));
	   indexMapMC.insert(std::make_pair(&(*p),ipar));
       	   ++ipar;
	   if (p->pdgId()==22 && p->status()==1 && abs(p->mother()->pdgId())<=21)
	     {
	       std::cout << "FSR PHOTON FOUND " << nMC << "\t" << p->pdgId() << "\t" << p->status() << "\t" << p->p4().pt() << "\t" << p->numberOfDaughters() << "\t" << p->numberOfMothers();
	       if (p->numberOfMothers()==1)
		 std::cout << "\t" << indexMapMC[(const GenParticle*)(&(*p->mother()))] << "\t" << p->mother()->pdgId();
	       std::cout << std::endl;
	     }
	 }

       if (nMC>=(nMaxMC-1)) {continue;}  // to reduce the root file size
       // need space for the mom, eventually

       
       //int motherIDMC = -1;
       if (p->numberOfMothers() > 0) { 
         const Candidate * mom = p->mother();
         for (size_t j = 0; j != genParticles->size(); ++j) {
	   const Candidate * ref = &((*genParticles)[j]);
	   //if (mom->px() == ref->px() && mom->py() == ref->py()
	   //&& mom->pz() == ref->pz() && mom->status() == ref->status()
	   //&& mom->pdgId()==ref->pdgId()) {
	   
           //assert(mom==ref); // address of the candidate is the same?
           //above works in about 99.7% of events
	   
	   if (mom==ref) {
	     motherIDMC[nMC] = j;
           //motherIDMC = j;
	   }
         }
       }

       // Select only a subset of particles to reduce size:
       // All the partons (8)
       // All the stable (status==1) particles within deltaR<1.0
       // around outgoing partons
       // Parents (status==2) of the stable photons and electrons in above
     
       double deltaR = 0;
       if (p->status() == 1) {
         double deltaR1 = oplus(delta_eta(p->eta(),etaMC[6]),
                          delta_phi(p->phi(),phiMC[6]));
         double deltaR2 = oplus(delta_eta(p->eta(),etaMC[7]),
                          delta_phi(p->phi(),phiMC[7]));
         deltaR = min(deltaR1, deltaR2);
       }

       // Neutral particles kept with >200 MeV (ECAL ZS threshold)
       // Charged particles kept with >75 MeV (tracking threshold)
       if (p->status()==3 || 
	   (p->status()==1 && deltaR<0.8 && (p->pt()>0.200 || (p->charge()!=0 && p->pt()>0.075))) || 
	   (p->pdgId()==22 && p->status()==1 && abs(p->mother()->pdgId())<=21 ) ||  //to select also all FSR photons
	   ( (abs(p->pdgId())>=11 && abs(p->pdgId())<=16 ) && p->status()==1 && (abs(p->mother()->pdgId())>=21 && abs(p->mother()->pdgId())<=25) ) //leptons sons of bosons
	   )
	 {
	 
         pdgIdMC[nMC] = p->pdgId();
         statusMC[nMC] = p->status();
         //massMC[nMC] = p->mass();
         ptMC[nMC] = p->pt();
         eMC[nMC] = p->energy();	 
         etaMC[nMC] = p->eta();	 
         phiMC[nMC] = p->phi();	 
         mapMC[&(*p)] = nMC;

	 //Parton and Particle Level Isolation
	 if (p->status()==3 ||
	     (p->pdgId()==22 && p->status()==1 && abs(p->mother()->pdgId())<=21 ) || 
	     ((abs(p->pdgId())>=11 && abs(p->pdgId())<=16 ) && p->status()==1 && (abs(p->mother()->pdgId())>=21 && abs(p->mother()->pdgId())<=25) )
	     )
	   {
	     std::map<std::string,float> isoParticleDR01=particleLevelIsolation(genParticles.product(),&(*p),0.1);
	     std::map<std::string,float> isoParticleDR02=particleLevelIsolation(genParticles.product(),&(*p),0.2);
	     std::map<std::string,float> isoParticleDR03=particleLevelIsolation(genParticles.product(),&(*p),0.3);
	     std::map<std::string,float> isoParticleDR04=particleLevelIsolation(genParticles.product(),&(*p),0.4);
	     std::map<std::string,float> isoParticleDR05=particleLevelIsolation(genParticles.product(),&(*p),0.5);
	     
	     std::map<std::string,float> isoPartonDR01=particleLevelIsolation(genParticles.product(),&(*p),0.1,3);
	     std::map<std::string,float> isoPartonDR02=particleLevelIsolation(genParticles.product(),&(*p),0.2,3);
	     std::map<std::string,float> isoPartonDR03=particleLevelIsolation(genParticles.product(),&(*p),0.3,3);
	     std::map<std::string,float> isoPartonDR04=particleLevelIsolation(genParticles.product(),&(*p),0.4,3);
	     std::map<std::string,float> isoPartonDR05=particleLevelIsolation(genParticles.product(),&(*p),0.5,3);
	     
	     isoParticleChargedDR01MC[nMC]=isoParticleDR01["Charged"];
	     isoParticleChargedDR02MC[nMC]=isoParticleDR02["Charged"];
	     isoParticleChargedDR03MC[nMC]=isoParticleDR03["Charged"];	     
	     isoParticleChargedDR04MC[nMC]=isoParticleDR04["Charged"];
	     isoParticleChargedDR05MC[nMC]=isoParticleDR05["Charged"];	     
	     
	     isoParticleEMNeutralDR01MC[nMC]=isoParticleDR01["EMNeutral"];
	     isoParticleEMNeutralDR02MC[nMC]=isoParticleDR02["EMNeutral"];
	     isoParticleEMNeutralDR03MC[nMC]=isoParticleDR03["EMNeutral"];	     
	     isoParticleEMNeutralDR04MC[nMC]=isoParticleDR04["EMNeutral"];
	     isoParticleEMNeutralDR05MC[nMC]=isoParticleDR05["EMNeutral"];	     
	     
	     isoParticleHADNeutralDR01MC[nMC]=isoParticleDR01["HADNeutral"];
	     isoParticleHADNeutralDR02MC[nMC]=isoParticleDR02["HADNeutral"];
	     isoParticleHADNeutralDR03MC[nMC]=isoParticleDR03["HADNeutral"];	     
	     isoParticleHADNeutralDR04MC[nMC]=isoParticleDR04["HADNeutral"];
	     isoParticleHADNeutralDR05MC[nMC]=isoParticleDR05["HADNeutral"];	     
	     
	     isoPartonDR01MC[nMC]=isoPartonDR01["Charged"]+isoPartonDR01["EMNeutral"]+isoPartonDR01["HADNeutral"];
	     isoPartonDR02MC[nMC]=isoPartonDR02["Charged"]+isoPartonDR02["EMNeutral"]+isoPartonDR02["HADNeutral"];
	     isoPartonDR03MC[nMC]=isoPartonDR03["Charged"]+isoPartonDR03["EMNeutral"]+isoPartonDR03["HADNeutral"];	     
	     isoPartonDR04MC[nMC]=isoPartonDR04["Charged"]+isoPartonDR04["EMNeutral"]+isoPartonDR04["HADNeutral"];
	     isoPartonDR05MC[nMC]=isoPartonDR05["Charged"]+isoPartonDR05["EMNeutral"]+isoPartonDR05["HADNeutral"];	     
	   }
	     
         ++nMC; 
	 
         // if stable photon/electron, find parent
         if (p->status() == 1 && motherIDMC[nMC] != -1
             //&& (pdgIdMC[nMC] == kPhoton || pdgIdMC[nMC] == kElectron)) {//bug
             && (p->pdgId() == kPhoton  || (abs(p->pdgId())>=11 && abs(p->pdgId())>=16)) ) {
	   
           //const Candidate * mom = p->mother();
           const GenParticle *mom = (const GenParticle*)p->mother();
           if (mom->status() == 2
               && (mom->pdgId()<81 || mom->pdgId()>100) // remove MC internal 
               && mothers.find(motherIDMC[nMC]) == mothers.end()) {
	     
	     mothers.insert(motherIDMC[nMC]);
	     
	     if (nMC>=nMaxMC) 
	       {
		 // std::cout << "No possibility to add  the mother" << std::endl;
		 motherIDMC[nMC-1]=-1;
		 continue;
	       }  // to reduce the root file size
	     //       std::cout << "Mother was not present. Fixing the decay tree" << std::endl;
	       motherIDMC[nMC-1]=nMC;
	       pdgIdMC[nMC] = mom->pdgId();
	       statusMC[nMC] = mom->status();
	       //massMC[nMC] = mom->mass();
	       ptMC[nMC] = mom->pt();
	       eMC[nMC] = mom->energy();
	       etaMC[nMC] = mom->eta();
	       phiMC[nMC] = mom->phi(); 
	       
	       mapMC[mom] = nMC;
	       ++nMC; 
	     
	   } // stable photon has parent
	 } // keep particle
	   
	 } // loop particles
     }
     //const double genjetptthr = 5.; // already implicit in GenJet reco
     //const int genjetnmin = 4;
     // Loop over gen Jets
     if (dumpKT4Jets_)
       {
	 for (GenJetCollection::const_iterator it = jetsgenkt4->begin(); 
	      it != jetsgenkt4->end(); ++it) {
	   
	   if (nJetGen_kt4>=100) {cout << "number of gen jets kt 04 is larger than 100. Skipping" << endl; continue;}
	   if (nJetGen_kt4 < genjetnmin_ || it->pt() > genjetptthr_) {
	     
	     ptJetGen_kt4[nJetGen_kt4] = it->pt();
	     eJetGen_kt4[nJetGen_kt4] = it->energy();	 
	     etaJetGen_kt4[nJetGen_kt4] = it->eta();	 
	     phiJetGen_kt4[nJetGen_kt4] = it->phi();	      
	     
	     nJetGen_kt4++;
	   }
	 }
       }

     if (dumpKT6Jets_)
       {
	 for (GenJetCollection::const_iterator it = jetsgenkt6->begin(); 
	      it != jetsgenkt6->end(); ++it) {
	   
	   if (nJetGen_kt6>=100) {cout << "number of gen jets kt 06 is larger than 100. Skipping" << endl; continue;}
	   if (nJetGen_kt6 < genjetnmin_ || it->pt() > genjetptthr_) {
	     
	     ptJetGen_kt6[nJetGen_kt6] = it->pt();
	     eJetGen_kt6[nJetGen_kt6] = it->energy();	 
	     etaJetGen_kt6[nJetGen_kt6] = it->eta();	 
	     phiJetGen_kt6[nJetGen_kt6] = it->phi();	      
	     
	     nJetGen_kt6++;
	   }
	 }
       }

     //----- Figure out the particle decays in tracker volume BEGIN ------

     // Associate grandParents to GenParticles
	 map<const GenParticle*, const SimTrack*> decayedSims;
	 map<const SimTrack*, const GenParticle*> decayedGens;
	 map<const SimTrack*, const SimTrack*> promptParent; // daughter->mother
	 map<const SimTrack*, set<const SimTrack*> > promptDecays; // m->ds
	 map<const SimTrack*, const SimVertex*> promptVertex; // daughter->vertex
     if (_debug)
       cout << Form("Figuring out the particle decays for event %d", event)
        << endl << flush;
     if( simTracks_h.isValid() )
       {
	 // Vertices only return trackID of their parent SimTrack
	 // Figure out the mapping from trackID to SimTrack
	 map<unsigned int, const SimTrack*> trackMap;
	 
	 for (SimTrackContainer::const_iterator iSim = simTracks->begin();
	      iSim != simTracks->end(); ++iSim) {
// 	   cout << iSim->trackId() << endl;
	   if (!iSim->noVertex()) {
// 	     cout << iSim->trackId() << endl;
	     assert(trackMap.find(iSim->trackId())==trackMap.end());
// 	     cout << iSim->trackId() << endl;
	     trackMap[iSim->trackId()] = &(*iSim);
	   }
	 }
	 
	 if (_debug)
	   cout << "Found mapping TrackID->ParentSimTrack" << endl << flush;
	 
	 // Find all SimTracks that come from decays before the ECAL
	 // and find their parent SimTracks
	 
	 for (SimTrackContainer::const_iterator iSim = simTracks->begin();
	      iSim != simTracks->end(); ++iSim) {
	   
	   if (!iSim->noVertex()) {
	     
	     // Find the parent vertex and see if it classifies as an early decay
	     // Exclude the primary vertex (noParent)
	     SimVertex const& vtx = (*simVertices)[iSim->vertIndex()];
	     if (!vtx.noParent() && vtx.position().Rho() < 129 &&
		 fabs(vtx.position().z()) < 304) {
	       
	       // Find parent SimParticle that produced this vertex
	       // vtx->parentIndex is NOT a vector index :( so use trackMap
	       assert(trackMap.find(vtx.parentIndex())!=trackMap.end());
	       const SimTrack* p = trackMap[vtx.parentIndex()];
	       promptParent[&(*iSim)] = p;
	       promptDecays[p].insert(&(*iSim));
	       //promptVertex[p] = &vtx;
	       promptVertex[&(*iSim)] = &vtx;
	     } // early decay
	   } // has vertex
	 } // for simTracks
	 
	 if (_debug)
	   cout << "Found mapping DaughterSimTracks<-(Vertex)->ParentSimTrack"
		<< endl << flush;
	 
	 // Find grandparent SimTracks all the way up the chain
	 map<const SimTrack*, const SimTrack*> chainParents;// g.daughter->grandma
	 map<const SimTrack*, set<const SimTrack*> > chainDecays; // gm->gds
	 
	 for (map<const SimTrack*, const SimTrack*>::const_iterator iSim
		= promptParent.begin(); iSim != promptParent.end(); ++iSim) {
	   
	   // Check that the SimTrack has no daughters itself (=grandchild)
	   if (promptDecays.find(iSim->first)==promptDecays.end()) {
	     // Find the first SimTrack in the parentage chain (=grandparent)
	     const SimTrack *p = iSim->second;
	     while (promptParent.find(p) != promptParent.end())
	       p = promptParent[p];
	     chainParents[iSim->first] = p;
	     chainDecays[p].insert(iSim->first);
	   } // is grandchild
	 } // for promptParent
	 
	 if (_debug)
	   cout << "Found mapping GrandDaughterSimTracks<->GrandParentSimTrack"
		<< " (i.e. identified decay chains)" << endl << flush;
	 
	 // Prune the decay chains to enrich useful information:
	 // - truncate (electron) brems
	 // - truncate non-primary photon conversions
	 
	 for (map<const SimTrack*, set<const SimTrack*> >::const_iterator iSim
		= chainDecays.begin(); iSim != chainDecays.end(); ++iSim) {
	   
	   // iteratively go down the chain and remove decays
	   pruneKids(iSim->first, promptDecays, promptParent, promptVertex, 0);
	 }
	 
	 if (_debug)
	   cout << "Pruned the decay chains" << endl << flush;
	 	 
	 for (map<const SimTrack*, set<const SimTrack*> >::const_iterator iSim
		= chainDecays.begin(); iSim != chainDecays.end(); ++iSim) {
	   
	   if (iSim->first->noGenpart()) {
	     if (_debug)
	       cout << Form("Error: no GenPart found for %d (%1.3g)",
			    iSim->first->type(),
			    iSim->first->momentum().pt()) << endl;
	     //assert(!iSim->first->noGenpart());
	     continue;
	   }
	   
	   // Make sure the decay chain wasn't already pruned out
	   if (promptDecays.find(iSim->first)!=promptDecays.end() &&
	       promptDecays[iSim->first].size()!=0) {
	     
	     // NB: genpartIndex offset by 1
	     const GenParticle* iGen =
	       &(*genParticles)[iSim->first->genpartIndex()-1];
	     assert(iGen->pdgId()==iSim->first->type());
	     decayedSims[iGen] = iSim->first;
	     decayedGens[iSim->first] = iGen;
	   }
	 } // for chainParents 
	 
	 if (_debug)
	   cout << "Found mappings GrandParentSimTracks<->GenParticles" << endl;
	 
	 // Save the particles (conversion) for the primary photon
	 for (map<const GenParticle*, const SimTrack*>::const_iterator iGen
		= decayedSims.begin(); iGen != decayedSims.end(); ++iGen) {
	   
	   const GenParticle *p = iGen->first;
	   if (p->pdgId()==22 && p->mother()->status()==3
	       && p->mother()->pdgId()==22) {
	     
	     if (_debug)
	       cout << "Decay chain for primary photon Gen id 22:" << endl;
	     bool saved = printChildren(decayedSims[p], promptDecays,
					promptVertex, 0, true);
	     if (saved && mapMC.find(p)!=mapMC.end()) {
	       statusMC[mapMC[p]] *= -1;
	     }
	   } // is primary photon
	 } // for iGen 
       }
     
     if (_debug)
       cout << "Found mappings for primary photon (if any)" << endl;
     
     //----- Figure out the particle decays in tracker volume END ------

     if( dumpAKT5Jets_ || dumpJPTAKT5Jets_ || dumpPFAKT5Jets_ )
       {
	 for (GenJetCollection::const_iterator it = jetsgenakt5->begin(); 
	      it != jetsgenakt5->end(); ++it) {
	   
	   if (nJetGen_akt5>=100) {cout << "number of gen jets kt 05 is larger than 100. Skipping" << endl; continue;}
	   if (nJetGen_akt5 < genjetnmin_ || it->pt() > genjetptthr_) {
	     
	     ptJetGen_akt5[nJetGen_akt5] = it->pt();	 
	     eJetGen_akt5[nJetGen_akt5] = it->energy();	 
	     etaJetGen_akt5[nJetGen_akt5] = it->eta();	 
	     phiJetGen_akt5[nJetGen_akt5] = it->phi();    
	     
	     if( simTracks_h.isValid() )
	       {
		 // Extra variables for PFlow
		 Int_t nMuonsGen = 0;
		 TLorentzVector p4MuonsGen;
		 
		 Int_t nElectronsGen = 0;
		 TLorentzVector p4ElectronsGen;
		 
		 Int_t nPhotonsGen = 0;
		 TLorentzVector p4PhotonsGen;
		 
		 Int_t nTracksGen = 0;
		 TLorentzVector p4TracksGen;
		 
		 Int_t nNeutralHadronsGen = 0;
		 TLorentzVector p4NeutralHadronsGen;
		 
		 Int_t nHFHadronsGen = 0;
		 TLorentzVector p4HFHadronsGen;
		 
		 Int_t nHFEMGen = 0;
		 TLorentzVector p4HFEMGen;
	     
		 Int_t nNeutronsGen = 0;
		 TLorentzVector p4NeutronsGen;
		 
		 Int_t nK0LGen = 0;
		 TLorentzVector p4K0LGen;
		 
		 Int_t nK0SGen = 0;
		 TLorentzVector p4K0SGen;
		 
		 Int_t nLambdasGen = 0;
		 TLorentzVector p4LambdasGen;
		 
		 Int_t nCsiGen = 0;
		 TLorentzVector p4CsiGen;
		 
		 Int_t nOtherNeutralHadronsGen = 0;
		 TLorentzVector p4OtherNeutralHadronsGen;
		 
		 
		 vector<const GenParticle*> jetParticles = it->getGenConstituents();
		 vector<const GenParticle*> shortPtcls;
		 
		 //----- Select the particle decays in tracker volume BEGIN ------      
		 
		 const int ijet = it - jetsgenakt5->begin();
		 if( simTracks_h.isValid() ){
		   if (_debug)
		     cout << Form("Found mappings GenParticles<->GrandParentSimTracks"
				  " for jet %d of pT=%1.3f eta=%1.3g",
				  ijet, it->pt(), it->eta())
			  << endl << flush;
		   
		   // Print out decay chains for this jet iteratively
		   // NB: this method will also save the SIM particles when asked
		   for (vector< const GenParticle* >::const_iterator iGen
			  = jetParticles.begin(); iGen != jetParticles.end(); ++iGen) {
		     
		     if (_debug)
		       cout << Form("Gen id %d (pT=%1.3g GeV, eta=%1.3g)",
				    (*iGen)->pdgId(), (*iGen)->pt(), (*iGen)->eta());
		     
		     if (decayedSims.find(*iGen) != decayedSims.end()) {
		       
		       if (_debug)
			 cout << " decay chain:" << endl;
		       assert((*iGen)->pdgId()==decayedSims[*iGen]->type());
		       bool saved = printChildren(decayedSims[*iGen],
						  promptDecays, promptVertex, 0,
						  (ijet<2 ? true : false)); // save to file?
		       if (saved && mapMC.find(*iGen)!=mapMC.end()) {
			 //assert(statusMC[mapMC[*iGen]]!=2);
			 //assert(mapMC.find(*iGen)!=mapMC.end());
			 statusMC[mapMC[*iGen]] *= -1;
		       }
		 }
		     else
		       if (_debug)
			 cout << endl;
		   }
		 }
		 //----- Select out the particle decays in tracker volume END------
		 
		 // Sum up the different types of GenParticle energy
		 for (vector< const GenParticle* >::iterator iPart = jetParticles.begin();
		      iPart != jetParticles.end(); ++iPart) {
		   
		   Int_t partPdgId = (*iPart)->pdgId();
		   // Convert particle momentum to normal TLorentzVector, wrong type :(
		   math::XYZTLorentzVectorD const& p4t = (*iPart)->p4();
		   TLorentzVector p4(p4t.px(), p4t.py(), p4t.pz(), p4t.energy());
		   
		   if (fabs(partPdgId)==13) { //muons
		     nMuonsGen += 1;
		     p4MuonsGen += p4;
		   } else if (fabs(partPdgId)==11) { //electrons
		     if( fabs((*iPart)->eta())<2.5) { //tracker acceptance
		       nElectronsGen += 1;
		       p4ElectronsGen += p4;
		     } else if( fabs((*iPart)->eta())<3. ) { //if no track, they're just photons
		       nPhotonsGen += 1;
		       p4PhotonsGen += p4;
		     } else { //HFEM in HF 
		       nHFEMGen += 1;
		       p4HFEMGen += p4;
		     }
		   } else if ((*iPart)->charge() != 0) { // charged hadrons
		     if( fabs((*iPart)->eta())<2.5) { //tracker acceptance
		       nTracksGen += 1;
		       p4TracksGen += p4;
		     } else { //if no track, they're just hadrons
		       if( fabs((*iPart)->eta())<3. ) { //HCAL acceptance
			 nNeutralHadronsGen += 1;
			 p4NeutralHadronsGen += p4;
		       } else { //HF
			 nHFHadronsGen += 1;
			 p4HFHadronsGen += p4;
		       }
		     } //if no track
		   } else if (partPdgId==22) { //photons
		     if( fabs((*iPart)->eta())<3.) { //tracker acceptance
		       nPhotonsGen += 1;
		       p4PhotonsGen += p4;
		       //save photons and later check for conversions:
		       //photon conversions switched OFF!
		   //shortPtcls.push_back(*iPart); //(only if eta<3)
		     } else {
		       nHFEMGen += 1;
		       p4HFEMGen += p4;
		     }
		   } else if ((fabs(partPdgId) != 12) && (fabs(partPdgId) != 14)
			      && (fabs(partPdgId) != 16)) { // veto neutrinos
		     
		     if( fabs((*iPart)->eta())<3. ) {
		       nNeutralHadronsGen += 1;
		       p4NeutralHadronsGen += p4;
		     } else {
		       nHFHadronsGen += 1;
		       p4HFHadronsGen += p4;
		     }
		     
		     // Decay K0S and Lambda later and correct fractions
		     if (abs(partPdgId)==310 || abs(partPdgId)==3122)
		       shortPtcls.push_back(*iPart);
		     
		     //save single neutral hadron components
		     if( abs(partPdgId)==2112 ) { //neutrons
		       nNeutronsGen += 1;
		       p4NeutronsGen += p4;
		     } else if( abs(partPdgId)==130 ) { //K0L's
		       nK0LGen += 1;
		       p4K0LGen += p4;
		     } else if( abs(partPdgId)==310 ) { //K0S's
		       nK0SGen += 1;
		       p4K0SGen += p4;
		     } else if( abs(partPdgId)==3122 ) { //lambda's
		       nLambdasGen += 1;
		       p4LambdasGen += p4;
		     } else if( abs(partPdgId)==3322 ) { //csi's
		       nCsiGen += 1;
		       p4CsiGen += p4;
		     } else { //not much else
		       nOtherNeutralHadronsGen += 1;
		       p4OtherNeutralHadronsGen += p4;
		     }
		 
		   } //if neutral hadrons
		   
		 } //for jetParticles
	     
		 // ------------------ BEGIN short decays (photons, K0S's and lambdas)
		 for (vector<const GenParticle*>::const_iterator iGen = shortPtcls.begin();
		      iGen != shortPtcls.end(); ++iGen) {
		   if( simTracks_h.isValid() ){
		     
		     // Check if GenParticle corresponds to a decayed SimTrack
		     if (decayedSims.find(*iGen) != decayedSims.end()) {
		       
		       // Check that the SimTrack Decay products were stored
		       const SimTrack *trk = decayedSims[*iGen];
		       if (promptDecays.find(trk) != promptDecays.end()) {
			 
			 // Convert track momentum to normal TLorentzVector, wrong type :(
			 math::XYZTLorentzVectorD const& p4t = (*iGen)->p4();
			 TLorentzVector p4mom(p4t.px(), p4t.py(), p4t.pz(), p4t.energy());
			 if ((*iGen)->pdgId()==22) {
			   nPhotonsGen -= 1;
			   p4PhotonsGen -= p4mom;
			 }
			 else {
			   if( fabs((*iGen)->eta())<3. ) {
			     nNeutralHadronsGen -= 1;
			     p4NeutralHadronsGen -= p4mom;
			   } else {
			     nHFHadronsGen -= 1;
			     p4HFHadronsGen -= p4mom;
			   }
			 }
			 
			 set<const SimTrack*> const& kids = promptDecays.find(trk)->second;
			 
			 for (set<const SimTrack*>::const_iterator iSim = kids.begin();
			      iSim != kids.end(); ++iSim) {
			   
			   // Convert track momentum to normal TLorentzVector, wrong type :(
			   math::XYZTLorentzVectorD const& p4t = (*iSim)->momentum();
			   TLorentzVector p4(p4t.px(), p4t.py(), p4t.pz(), p4t.energy());
			   
			   const SimVertex *vtx = promptVertex[*iSim];
			   
			   double vertR = vtx->position().Rho();
			   double vertZ = vtx->position().z();
			   //double trkEta = trk->momentum().eta();
			   bool decayedBeforeCalo =  ((vertR < 129.) && ( vertZ < 304.));
		       
			   // Check if the decay happened early enough for the
			   // charged track to be reconstructed
			   bool decayTrackable = (vertR < 30. && fabs(p4.Eta()) < 2.5
						  && (*iSim)->charge() != 0 && p4.Pt() > 0.075);
			   
			   if (decayedBeforeCalo && (*iSim)->type()==111) { //pizeros
			     if( fabs(p4.Eta())<3. ) {
			       nPhotonsGen += 2;  //both
			       p4PhotonsGen += p4;
			     } else {
			       nHFEMGen += 2;  //both
			       p4HFEMGen += p4;
			     }
			   }
			   else if (decayTrackable) {
			     if( fabs((*iSim)->type())==11 ) { //electrons
			       nElectronsGen += 1;
			       p4ElectronsGen += p4;
			     } else {
			       nTracksGen += 1;
			       p4TracksGen += p4;
			     }
			   } else {
			     //if( (*iGen)->pdgId()==22 ) { //photons //BEFORE: wrong??
			     if( (*iSim)->type()==22 || fabs((*iSim)->type())==11 ) { //photons or electrons
			       if( fabs(p4.Eta())<3. ) {
				 nPhotonsGen += 1;
				 p4PhotonsGen += p4;
			       } else {
				 nHFEMGen += 1;
				 p4HFEMGen += p4;
			       }
			     } else {
			       if( fabs(p4.Eta())<3. ) {
				 nNeutralHadronsGen += 1;
				 p4NeutralHadronsGen += p4;
			       } else {
				 nHFHadronsGen += 1;
				 p4HFHadronsGen += p4;
			       }
			     } //if-else photons
			   } //if-else decay trackable
			 } // for iSim loop on kids
		       } // has promptDecays
		     } // has decayedSims
		   }
		 } // for iGen
	  
		 // ------------------ END short decays
		 
		 const TLorentzVector *p = 0;
	     
		 nMuonsGen_akt5[nJetGen_akt5] = nMuonsGen;
		 p = &p4MuonsGen;
		 eMuonsGen_akt5[nJetGen_akt5] = p->E() / it->energy();
		 
		 nElectronsGen_akt5[nJetGen_akt5] = nElectronsGen;
		 p = &p4ElectronsGen;
		 eElectronsGen_akt5[nJetGen_akt5] = p->E() / it->energy();
		 
		 nPhotonsGen_akt5[nJetGen_akt5] = nPhotonsGen;
		 p = &p4PhotonsGen;
		 ePhotonsGen_akt5[nJetGen_akt5] = p->E() / it->energy();
		 
		 nTracksGen_akt5[nJetGen_akt5] = nTracksGen;
		 p = &p4TracksGen;
		 eTracksGen_akt5[nJetGen_akt5] = p->E() / it->energy();
		 
		 nNeutralHadronsGen_akt5[nJetGen_akt5] = nNeutralHadronsGen;
		 p = &p4NeutralHadronsGen;
		 eNeutralHadronsGen_akt5[nJetGen_akt5] = p->E() / it->energy();
		 
		 nHFHadronsGen_akt5[nJetGen_akt5] = nHFHadronsGen;
		 p = &p4HFHadronsGen;
		 eHFHadronsGen_akt5[nJetGen_akt5] = p->E() / it->energy();
		 
		 nHFEMGen_akt5[nJetGen_akt5] = nHFEMGen;
		 p = &p4HFEMGen;
		 eHFEMGen_akt5[nJetGen_akt5] = p->E() / it->energy();
		 
		 nNeutronsGen_akt5[nJetGen_akt5] = nNeutronsGen;
		 p = &p4NeutronsGen;
		 eNeutronsGen_akt5[nJetGen_akt5] = p->E() / it->energy();
		 
		 nK0LGen_akt5[nJetGen_akt5] = nK0LGen;
		 p = &p4K0LGen;
		 eK0LGen_akt5[nJetGen_akt5] = p->E() / it->energy();
		 
		 nK0SGen_akt5[nJetGen_akt5] = nK0SGen;
		 p = &p4K0SGen;
		 eK0SGen_akt5[nJetGen_akt5] = p->E() / it->energy();
		 
		 nLambdasGen_akt5[nJetGen_akt5] = nLambdasGen;
		 p = &p4LambdasGen;
		 eLambdasGen_akt5[nJetGen_akt5] = p->E() / it->energy();
		 
		 nCsiGen_akt5[nJetGen_akt5] = nCsiGen;
		 p = &p4CsiGen;
		 eCsiGen_akt5[nJetGen_akt5] = p->E() / it->energy();
		 
		 nOtherNeutralHadronsGen_akt5[nJetGen_akt5] = nOtherNeutralHadronsGen;
		 p = &p4OtherNeutralHadronsGen;
		 eOtherNeutralHadronsGen_akt5[nJetGen_akt5] = p->E() / it->energy();
	       }
	     nJetGen_akt5++;
	   } // if >genjetptthr
	 } // gen_akt5
       }

     if( dumpAKT7Jets_ || dumpJPTAKT7Jets_ || dumpPFAKT7Jets_ )
       {     
	 for (GenJetCollection::const_iterator it = jetsgenakt7->begin(); 
	      it != jetsgenakt7->end(); ++it) {
	   
	   if (nJetGen_akt7>=100) {cout << "number of gen jets antikt 07 is larger than 100. Skipping" << endl; continue;}
	   ptJetGen_akt7[nJetGen_akt7] = it->pt();	 
	   eJetGen_akt7[nJetGen_akt7] = it->energy();	 
	   etaJetGen_akt7[nJetGen_akt7] = it->eta();	 
	   phiJetGen_akt7[nJetGen_akt7] = it->phi();	      
	   
	   nJetGen_akt7++;
	 }
       }
   } //if(isMC)
   
     // Get the primary vertex coordinates
   for (VertexCollection::const_iterator it = VertexHandle->begin(); 
	it != VertexHandle->end(); ++it) {
     
     vx[nvertex] = (it->isValid()) ? it->x() : 999.;
     vy[nvertex] = (it->isValid()) ? it->y() : 999.;
     vz[nvertex] = (it->isValid()) ? it->z() : 999.;
     
     vntracks[nvertex] = (it->isValid()) ? it->tracksSize() : 0;
     vchi2[nvertex] = (it->isValid()) ? it->normalizedChi2() : 100.;
     vndof[nvertex] = (it->isValid()) ? it->ndof() : 0.;
     
     nvertex++;
   }
   

   //////////////////////////////////////
   /// ISOLATION TOOLS
   //////////////////////////////////////

    // compute track isolation w/o dz cut
    float isolationtrackThresholdA = 0.0;
    float TrackConeOuterRadiusA    = 0.4;
    float TrackConeInnerRadiusA    = 0.04;
    float isolationtrackEtaSliceA  = 0.015;
    float longImpactParameterA     = 99999.;
    float transImpactParameterA    = 0.1;
    float isolationtrackThresholdB = 0.0;
    float TrackConeOuterRadiusB    = 0.3;
    float TrackConeInnerRadiusB    = 0.04;
    float isolationtrackEtaSliceB  = 0.015;
    float longImpactParameterB     = 99999.;
    float transImpactParameterB    = 0.1;

   PhotonTkIsolation tkIsoCone04( TrackConeOuterRadiusA,
                                  TrackConeInnerRadiusA,
				  isolationtrackEtaSliceA,
				  isolationtrackThresholdA,
				  longImpactParameterA,
				  transImpactParameterA,
				  tracksHP.product(),
				  math::XYZPoint(vertexBeamSpot.x0(),vertexBeamSpot.y0(),vertexBeamSpot.z0())
				  );
        
   PhotonTkIsolation tkIsoCone03( TrackConeOuterRadiusB,
                                  TrackConeInnerRadiusB,
				  isolationtrackEtaSliceB,
				  isolationtrackThresholdB,
				  longImpactParameterB,
				  transImpactParameterB,
				  tracksHP.product(),
				  math::XYZPoint(vertexBeamSpot.x0(),vertexBeamSpot.y0(),vertexBeamSpot.z0())
				  );
        

    float isolationtrackThresholdA_newOpt = 0.0;
    float TrackConeOuterRadiusA_newOpt    = 0.4;
    float TrackConeInnerRadiusA_newOpt    = 0.02;
    float isolationtrackEtaSliceA_newOpt  = 0.0;
    float longImpactParameterA_newOpt     = 1.0;
    float transImpactParameterA_newOpt    = 0.1;
    float isolationtrackThresholdB_newOpt = 0.0;
    float TrackConeOuterRadiusB_newOpt    = 0.3;
    float TrackConeInnerRadiusB_newOpt    = 0.02;
    float isolationtrackEtaSliceB_newOpt  = 0.0;
    float longImpactParameterB_newOpt     = 1.0;
    float transImpactParameterB_newOpt    = 0.1;

   PhotonTkIsolation tkIsoCone04_newOpt( TrackConeOuterRadiusA_newOpt,
                                  TrackConeInnerRadiusA_newOpt,
        isolationtrackEtaSliceA_newOpt,
        isolationtrackThresholdA_newOpt,
        longImpactParameterA_newOpt,
        transImpactParameterA_newOpt,
        tracks.product(),
        math::XYZPoint(vertexBeamSpot.x0(),vertexBeamSpot.y0(),vertexBeamSpot.z0()) );
        
   PhotonTkIsolation tkIsoCone03_newOpt( TrackConeOuterRadiusB_newOpt,
                                  TrackConeInnerRadiusB_newOpt,
        isolationtrackEtaSliceB_newOpt,
        isolationtrackThresholdB_newOpt,
        longImpactParameterB_newOpt,
        transImpactParameterB_newOpt,
        tracks.product(),
        math::XYZPoint(vertexBeamSpot.x0(),vertexBeamSpot.y0(),vertexBeamSpot.z0()) );


   edm::ParameterSet remover02_config;
   remover02_config.addUntrackedParameter<double>("isolation_cone_size_forSCremoval",0.2);
   remover02_config.addUntrackedParameter<edm::InputTag>("tag_jets",edm::InputTag("ak5PFJetsCorrected"));
   SuperClusterFootprintRemoval remover02(iEvent,iSetup,remover02_config); 

   edm::ParameterSet remover03_config;
   remover03_config.addUntrackedParameter<double>("isolation_cone_size_forSCremoval",0.3);
   remover03_config.addUntrackedParameter<edm::InputTag>("tag_jets",edm::InputTag("ak5PFJetsCorrected"));
   SuperClusterFootprintRemoval remover03(iEvent,iSetup,remover03_config); 

   edm::ParameterSet remover04_config;
   remover04_config.addUntrackedParameter<double>("isolation_cone_size_forSCremoval",0.4);
   remover04_config.addUntrackedParameter<edm::InputTag>("tag_jets",edm::InputTag("ak5PFJetsCorrected"));
   SuperClusterFootprintRemoval remover04(iEvent,iSetup,remover04_config); 

   EcalClusterLazyTools lazyTools(iEvent, iSetup, edm::InputTag("reducedEcalRecHitsEB"),
				  edm::InputTag("reducedEcalRecHitsEE"), edm::InputTag("reducedEcalRecHitsES") );
    CaloSubdetectorTopology *topology_p = new EcalPreshowerTopology(geoHandle);

   // pass the collection to the ID calculator
   cicPhotonId->configure(VertexHandle, tracks, ElectronHandle, PFCandidates, rhoAllJets); 

   for (PhotonCollection::const_iterator it = PhotonHandle->begin(); 
	it != PhotonHandle->end(); ++it) {
     
     if (it->energy()<3.) continue;
     if (nPhot>=40) {cout << "number of photons is larger than 40. Skipping" << endl; continue;}
     
     ptPhot[nPhot] = it->pt();
     ePhot[nPhot] = it->energy();	 
     escPhot[nPhot] = it->superCluster()->energy();	 

     // new for 52X
     std::pair<double,double> cor = ecorr_->CorrectedEnergyWithErrorV3(*it,*VertexHandle,*rhoAllH,lazyTools,iSetup,true); //rescaling input variables for sigmaE calculation
     escRegrPhot[nPhot] = cor.first;        
     escRegrPhotError[nPhot] = cor.second;  

     escPhFixPhot[nPhot] = -999.;
     escPhFixPhotError[nPhot] = -999.;

     escRawPhot[nPhot] = it->superCluster()->rawEnergy();	 
     eseedPhot[nPhot] = it->superCluster()->seed()->energy();	 

     etaPhot[nPhot] = it->eta();	 
     phiPhot[nPhot] = it->phi();	      
//      etascPhot[nPhot] = it->superCluster()->eta();	 
//      phiscPhot[nPhot] = it->superCluster()->phi();	      
//      xscPhot[nPhot] = it->superCluster()->position().x();	 
//      yscPhot[nPhot] = it->superCluster()->position().y();	      
//      zscPhot[nPhot] = it->superCluster()->position().z();	      

     etascPhot[nPhot] = it->superCluster()->position().eta();	 
     phiscPhot[nPhot] = it->superCluster()->position().phi();	      

     xscPhot[nPhot] = it->superCluster()->position().x();	 
     yscPhot[nPhot] = it->superCluster()->position().y();	      
     zscPhot[nPhot] = it->superCluster()->position().z();	      

     xcaloPhot[nPhot] = it->caloPosition().x();	 
     ycaloPhot[nPhot] = it->caloPosition().y();	      
     zcaloPhot[nPhot] = it->caloPosition().z();	      

     hasPixelSeedPhot[nPhot] = it->hasPixelSeed();
     isEBPhot[nPhot] = it->isEB();
     isEEPhot[nPhot] = it->isEE();
     isEBEEGapPhot[nPhot] = it->isEBEEGap();

     reco::ConversionRef conv = 
       ConversionTools::matchedConversion(*(it->superCluster()),hConversions,vertexBeamSpot.position());
     hasMatchedConvPhot[nPhot] = conv.isNonnull();

     isValidVtxConvPhot[nPhot] = 0;
     nTracksConvPhot[nPhot]                    =   0; 
     pairInvariantMassConvPhot[nPhot]          =   -9999.; 
     pairCotThetaSeparationConvPhot[nPhot]     =   -9999.; 
     pairMomentum_xConvPhot[nPhot]             =   -9999.; 
     pairMomentum_yConvPhot[nPhot]             =   -9999.; 
     pairMomentum_zConvPhot[nPhot]             =   -9999.; 
     chi2ConvPhot[nPhot]                       =   -9999.; 
     nDofConvPhot[nPhot]                       =   -9999.; 
     conv_vxConvPhot[nPhot]                    =   -9999.; 
     conv_vyConvPhot[nPhot]                    =   -9999.; 
     conv_vzConvPhot[nPhot]                    =   -9999.; 
     eOverPConvPhot[nPhot]                     =   -9999.; 
     distOfMinimumApproachConvPhot[nPhot]      =   -9999.; 
     dPhiTracksAtVtxConvPhot[nPhot]            =   -9999.; 
//      dPhiTracksAtEcalConvPhot[nPhot]           =   -9999.; 
//      dEtaTracksAtEcalConvPhot[nPhot]           =   -9999.; 

     if ( conv.isNonnull() ){
       isValidVtxConvPhot[nPhot] = conv->conversionVertex().isValid();
       nTracksConvPhot[nPhot]                    = conv->nTracks();
       pairInvariantMassConvPhot[nPhot]          = conv->pairInvariantMass();
       pairCotThetaSeparationConvPhot[nPhot]     = conv->pairCotThetaSeparation();
       pairMomentum_xConvPhot[nPhot]             = conv->pairMomentum().x();
       pairMomentum_yConvPhot[nPhot]             = conv->pairMomentum().y();
       pairMomentum_zConvPhot[nPhot]             = conv->pairMomentum().z();
       chi2ConvPhot[nPhot]                       = conv->conversionVertex().chi2();
       nDofConvPhot[nPhot]                       = conv->conversionVertex().ndof();
       conv_vxConvPhot[nPhot]                    = conv->conversionVertex().x();
       conv_vyConvPhot[nPhot]                    = conv->conversionVertex().y();
       conv_vzConvPhot[nPhot]                    = conv->conversionVertex().z();
       eOverPConvPhot[nPhot]                     = conv->EoverP();
       distOfMinimumApproachConvPhot[nPhot]      = conv->distOfMinimumApproach();
       dPhiTracksAtVtxConvPhot[nPhot]            = conv->dPhiTracksAtVtx();
//        dPhiTracksAtEcalConvPhot[nPhot]           = conv->dPhiTracksAtEcal();
//        dEtaTracksAtEcalConvPhot[nPhot]           = conv->dEtaTracksAtEcal();
     }

     hasMatchedPromptElePhot[nPhot]    = ConversionTools::hasMatchedPromptElectron(it->superCluster(), 
                                          ElectronHandle, hConversions, vertexBeamSpot.position());
     
     
     const Ptr<CaloCluster> theSeed = it->superCluster()->seed(); 
     
     const EBRecHitCollection* rechits = ( it->isEB()) ? rhitseb : rhitsee;

     // photon ID (spike ID) related info:
     // timing:
     std::pair<DetId, float> maxRH = EcalClusterTools::getMaximum( *theSeed, &(*rechits) );
     DetId seedCrystalId = maxRH.first;
     EcalRecHitCollection::const_iterator seedRH = rechits->find(seedCrystalId);
     if(maxRH.second) timePhot[nPhot] = (float)seedRH->time();
     else timePhot[nPhot] = 0;
     // swiss cross:
     e4SwissCrossPhot[nPhot] = (!it->isEB()) ? 0. :
       ( EcalClusterTools::eLeft( *theSeed, &(*rechits), topology ) +
	 EcalClusterTools::eRight( *theSeed, &(*rechits), topology ) +
	 EcalClusterTools::eTop( *theSeed, &(*rechits), topology ) +
	 EcalClusterTools::eBottom( *theSeed, &(*rechits), topology ) );
     
     double ptiso004(0.), ptiso035(0.), ptiso04(0.);
     int ntrkiso004(0), ntrkiso035(0), ntrkiso04(0);



     // Fill default photon ID variables
     
     pid_jurECAL[nPhot] = it->ecalRecHitSumEtConeDR04();//isolationEcalRecHit
     //     pid_twrHCAL[nPhot] = it->hcalTowerSumEtConeDR04();//isolationHcalRecHit
     pid_twrHCAL[nPhot] = it->hcalTowerSumEtConeDR04() + (it->hadronicOverEm() - it->hadTowOverEm())*it->superCluster()->energy()/cosh(it->superCluster()->eta()); //new H/E 2012
     pid_HoverE[nPhot] = it->hadTowOverEm(); //new H/E for 2012 data
     //     pid_HoverE[nPhot] = it->hadronicOverEm();
     pid_hlwTrack[nPhot] = it->trkSumPtHollowConeDR04();//isolationHollowTrkCone
     pid_etawid[nPhot] = it->sigmaIetaIeta();//sqrt(it->covEtaEta());
     pid_hlwTrackNoDz[nPhot] = tkIsoCone04.getPtTracks( &(*it) );

     pid_scetawid[nPhot]       = it->superCluster()->etaWidth();
     pid_scphiwid[nPhot]        = it->superCluster()->phiWidth();

     pid_jurECAL03[nPhot] = it->ecalRecHitSumEtConeDR03();//isolationEcalRecHit
     //pid_twrHCAL03[nPhot] = it->hcalTowerSumEtConeDR03();//isolationHcalRecHit
     pid_twrHCAL03[nPhot] = it->hcalTowerSumEtConeDR03() + (it->hadronicOverEm() - it->hadTowOverEm())*it->superCluster()->energy()/cosh(it->superCluster()->eta()); //new H/E 2012
     pid_hlwTrack03[nPhot] = it->trkSumPtHollowConeDR03();//isolationHollowTrkCone
     pid_hlwTrack03NoDz[nPhot] = tkIsoCone03.getPtTracks( &(*it) );

     //Dumping trackIsolation wrt vtx hypothesis 
     for (unsigned int ivtx=0;ivtx<VertexHandle->size();++ivtx)
       {
	 reco::Photon aPhot((*it));
	 aPhot.setVertex(reco::Vertex::Point((*VertexHandle)[ivtx].x(),(*VertexHandle)[ivtx].y(),(*VertexHandle)[ivtx].z()));
	 pid_hlwTrackForCiC[nPhot][ivtx]=tkIsoCone04_newOpt.getPtTracks( &aPhot );
	 pid_hlwTrack03ForCiC[nPhot][ivtx]=tkIsoCone03_newOpt.getPtTracks( &aPhot );
       }
 
     // fill default ID variables from
     // https://twiki.cern.ch/twiki/bin/view/CMS/PhotonIDAnalysis for 31X
     // (followed from https://twiki.cern.ch/twiki/bin/view/CMS/PhotonID)
     // using reco equivalents to pat::Photon methods

     if (it->isEB()) {

       pid_isEM[nPhot] =    (pid_jurECAL[nPhot]  < 5 + 0.15*ptPhot[nPhot] &&
			     pid_twrHCAL[nPhot]  < 5 &&
			     pid_HoverE[nPhot]   < 0.5);
       pid_isLoose[nPhot] = (pid_isEM[nPhot] &&
			     pid_jurECAL[nPhot]  < 5 + 0.0045*ptPhot[nPhot] &&
			     pid_HoverE[nPhot]   < 0.15 &&
			     pid_hlwTrack[nPhot] < 9);
       pid_isTight[nPhot] = (pid_isLoose[nPhot] &&
			     pid_etawid[nPhot]   < 0.013);
     }
     else {

       pid_isEM[nPhot] =    (pid_jurECAL[nPhot]  < 5 + 0.15*ptPhot[nPhot] &&
			     pid_twrHCAL[nPhot]  < 7 &&
			     pid_HoverE[nPhot]   < 0.5);
       pid_isLoose[nPhot] = (pid_isEM[nPhot] &&
			     pid_jurECAL[nPhot]  < 5 + 0.02*ptPhot[nPhot] &&
			     pid_HoverE[nPhot]   < 0.15 &&
			     pid_hlwTrack[nPhot] < 9);
       pid_isTight[nPhot] = (pid_isLoose[nPhot] &&
			     pid_etawid[nPhot]   < 0.026); // fixed default
     }

     // Need to reconstruct PAT photon to get the decisions
     // nah, doesn't fill isolation automatically, need to run some produces
     //pat::Photon pho(*it);
     //const vector<pat::Photon::IdPair> &phoids = pho.photonIDs();
     //pid_isEM[nPhot] = false;
     //pid_isLoose[nPhot] = false;
     //pid_isTight[nPhot] = false;
     //for (unsigned int i = 0; i != phoids.size(); ++i) {
     //if (phoids[i].first=="EM") pid_isEM[nPhot] = phoids[i].second;
     //if (phoids[i].first=="Loose") pid_isLoose[nPhot] = phoids[i].second;
     //if (phoids[i].first=="Tight") pid_isTight[nPhot] = phoids[i].second;
     //}

     // PAT equivalents of photon ID variables
     //if (_debug) {
     //if (!(fabs(pid_jurECAL[nPhot]-pho.ecalIso())<1e-3)) {
     // cerr << Form("jurECAL %1.4g, ecaliso %1.4g",
     //	      pid_jurECAL[nPhot],pho.ecalIso()) << endl << flush;
     // //assert(fabs(pid_jurECAL[nPhot]-pho.ecalIso())<1e-3);
     //}
     //if (!(fabs(pid_twrHCAL[nPhot]-pho.hcalIso())<1e-3)) {
     // cerr << Form("twrHCAL %1.4g, hcaliso %1.4g",
     //	      pid_twrHCAL[nPhot],pho.hcalIso()) << endl << flush;
     // //assert(fabs(pid_twrHCAL[nPhot]-pho.hcalIso())<1e-3);
     //}
     //if (!(fabs(pid_hlwTrack[nPhot]-pho.trackIso())<1e-3)) {
     // cerr << Form("hlwTrack %1.4g, trackiso %1.4g",
     //	      pid_hlwTrack[nPhot],pho.trackIso()) << endl << flush;
     // //assert(fabs(pid_hlwTrack[nPhot]-pho.trackIso())<1e-3);
     //}
     //}

     // calculate track isolation for different cone values


     for (TrackCollection::const_iterator itTrack = tracks->begin();
	 itTrack != tracks->end(); ++itTrack) {

       double etaTrack = itTrack->eta();
       double phiTrack = itTrack->phi();

       double deltaPhi = phiTrack-it->phi();
       double deltaEta = etaTrack-it->eta();
       if (deltaPhi > Geom::pi()) deltaPhi -= 2.*Geom::pi();
       if (deltaPhi < -Geom::pi()) deltaPhi += 2.*Geom::pi();
       double deltaR = std::sqrt(deltaEta*deltaEta+deltaPhi*deltaPhi);

       if (deltaR < .04)  {ptiso004  += sqrt(itTrack->pt()); ntrkiso004++; }
       if (deltaR < .35)   {ptiso035 += sqrt(itTrack->pt()); ntrkiso035++;}
       if (deltaR < .4)   {ptiso04 += sqrt(itTrack->pt()); ntrkiso04++;}
    
     }

     ptiso004Phot[nPhot] = ptiso004;
     ntrkiso004Phot[nPhot] = ntrkiso004;
     ptiso035Phot[nPhot] = ptiso035;
     ntrkiso035Phot[nPhot] = ntrkiso035;
     ptiso04Phot[nPhot] = ptiso04;
     ntrkiso04Phot[nPhot] = ntrkiso04;

     // calculate HCAL isolation
  
     double hcalEnergy = 0;
     reco::SuperClusterRef sc = it->get<reco::SuperClusterRef>();
     CaloConeSelector selector4(0.4, geometry, DetId::Hcal); 
     //     cout << "sc eta: " << sc->eta() << "    and phi : " << sc->phi() << endl;
     std::auto_ptr<CaloRecHitMetaCollectionV> selected = selector4.select(sc->eta(),sc->phi(),mhbhe); 
     for (CaloRecHitMetaCollectionV::const_iterator hit=selected->begin(); hit != selected->end(); ++hit){
       hcalEnergy += hit->energy(); 
     }
     hcalovecal04Phot[nPhot] = hcalEnergy/it->energy(); 

   h1_etaPhot->Fill( sc->eta() );
   for( CaloRecHitMetaCollectionV::const_iterator ihbherh=selected->begin(); ihbherh!=selected->end(); ++ihbherh ) {
     h1_hbherh_detid->Fill( ihbherh->detid().subdetId() );
   }
   h2_n_vs_eta->Fill( fabs(sc->eta()) , selected->size() );
   
   CaloClusterPtr SCseed = it->superCluster()->seed();
   
   // calculate ECAL isolation 
   
   // ecal isolation with SC seed rechits removal
   SuperClusterHitsEcalIsolation scBasedIsolation(rhitseb,rhitsee);
   scBasedIsolation.setExtRadius(0.1);
   scBasedIsolation.excludeHalo(false);
   scBasedIsolation.setExtRadius(0.4);
   ecaliso04Phot[nPhot]  = scBasedIsolation.getSum(iEvent,iSetup,&(*SCseed));
   
   // cluster shape variables
   
   //     if (it->isEB()){
   if(maxRH.second) {
     Cluster2ndMoments moments = EcalClusterTools::cluster2ndMoments(*SCseed, *rechits);
     std::vector<float> etaphimoments = EcalClusterTools::localCovariances(*SCseed, &(*rechits), &(*topology));
     sMajMajPhot[nPhot]=moments.sMaj;
     sMinMinPhot[nPhot]=moments.sMin;
     alphaPhot[nPhot]=moments.alpha;
     sEtaEtaPhot[nPhot]=etaphimoments[0];
     sEtaPhiPhot[nPhot]=etaphimoments[1];
     sPhiPhiPhot[nPhot]=etaphimoments[2];
     float lambdaMinus         = (etaphimoments[0] + etaphimoments[2] - sqrt(pow(etaphimoments[0] - etaphimoments[2], 2) + 4*pow(etaphimoments[1], 2)));
     float lambdaPlus          = (etaphimoments[0] + etaphimoments[2] + sqrt(pow(etaphimoments[0] - etaphimoments[2], 2) + 4*pow(etaphimoments[1], 2)));
     pid_lambdaRatio[nPhot]    = lambdaMinus/lambdaPlus;
   }else{
     sMajMajPhot[nPhot]=-100.;
     sMinMinPhot[nPhot]=-100.;
     alphaPhot[nPhot]=-100.;
     sEtaEtaPhot[nPhot] = it->sigmaEtaEta();//-100.;
     sEtaPhiPhot[nPhot]=-100.;
     sPhiPhiPhot[nPhot]=-100.;
     pid_lambdaRatio[nPhot]   = -100.;
   }
   
   // ES variables
   pid_esXwidth[nPhot] = 0.;
   pid_esYwidth[nPhot] = 0.;
   
   if (ecalhitses.isValid() && ( fabs(it->superCluster()->eta()) > 1.6 && fabs(it->superCluster()->eta()) < 3 ) ) {
     std::vector<float> phoESHits0 = getESHits(it->superCluster()->x(), it->superCluster()->y(), it->superCluster()->z(), rechits_map_, *geometry, topology_p, 0);
     std::vector<float> phoESShape = getESShape(phoESHits0, &(pid_esXShape[nPhot][0]), &(pid_esYShape[nPhot][0]) );
     pid_esXwidth[nPhot] = phoESShape[0];
     pid_esYwidth[nPhot] = phoESShape[1];
   }

   E1Phot[nPhot] = maxRH.second;//-100.;
   E2OverE9Phot[nPhot]=-99.;
   if (it->isEB())
     E2OverE9Phot[nPhot] = GetE2OverE9(seedCrystalId,*rechits);
   E9Phot[nPhot] = it->e3x3();//-100.;
   E25Phot[nPhot] = it->e5x5();//-100.;
   E4Phot[nPhot] = lazyTools.e2x2(*(it->superCluster()->seed()));//-100.;
   
   if (_debug && it->isEB() &&
       fabs(it->maxEnergyXtal()-387)>1 && // weird bug in 31X?
       (fabs(sEtaEtaPhot[nPhot]-it->sigmaEtaEta())>1e-4 ||
	//fabs(E1Phot[nPhot]-it->maxEnergyXtal())>1e-4 ||
	fabs(E1Phot[nPhot]-SCseed->energy())>1e-4 ||
	fabs(E9Phot[nPhot]-it->e3x3())>1e-4 ||
	fabs(E25Phot[nPhot]-it->e5x5())>1e-4)) {
     cout << Form("Eta %1.3g\n"
		  "sEtaEta %1.3g vs %1.3f\n"
		  "E1      %1.3g vs %1.3g\n"
		  "E9      %1.3g vs %1.3g\n"
		  "E25     %1.3g vs %1.3g\n",
		  SCseed->eta(),
		  sEtaEtaPhot[nPhot], it->sigmaEtaEta(),
		  E1Phot[nPhot], SCseed->energy(),
		  E9Phot[nPhot], it->e3x3(),
		  E25Phot[nPhot], it->e5x5());
   }
   
   ieleassocPhot[nPhot] = -999; 
   pid_deltaRToTrackPhot[nPhot] = 99.; 
   
   reco::PhotonRef localPho(PhotonHandle, nPhot);
   // isolation variables
   std::vector<reco::PFCandidate::ParticleType> temp;
   temp.push_back(reco::PFCandidate::gamma);
//    pid_pfIsoPhotons01ForCiC[nPhot]  = cicPhotonId->pfEcalIso(localPho, 0.1, 0.0, 0.070, 0.015, 0.0, 0.0, 0.0, reco::PFCandidate::gamma);
   pid_pfIsoPhotons02ForCiC[nPhot]  = cicPhotonId->pfEcalIso(localPho, 0.2, 0.0, 0.070, 0.015, 0.0, 0.0, 0.0, reco::PFCandidate::gamma);
   pid_pfIsoPhotons03ForCiC[nPhot]  = cicPhotonId->pfEcalIso(localPho, 0.3, 0.0, 0.070, 0.015, 0.0, 0.0, 0.0, reco::PFCandidate::gamma);
   pid_pfIsoPhotons04ForCiC[nPhot]  = cicPhotonId->pfEcalIso(localPho, 0.4, 0.0, 0.070, 0.015, 0.0, 0.0, 0.0, reco::PFCandidate::gamma);
//    pid_pfIsoPhotons05ForCiC[nPhot]  = cicPhotonId->pfEcalIso(localPho, 0.5, 0.0, 0.070, 0.015, 0.0, 0.0, 0.0, reco::PFCandidate::gamma);
//    pid_pfIsoPhotons06ForCiC[nPhot]  = cicPhotonId->pfEcalIso(localPho, 0.6, 0.0, 0.070, 0.015, 0.0, 0.0, 0.0, reco::PFCandidate::gamma);
   
   temp.clear();
   temp.push_back(reco::PFCandidate::h0);
   // Custom Egamma and noveto are the same
//    pid_pfIsoNeutrals01ForCiC[nPhot] = cicPhotonId->pfHcalIso(localPho, 0.1, 0.00, reco::PFCandidate::h0);
   pid_pfIsoNeutrals02ForCiC[nPhot] = cicPhotonId->pfHcalIso(localPho, 0.2, 0.00, reco::PFCandidate::h0);
   pid_pfIsoNeutrals03ForCiC[nPhot] = cicPhotonId->pfHcalIso(localPho, 0.3, 0.00, reco::PFCandidate::h0);
   pid_pfIsoNeutrals04ForCiC[nPhot] = cicPhotonId->pfHcalIso(localPho, 0.4, 0.00, reco::PFCandidate::h0);
//    pid_pfIsoNeutrals05ForCiC[nPhot] = cicPhotonId->pfHcalIso(localPho, 0.5, 0.00, reco::PFCandidate::h0);
//    pid_pfIsoNeutrals06ForCiC[nPhot] = cicPhotonId->pfHcalIso(localPho, 0.6, 0.00, reco::PFCandidate::h0);
   
   temp.clear();
   temp.push_back(reco::PFCandidate::h);
//    std::vector<float> pid_pfIsoCharged01;
   std::vector<float> pid_pfIsoCharged02;
   std::vector<float> pid_pfIsoCharged03;
   std::vector<float> pid_pfIsoCharged04;
//    std::vector<float> pid_pfIsoCharged05;
//    std::vector<float> pid_pfIsoCharged06;
//    pid_pfIsoCharged01=cicPhotonId->pfTkIsoWithVertex(localPho, 0.1, 0.02, 0.02, 0.0, 0.2, 0.1, reco::PFCandidate::h); 
   pid_pfIsoCharged02=cicPhotonId->pfTkIsoWithVertex(localPho, 0.2, 0.02, 0.02, 0.0, 0.2, 0.1, reco::PFCandidate::h); 
   pid_pfIsoCharged03=cicPhotonId->pfTkIsoWithVertex(localPho, 0.3, 0.02, 0.02, 0.0, 0.2, 0.1, reco::PFCandidate::h); 
   pid_pfIsoCharged04=cicPhotonId->pfTkIsoWithVertex(localPho, 0.4, 0.02, 0.02, 0.0, 0.2, 0.1, reco::PFCandidate::h); 
//    pid_pfIsoCharged05=cicPhotonId->pfTkIsoWithVertex(localPho, 0.5, 0.02, 0.02, 0.0, 0.2, 0.1, reco::PFCandidate::h); 
//    pid_pfIsoCharged06=cicPhotonId->pfTkIsoWithVertex(localPho, 0.6, 0.02, 0.02, 0.0, 0.2, 0.1, reco::PFCandidate::h); 

   assert (
//  	   (pid_pfIsoCharged01.size() == VertexHandle->size()) &&
  	   (pid_pfIsoCharged02.size() == VertexHandle->size()) &&
	   (pid_pfIsoCharged03.size() == VertexHandle->size()) &&
	   (pid_pfIsoCharged04.size() == VertexHandle->size()) //&&
//  	   (pid_pfIsoCharged05.size() == VertexHandle->size()) &&
//  	   (pid_pfIsoCharged06.size() == VertexHandle->size())
	   );
   
   for (unsigned int iiv=0;iiv<VertexHandle->size();++iiv)
     {
//        pid_pfIsoCharged01ForCiC[nPhot][iiv]=pid_pfIsoCharged01[iiv]; 
       pid_pfIsoCharged02ForCiC[nPhot][iiv]=pid_pfIsoCharged02[iiv]; 
       pid_pfIsoCharged03ForCiC[nPhot][iiv]=pid_pfIsoCharged03[iiv]; 
       pid_pfIsoCharged04ForCiC[nPhot][iiv]=pid_pfIsoCharged04[iiv]; 
//        pid_pfIsoCharged05ForCiC[nPhot][iiv]=pid_pfIsoCharged05[iiv]; 
//        pid_pfIsoCharged06ForCiC[nPhot][iiv]=pid_pfIsoCharged06[iiv]; 
     }
   
   double minDR = 1000.;       

   PFIsolation_struct pfIsoFPR02=remover02.PFIsolation(it->superCluster(),edm::Ptr<reco::Vertex>(VertexHandle,0));
   pid_pfIsoFPRCharged02[nPhot] = pfIsoFPR02.chargediso_primvtx; 
   pid_pfIsoFPRNeutral02[nPhot] = pfIsoFPR02.neutraliso; 
   pid_pfIsoFPRPhoton02[nPhot] = pfIsoFPR02.photoniso; 

   PFIsolation_struct pfIsoFPR03=remover03.PFIsolation(it->superCluster(),edm::Ptr<reco::Vertex>(VertexHandle,0));
   pid_pfIsoFPRCharged03[nPhot] = pfIsoFPR03.chargediso_primvtx; 
   pid_pfIsoFPRNeutral03[nPhot] = pfIsoFPR03.neutraliso; 
   pid_pfIsoFPRPhoton03[nPhot] = pfIsoFPR03.photoniso; 

   PFIsolation_struct pfIsoFPR04=remover04.PFIsolation(it->superCluster(),edm::Ptr<reco::Vertex>(VertexHandle,0));
   pid_pfIsoFPRCharged04[nPhot] = pfIsoFPR04.chargediso_primvtx; 
   pid_pfIsoFPRNeutral04[nPhot] = pfIsoFPR04.neutraliso; 
   pid_pfIsoFPRPhoton04[nPhot] = pfIsoFPR04.photoniso; 

   pid_pfIsoFPRRandomConeCharged02[nPhot] = pfIsoFPR02.chargediso_primvtx_rcone; 
   pid_pfIsoFPRRandomConeNeutral02[nPhot] = pfIsoFPR02.neutraliso_rcone;
   pid_pfIsoFPRRandomConePhoton02[nPhot] = pfIsoFPR02.photoniso_rcone;
   pid_pfIsoFPRRandomConeEta02[nPhot] = pfIsoFPR02.eta_rcone;
   pid_pfIsoFPRRandomConePhi02[nPhot] = pfIsoFPR02.phi_rcone;

   pid_pfIsoFPRRandomConeCharged03[nPhot] = pfIsoFPR03.chargediso_primvtx_rcone; 
   pid_pfIsoFPRRandomConeNeutral03[nPhot] = pfIsoFPR03.neutraliso_rcone;
   pid_pfIsoFPRRandomConePhoton03[nPhot] = pfIsoFPR03.photoniso_rcone;
   pid_pfIsoFPRRandomConeEta03[nPhot] = pfIsoFPR03.eta_rcone;
   pid_pfIsoFPRRandomConePhi03[nPhot] = pfIsoFPR03.phi_rcone;

   pid_pfIsoFPRRandomConeCharged04[nPhot] = pfIsoFPR04.chargediso_primvtx_rcone; 
   pid_pfIsoFPRRandomConeNeutral04[nPhot] = pfIsoFPR04.neutraliso_rcone;
   pid_pfIsoFPRRandomConePhoton04[nPhot] = pfIsoFPR04.photoniso_rcone;
   pid_pfIsoFPRRandomConeEta04[nPhot] = pfIsoFPR04.eta_rcone;
   pid_pfIsoFPRRandomConePhi04[nPhot] = pfIsoFPR04.phi_rcone;


     const reco::GsfElectron* bestMatchedPromptEle=0;
     for (GsfElectronCollection::const_iterator itElectron = ElectronHandle->begin();
	  itElectron != ElectronHandle->end(); ++itElectron) {
       

       if (itElectron->superCluster() != it->superCluster())
	 continue;
       
       if (itElectron->pt()<2.5)
	 continue;
       if (itElectron->gsfTrack()->trackerExpectedHitsInner().numberOfHits()>0)
	 continue;
       
//        reco::ConversionRef conv = 
// 	 ConversionTools::matchedConversion(*itElectron,hConversions,vertexBeamSpot.position());
       
//        if (conv.isNonnull())
// 	 continue;
       
       float deta=itElectron->deltaEtaSuperClusterTrackAtVtx(); 
       float dphi=itElectron->deltaPhiSuperClusterTrackAtVtx(); 
       float DR = sqrt ( deta*deta + dphi*dphi);
       if (DR<minDR)
	 {
	   minDR=DR;
	   bestMatchedPromptEle = &(*itElectron);
	 }
       
     }
     
     
    if(bestMatchedPromptEle)
       {
	 ieleassocPhot[nPhot] = nElePhot;
	 pid_deltaRToTrackPhot[nPhot] = minDR; 
      
	 pid_jurECALElePhot[nElePhot] = bestMatchedPromptEle->dr03EcalRecHitSumEt(); 
	 //	 pid_twrHCALElePhot[nElePhot] = bestMatchedPromptEle->dr03HcalTowerSumEt(); 
	 pid_twrHCALElePhot[nElePhot] = bestMatchedPromptEle->dr03HcalDepth1TowerSumEtBc();
	 //pid_HoverEElePhot[nElePhot] = bestMatchedPromptEle->hadronicOverEm(); 
	 pid_HoverEElePhot[nElePhot] = bestMatchedPromptEle->hcalOverEcalBc();
	 pid_hlwTrackElePhot[nElePhot] = bestMatchedPromptEle->dr03TkSumPt(); 
	 pid_etawidElePhot[nElePhot] = bestMatchedPromptEle->sigmaIetaIeta(); 
	 pid_dphivtxElePhot[nElePhot] = bestMatchedPromptEle->deltaPhiSuperClusterTrackAtVtx(); 
	 pid_detavtxElePhot[nElePhot] = bestMatchedPromptEle->deltaEtaSuperClusterTrackAtVtx(); 
	 pid_mishitsElePhot[nElePhot] = bestMatchedPromptEle->gsfTrack()->trackerExpectedHitsInner().numberOfHits(); 
	 ConversionInfo convInfo = convFinder.getConversionInfo(*bestMatchedPromptEle,tracks, bfield);      
	 pid_distElePhot[nElePhot] = (convInfo.dist() == -9999.? 9999:convInfo.dist());
	 pid_dcotElePhot[nElePhot] = (convInfo.dcot() == -9999.? 9999:convInfo.dcot());
	 pid_ptElePhot[nElePhot] = bestMatchedPromptEle->gsfTrack()->pt(); 

	 nElePhot++;       
       }	 
    ++nPhot;
   }



   std::vector<int> rankprod;  
   // Loop over reco photons to identify photons for vertex selection
   if (_debug)
     std::cout << "============= Starting vertex Analysis =================" << std::endl;

   std::map<int,PhotonInfo> photonsForVertexComputation;
   int iphot=0;
   for (PhotonCollection::const_iterator it = PhotonHandle->begin(); 
	it != PhotonHandle->end(); ++it) 
     {
       //pass preselection cuts
       
       //        bool hasMatchedPromptEle = ConversionTools::hasMatchedPromptElectron(it->superCluster(), 
       // 									    ElectronHandle, hConversions, vertexBeamSpot.position());
       //        if (hasMatchedPromptEle)
       // 	 continue;
       
       if(  it->pt()<25.  || ( fabs(it->superCluster()->position().eta())>1.4442 && fabs(it->superCluster()->position().eta())<1.566 ) || fabs(it->superCluster()->position().eta()) >2.5)
	 {
// 	   std::cout << "Photon skipped because of acceptance " << std::endl;
	   ++iphot;
	   continue;         
	 }
//        bool isEB = it->isEB(); 


//        //       if( it->ecalRecHitSumEtConeDR03() >= 8. || ((isEB && it->sigmaIetaIeta()>0.013) || (!isEB && it->sigmaIetaIeta()>0.031))  || it->hadronicOverEm() > 0.15 ) 
//        if( ((isEB && it->sigmaIetaIeta()>0.013) || (!isEB && it->sigmaIetaIeta()>0.03))  || it->hadronicOverEm() > 0.15 ) 
// 	 {
// // 	   std::cout << "Photon skipped because of preselection" << std::endl;
// 	   ++iphot;
// 	   continue;         
// 	 }

//     New MIT preselection that matches the trigger selection
       if (!PhotonMITPreSelection(iphot,1))
	 {
	   ++iphot;
	   continue;         
	 }

       reco::ConversionRef conv = 
	 ConversionTools::matchedConversion(*(it->superCluster()),hConversions,vertexBeamSpot.position());
       if (conv.isNonnull())
	 {
	   
	   reco::Vertex vtx=conv->conversionVertex();
	   // photonsForVertexComputation.insert(std::make_pair<int,PhotonInfo>(iphot,PhotonInfo(it - PhotonHandle->begin(),
	   photonsForVertexComputation.insert(std::make_pair(iphot,PhotonInfo(it - PhotonHandle->begin(),
									      TVector3(it->superCluster()->position().x(),it->superCluster()->position().y(),it->superCluster()->position().z()),
									      TVector3(vertexBeamSpot.x0(),vertexBeamSpot.y0(),vertexBeamSpot.z0()),
									      TVector3(vtx.x(), vtx.y(), vtx.z()),
									      TVector3(conv->refittedPairMomentum().x(), conv->refittedPairMomentum().y(), conv->refittedPairMomentum().z()),
									      it->energy(),
									      it->isEB(),
									      conv->nTracks(),
									      conv->conversionVertex().isValid(),
									      TMath::Prob(vtx.chi2(), vtx.ndof()),
									      -999.) ) );
	 }
       else
	 {
	   // photonsForVertexComputation.insert(std::make_pair<int,PhotonInfo>(iphot,PhotonInfo(it - PhotonHandle->begin(),
	   photonsForVertexComputation.insert(std::make_pair(iphot,PhotonInfo(it - PhotonHandle->begin(),
									      TVector3(it->superCluster()->position().x(),it->superCluster()->position().y(),it->superCluster()->position().z()),
									      TVector3(vertexBeamSpot.x0(),vertexBeamSpot.y0(),vertexBeamSpot.z0()),
									      TVector3(-999.,-999.,-999),
									      TVector3(-999.,-999.,-999),
									      it->energy(),
									      it->isEB(),
									      -999.,
									      0,
									      -999.,
									      -999.) ) );
	 }
       ++iphot;
     }
	 
   //   if (_debug)
   //   std::cout << "Preselected " << photonsForVertexComputation.size() << " photons " << std::endl;
   std::map<int,PhotonInfo>::const_iterator iphot1;
   std::map<int,PhotonInfo>::const_iterator iphot2;


   
   GlobalVertexInfo aVtx;
   aVtx.fillInfo(iEvent,iSetup);     
   nPreselPhotonPairs=0;
   vtxAna->clear();

   if (photonsForVertexComputation.size()>=2)
     {
       for (iphot1=photonsForVertexComputation.begin();iphot1!=(--photonsForVertexComputation.end());++iphot1)
	 {
	   iphot2=iphot1;
	   ++iphot2;
	   for (;iphot2!=photonsForVertexComputation.end();++iphot2)
	     {
	       if (nPreselPhotonPairs>19)
		 {
		   std::cout << "MAXIMUM NUMBER OF ALLOWED PHOTON PAIRS EXCEEDED!!" << std::endl;
		   continue;
	     }
	       
	       indexPreselPhot1[nPreselPhotonPairs]=(*iphot1).first;
	       indexPreselPhot2[nPreselPhotonPairs]=(*iphot2).first;
	       
	       //	       std::cout << "+++ " << indexPreselPhot1[nPreselPhotonPairs] << "," << indexPreselPhot2[nPreselPhotonPairs] << std::endl;
	       //Vertex Analysis	   						 
	       //    for(int ii=0; ii<nvertex; ++ii) 
	       //      {
	       //        vrank[ii]=-1;
	       //        vptbal[ii]=-9999.;
	       //        vptasym[ii]=-9999.;
	       //        vlogsumpt2[ii]=-9999.;
	       //      }
	       
	       
	       //   //Prepare vertex information
	       //    if (photonsForVertexComputation.size()>=2)
	       //     {
	       

	       //Run HggVertexAnalysis
	       vtxAna->analyze(aVtx,(*iphot1).second,(*iphot2).second);

	       int p1 = (*iphot1).second.id(), p2 = (*iphot2).second.id();
	       // assert( p1 == vtxAna.pho1() && p2 == vtxAna.pho2() );
	       vtxAna->setPairID(p1,p2);
	       /*

	       
	       //Choiche of the vertex
	       // preselect vertices : all vertices

	       std::vector<int> preselAll;
	       for(int i=0; i<aVtx.nvtx() ; i++) {
		 preselAll.push_back(i); 
	       }
	       
	       float zconv = 0; 
	       float dzconv = 0;
	       bool isConv = false;
	       
	       if ( ((*iphot1).second.isAConversion() || (*iphot2).second.isAConversion() ) )  {
		 isConv=true;
		 
		 if ((*iphot1).second.isAConversion()  && !(*iphot2).second.isAConversion() ){
		   zconv  = vtxAnaFromConv->vtxZ((*iphot1).second);
		   dzconv = vtxAnaFromConv->vtxdZ((*iphot1).second);
		 }
		 
		 if ((*iphot2).second.isAConversion() && !(*iphot1).second.isAConversion()){
		   zconv  = vtxAnaFromConv->vtxZ((*iphot2).second);
		   dzconv = vtxAnaFromConv->vtxdZ((*iphot2).second);
		 }
		 
		 if ( (*iphot1).second.isAConversion() && (*iphot2).second.isAConversion()){
		   float z1  = vtxAnaFromConv->vtxZ((*iphot1).second);
		   float dz1 = vtxAnaFromConv->vtxdZ((*iphot1).second);
		   
		   float z2  = vtxAnaFromConv->vtxZ((*iphot2).second);
		   float dz2 = vtxAnaFromConv->vtxdZ((*iphot2).second);
		   zconv  = (z1/dz1/dz1 + z2/dz2/dz2)/(1./dz1/dz1 + 1./dz2/dz2 );  // weighted average
		   dzconv = sqrt( 1./(1./dz1/dz1 + 1./dz2/dz2)) ;		

		 }
		 
		 // preselect vertices : only vertices in a window zconv +/- dzconv
		 for(int i=0; i < aVtx.nvtx(); i++) {
		 // 	     //	    std::cout << "zconv " << zconv << " fabs(zconv - aVtx.vtxz(i) ) " << fabs(zconv - aVtx.vtxz(i) ) << " dzconv " << dzconv << std::endl;
		   if ( fabs(zconv - aVtx.vtxz(i) ) < dzconv ) 
		     preselConv.push_back(i); 
		 }
		 
	       }

	       // preselection 
	       // 	if ( preselConv.size()==0 ) 
	       //           vtxAna->preselection(preselAll);
	       //         else 
	       // 	  {
	       // 	    //	    std::cout << "&&&&&& Using conversions " << preselConv.size()  << " vertices instead of " <<  preselAll.size() << std::endl;
	       // 	    vtxAna->preselection(preselConv);
	       // 	  }
	       
	       //  	rankprod = vtxAna->rankprod(rankVariables);
	       
	       */

	       std::vector<int> rankprodAll = vtxAna->rank(*tmvaPerVtxReader_,tmvaPerVtxMethod); 
	       float vtxEvtMva = vtxAna->perEventMva( *tmvaPerEvtReader_, tmvaPerEvtMethod, rankprodAll );
	       float vtxProbability= vtxAna->vertexProbability(vtxEvtMva, nvertex); 
	       //	       std::cout << "*** " << vtxEvtMva << " , " << nPreselPhotonPairs << std::endl;
//  	       std::vector<int> rankprodAll = vtxAna->rankprod(rankVariables);
	       
	       
// 	       int iClosestConv = -1;
// 	       float dminconv = 9999999;
	       
// 	       //Using here pt wrt sumpt2 (default) to assess diPhotonPt
// 	       TLorentzVector dipho =  (*iphot1).second.p4(vx[0],vy[0],vz[0]) +  (*iphot2).second.p4(vx[0],vy[0],vz[0]);
	       
// 	       unsigned int nbest ;
// 	       if (  dipho.Pt() < 30 ) nbest = 5;
// 	       else nbest = 3;
// 	       if (rankprodAll.size() < nbest ) nbest = rankprodAll.size();
	       
// 	       for (unsigned int ii = 0; ii < nbest; ii++ ){
// 		 TVector3 vtxpos(vx[rankprodAll[ii]],vy[rankprodAll[ii]],vz[rankprodAll[ii]]);
// 		 if ( isConv && fabs( vtxpos.Z()-zconv ) < dzconv && fabs(vtxpos.Z() - zconv ) < dminconv){
// 		   iClosestConv = rankprodAll[ii];
// 		   dminconv = fabs(vtxpos.Z()-zconv );
// 		 }
// 	       }
	       
	 
// 	       std::vector<int> rankprod;
// 	       rankprod.clear();
 	       int bestVtx=-1;
// 	       if (iClosestConv!=-1 ) 
// 		 bestVtx = iClosestConv;
// 	       else
	       bestVtx = rankprodAll[0];

	       //	       if (_debug)
	       //	       cout << "best vertex for couple " << nPreselPhotonPairs << " iphot1 " << indexPreselPhot1[nPreselPhotonPairs] << " iphot2 " << indexPreselPhot2[nPreselPhotonPairs] 
	       //   << " is " << bestVtx << " ptbal " << vtxAna->ptbal(bestVtx) << " ptasymm " << vtxAna->ptasym(bestVtx) <<  std::endl;

	       //Filling 
	       vrankPhotonPairs[nPreselPhotonPairs]=bestVtx;
	       vevtMvaPhotonPairs[nPreselPhotonPairs]=vtxEvtMva;
	       vevtProbPhotonPairs[nPreselPhotonPairs]=vtxProbability;
	       vptbalPhotonPairs[nPreselPhotonPairs]=vtxAna->ptbal(bestVtx);
	       vptasymPhotonPairs[nPreselPhotonPairs]=vtxAna->ptasym(bestVtx);
	       ++nPreselPhotonPairs;
	     }
	 }
     }


   
   // chiara
   unsigned index_gsf = 0;
   for (GsfElectronCollection::const_iterator itElectron = ElectronHandle->begin();
	itElectron != ElectronHandle->end(); ++itElectron, ++index_gsf) {

     // chiara
     GsfElectronRef eleRef = GsfElectronRef(ElectronHandle, index_gsf);
     
     const EBRecHitCollection* rechits = ( itElectron->isEB()) ? rhitseb : rhitsee;     
     const Ptr<CaloCluster> theSeed = itElectron->superCluster()->seed();
     float e9=EcalClusterTools::e3x3( *theSeed, &(*rechits), topology );
     float sigIPhiIPhi=sqrt((EcalClusterTools::localCovariances( *theSeed, &(*rechits), topology ))[2]);
     electron_pt[nEle]          = itElectron->pt();
     electron_energy[nEle]      = itElectron->energy();
     electron_ecalEnergy[nEle]  = itElectron->ecalEnergy();                 // chiara
     electron_trackPatVtx[nEle] = itElectron->trackMomentumAtVtx().R();     // chiara
     electron_px[nEle]       = itElectron->px();
     electron_py[nEle]       = itElectron->py();
     electron_pz[nEle]       = itElectron->pz();
     electron_vx[nEle]       = itElectron->vx();
     electron_vy[nEle]       = itElectron->vy();
     electron_vz[nEle]       = itElectron->vz();
     electron_phi[nEle]      = itElectron->phi();
     electron_eta[nEle]      = itElectron->eta();
     electron_charge[nEle]       = itElectron->charge();
     electron_fBrem[nEle]       = itElectron->fbrem();
     electron_EoP[nEle]       = itElectron->eSuperClusterOverP();
     electron_OneOverEMinusOneOverP[nEle]       = (1/itElectron->caloEnergy())-(1/itElectron->trackMomentumAtVtx().R());
     electron_r9[nEle]       = e9/itElectron->superCluster()->rawEnergy();
     electron_misHits[nEle]       = itElectron->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
     ConversionInfo convInfo = convFinder.getConversionInfo(*itElectron,tracks, bfield);      
     electron_dist[nEle] = (convInfo.dist() == -9999.? 9999:convInfo.dist());
     electron_dcot[nEle] = (convInfo.dcot() == -9999.? 9999:convInfo.dcot());
     bool matchesConv = ConversionTools::hasMatchedConversion(*eleRef,hConversions,recoBeamSpotHandle->position());
     electron_matchedConv[nEle] = matchesConv;
     electron_seedType[nEle]       = itElectron->ecalDrivenSeed()+2*itElectron->trackerDrivenSeed();
     electron_nSubClusters[nEle]       = itElectron->numberOfBrems();
     //     electron_HoE[nEle]          = itElectron->hadronicOverEm();
     electron_HoE[nEle]          = itElectron->hcalOverEcalBc(); //new H/E 2012
     electron_pFlowMVA[nEle]          = itElectron->mvaOutput().mva;
     electron_SigmaIetaIeta[nEle]= itElectron->sigmaIetaIeta();
     electron_SigmaIphiIphi[nEle]= sigIPhiIPhi;
     electron_trkIso[nEle]       = itElectron->dr04TkSumPt() ;
     electron_ecalIso[nEle]      = itElectron->dr04EcalRecHitSumEt();
     //electron_hcalIso[nEle]      = itElectron->dr04HcalTowerSumEt();
     electron_hcalIso[nEle]      = itElectron->dr04HcalDepth1TowerSumEtBc(); // new H/E 2012
     electron_trkIso03[nEle]       = itElectron->dr03TkSumPt() ;
     electron_ecalIso03[nEle]      = itElectron->dr03EcalRecHitSumEt();
     //electron_hcalIso03[nEle]      = itElectron->dr03HcalTowerSumEt();
     electron_hcalIso03[nEle]      = itElectron->dr03HcalDepth1TowerSumEtBc(); // new H/E 2012
     electron_dEtaIn[nEle]       = itElectron->deltaEtaSuperClusterTrackAtVtx();
     electron_dPhiIn[nEle]       = itElectron->deltaPhiSuperClusterTrackAtVtx();
     electron_sc_energy[nEle]    = itElectron->superCluster()->energy();
     electron_sc_eta[nEle]       = itElectron->superCluster()->eta();
     electron_sc_phi[nEle]       = itElectron->superCluster()->phi();
    
     // chiara ----------------------------
     if (IS2012) {
       const reco::GsfElectron *ele = &(*itElectron);
       
       double mvaidnontrig = myMVANonTrig->mvaValue(*ele, 
						    VertexHandleJetId->front(), 
						    thebuilder,
						    lazyTools,
						    false);
       
       double mvaidtrig = myMVATrig->mvaValue(*ele,
					      VertexHandleJetId->front(),
					      thebuilder,
					      lazyTools,
					      false);
       
       electron_mvaNonTrig[nEle] = mvaidnontrig;
       electron_mvaTrig[nEle]    = mvaidtrig;
       
       const isoFromPFCandsMap & electronsPfCandChHad03IsoVal  = *( (*eIsoFromPFCandsValueMap_)[0] );
       const isoFromPFCandsMap & electronsPfCandNHad03IsoVal   = *( (*eIsoFromPFCandsValueMap_)[1] );
       const isoFromPFCandsMap & electronsPfCandPhoton03IsoVal = *( (*eIsoFromPFCandsValueMap_)[2] );
       const isoFromPFCandsMap & electronsPfCandChHad04IsoVal  = *( (*eIsoFromPFCandsValueMap_)[3] );
       const isoFromPFCandsMap & electronsPfCandNHad04IsoVal   = *( (*eIsoFromPFCandsValueMap_)[4] );
       const isoFromPFCandsMap & electronsPfCandPhoton04IsoVal = *( (*eIsoFromPFCandsValueMap_)[5] );
       const isoFromPFCandsMap & electronsPfCandChHad05IsoVal  = *( (*eIsoFromPFCandsValueMap_)[6] );
       const isoFromPFCandsMap & electronsPfCandNHad05IsoVal   = *( (*eIsoFromPFCandsValueMap_)[7] );
       const isoFromPFCandsMap & electronsPfCandPhoton05IsoVal = *( (*eIsoFromPFCandsValueMap_)[8] );
       electron_chHad03Iso[nEle] = electronsPfCandChHad03IsoVal[eleRef]; 
       electron_nHad03Iso[nEle]  = electronsPfCandNHad03IsoVal[eleRef];
       electron_phot03Iso[nEle]  = electronsPfCandPhoton03IsoVal[eleRef];
       electron_chHad04Iso[nEle] = electronsPfCandChHad04IsoVal[eleRef]; 
       electron_nHad04Iso[nEle]  = electronsPfCandNHad04IsoVal[eleRef];
       electron_phot04Iso[nEle]  = electronsPfCandPhoton04IsoVal[eleRef];
       electron_chHad05Iso[nEle] = electronsPfCandChHad05IsoVal[eleRef]; 
       electron_nHad05Iso[nEle]  = electronsPfCandNHad05IsoVal[eleRef];
       electron_phot05Iso[nEle]  = electronsPfCandPhoton05IsoVal[eleRef];
     } else {
       electron_mvaNonTrig[nEle] = -999.;
       electron_mvaTrig[nEle]    = -999.;
       electron_chHad03Iso[nEle] = -999.;
       electron_nHad03Iso[nEle]  = -999.; 
       electron_phot03Iso[nEle]  = -999.; 
       electron_chHad04Iso[nEle] = -999.; 
       electron_nHad04Iso[nEle]  = -999.; 
       electron_phot04Iso[nEle]  = -999.; 
       electron_chHad05Iso[nEle] = -999.; 
       electron_nHad05Iso[nEle]  = -999.;
       electron_phot05Iso[nEle]  = -999.;
     }
     // chiara ----------------------------

     nEle++;
   }//end of for loop


  if (dumpKT4Jets_)
    {
      for (CaloJetCollection::const_iterator it = jetskt4->begin(); 
	   it != jetskt4->end(); ++it) {
	
	if (nJet_kt4>=100) {cout << "number of reco jets kt 04 is larger than 100. Skipping" << endl; continue;}
	if (it->pt() > calojetptthr_) {
	  
	  ptJet_kt4[nJet_kt4] = it->pt();
	  eJet_kt4[nJet_kt4] = it->energy();	 
	  etaJet_kt4[nJet_kt4] = it->eta();	 
	  phiJet_kt4[nJet_kt4] = it->phi();
	  //emfJet_kt4[nJet_kt4] = it->emEnergyFraction();	      
	  emfJet_kt4[nJet_kt4] = fixEMF(it->emEnergyFraction(), it->eta());
	  
	  nJet_kt4++;
	}
      }
    }
   
  if (dumpKT6Jets_)
    {
      for (CaloJetCollection::const_iterator it = jetskt6->begin(); 
	   it != jetskt6->end(); ++it) {
	
	if (nJet_kt6>=100) {cout << "number of reco jets kt 06 is larger than 100. Skipping" << endl; continue;}
	if (it->pt() > calojetptthr_) {
	  
	  ptJet_kt6[nJet_kt6] = it->pt(); 
	  eJet_kt6[nJet_kt6] = it->energy();	 
	  etaJet_kt6[nJet_kt6] = it->eta();	 
	  phiJet_kt6[nJet_kt6] = it->phi();	      
	  //emfJet_kt6[nJet_kt6] = it->emEnergyFraction();     
	  emfJet_kt6[nJet_kt6] = fixEMF(it->emEnergyFraction(), it->eta());
	  
	  nJet_kt6++;
	}
      }
    }

  if (dumpAKT5Jets_)
    {
      for (CaloJetCollection::const_iterator it = jetsakt5->begin(); 
	   it != jetsakt5->end(); ++it) {
	
	if (nJet_akt5>=100) {cout << "number of reco jets akt 05 is larger than 100. Skipping" << endl; continue;}
	if (it->pt() > calojetptthr_) {
	  //        jetID_->calculate(iEvent, *it);
	  
	  ptJet_akt5[nJet_akt5] = it->pt();	 
	  
	  // Jet Energy Scale Corrections on-the-fly     
	  CaloJet  correctedJet = *it;
	  edm::RefToBase<reco::Jet> jetRef(edm::Ref<CaloJetCollection>(jetsakt5,nJet_akt5));
	  double scale = corrector_akt5->correction(*it,iEvent,iSetup);
	  correctedJet.scaleEnergy(scale);
	  ptCorrJet_akt5[nJet_akt5] = correctedJet.pt();
	  
	  eJet_akt5[nJet_akt5] = it->energy();	 
	  etaJet_akt5[nJet_akt5] = it->eta();	 
	  phiJet_akt5[nJet_akt5] = it->phi();	      
	  //emfJet_akt5[nJet_akt5] = it->emEnergyFraction();
	  emfJet_akt5[nJet_akt5] = fixEMF(it->emEnergyFraction(), it->eta());
	  //        n90Jet_akt5[nJet_akt5] = jetID_->hitsInN90();	      
	  //        n90HitsJet_akt5[nJet_akt5] = jetID_->n90Hits();	      
	  //        fHPDJet_akt5[nJet_akt5] = jetID_->fHPD();	      
	  //        fRBXJet_akt5[nJet_akt5] = jetID_->fRBX();	      
	  
	  nJet_akt5++;
	}
      }
    }

  if (dumpAKT7Jets_)
    {
      for (CaloJetCollection::const_iterator it = jetsakt7->begin(); 
	   it != jetsakt7->end(); ++it) {
	
	if (nJet_akt7>=100) {cout << "number of reco jets akt 07 is larger than 100. Skipping" << endl; continue;}
	if (it->pt() > calojetptthr_) {
	  //       jetID_->calculate(iEvent, *it);
	  ptJet_akt7[nJet_akt7] = it->pt();	 
	  
	  // Jet Energy Scale Corrections on-the-fly     
	  CaloJet  correctedJet = *it;
	  edm::RefToBase<reco::Jet> jetRef(edm::Ref<CaloJetCollection>(jetsakt7,nJet_akt7));
	  double scale = corrector_akt5->correction(*it,iEvent,iSetup);
	  correctedJet.scaleEnergy(scale);
	  ptCorrJet_akt7[nJet_akt7] = correctedJet.pt();
	  
	  eJet_akt7[nJet_akt7] = it->energy();	 
	  etaJet_akt7[nJet_akt7] = it->eta();	 
	  phiJet_akt7[nJet_akt7] = it->phi();	      
	  //emfJet_akt7[nJet_akt7] = it->emEnergyFraction();
	  emfJet_akt7[nJet_akt7] = fixEMF(it->emEnergyFraction(), it->eta());     
	  //        n90Jet_akt7[nJet_akt7] = jetID_->hitsInN90();	      
	  //        n90HitsJet_akt7[nJet_akt7] = jetID_->n90Hits();	      
	  //        fHPDJet_akt7[nJet_akt7] = jetID_->fHPD();	      
	  //        fRBXJet_akt7[nJet_akt7] = jetID_->fRBX();	      
	  
	  nJet_akt7++;
	}
      }
    }
   
  if (dumpJPTAKT5Jets_)
    {
      for (JPTJetCollection::const_iterator it = jptjetsak5->begin(); 
	   it != jptjetsak5->end(); ++it) {
	
	if (nJet_jptak5>=100) {cout << "number of reco jets jptakt 05 is larger than 100. Skipping" << endl; continue;}
	if (nJet_jptak5 < jptjetnmin_ || it->pt() > jptjetptthr_) {
	  
	  // Jet Energy Scale Corrections on-the-fly     
	  JPTJet  correctedJet = *it;
	  edm::RefToBase<reco::Jet> jetRef(edm::Ref<JPTJetCollection>(jptjetsak5,nJet_jptak5));
	  double scale = corrector_jptak5->correction(*it,iEvent,iSetup);
	  correctedJet.scaleEnergy(scale);
	  ptCorrJet_jptak5[nJet_jptak5] = correctedJet.pt();
	  
	  ptJet_jptak5[nJet_jptak5] = it->pt();	 
	  eJet_jptak5[nJet_jptak5] = it->energy();	 
	  etaJet_jptak5[nJet_jptak5] = it->eta();	 
	  phiJet_jptak5[nJet_jptak5] = it->phi();	      
	  //emfJet_jptak5[nJet_jptak5] = fixEMF(it->emEnergyFraction(), it->eta());
	  
	  nJet_jptak5++;
	}
      }
    }

  if (dumpKT4Jets_)
    {
      for (PFJetCollection::const_iterator it = pfjetskt4->begin(); 
	   it != pfjetskt4->end(); ++it) {
	
	if (nJet_pfkt4>=100) {cout << "number of reco jets pfkt4 is larger than 100. Skipping" << endl; continue;}
	if (nJet_pfkt4 < pfjetnmin_ || it->pt() > pfjetptthr_) {
	  
	  ptJet_pfkt4[nJet_pfkt4] = it->pt();
	  eJet_pfkt4[nJet_pfkt4] = it->energy();	 
	  etaJet_pfkt4[nJet_pfkt4] = it->eta();	 
	  phiJet_pfkt4[nJet_pfkt4] = it->phi();	      
	  
	  nJet_pfkt4++;
	}
      }
    }
  Int_t nPfCand_all = 0;
  Int_t nPfCand_injet = 0;
  Int_t nChargedHadrons = 0;
  Int_t nChargedHadronsgoodvtx = 0;
  Int_t nChargedHadronsnoothervtx = 0;
  Int_t nPhotons = 0;
  Int_t nNeutralHadrons = 0;
  Int_t nElectrons = 0;
  Int_t nMuons = 0;
  Int_t nHFHadrons = 0;
  Int_t nHFEM = 0;
  
  TLorentzVector p4PfCand_all;
  TLorentzVector p4ChargedHadrons_uncl;
  TLorentzVector p4ChargedHadronsgoodvtx_uncl;
  TLorentzVector p4ChargedHadronsnoothervtx_uncl;
  TLorentzVector p4Photons_uncl;
  TLorentzVector p4NeutralHadrons_uncl;
  TLorentzVector p4Electrons_uncl;
  TLorentzVector p4Muons_uncl;
  TLorentzVector p4HFHadrons_uncl;
  TLorentzVector p4HFEM_uncl;

  sumptpfcand_all = 0;
  sumptChargedHadrons_uncl = 0;
  sumptChargedHadronsgoodvtx_uncl = 0;
  sumptChargedHadronsnoothervtx_uncl = 0;
  sumptPhotons_uncl = 0;
  sumptNeutralHadrons_uncl = 0;
  sumptElectrons_uncl = 0;
  sumptMuons_uncl = 0;
  sumptHFHadrons_uncl = 0;
  sumptHFEM_uncl = 0;
	  
  for (PFCandidateCollection::const_iterator jt = PFCandidates->begin();
       jt != PFCandidates->end(); ++jt) {

    PFCandidate::ParticleType id = jt->particleId();
    // Convert particle momentum to normal TLorentzVector, wrong type :(
    math::XYZTLorentzVectorD const& p4t = jt->p4();
    TLorentzVector p4(p4t.px(), p4t.py(), p4t.pz(), p4t.energy());
    TLorentzVector jetp4;
    nPfCand_all += 1;
    p4PfCand_all += p4;
    sumptpfcand_all += p4.Pt();
 
    bool injet(0);
    int countjet(0);

    for (PFJetCollection::const_iterator it = pfjetsakt5->begin(); 
	 it != pfjetsakt5->end(); ++it) {
      
      vector<PFCandidatePtr> pfCandidates = it->getPFConstituents();
      
      for (vector<PFCandidatePtr>::const_iterator jtjet = pfCandidates.begin();
	   jtjet != pfCandidates.end(); ++jtjet) {
	

	if (jt->pt() == (*jtjet)->pt() && it->pt() > pfjetptthr_ && countjet<100) 
	  {
	    injet = 1;
	    nPfCand_injet += 1;    
	  }
      }     
      countjet++;
      
    }

    if(!injet) {
      PFCandidate::ParticleType id = jt->particleId();
      
      if (id==PFCandidate::h) { // charged hadrons
	nChargedHadrons += 1;
	p4ChargedHadrons_uncl += p4;
	sumptChargedHadrons_uncl += p4.Pt();
	int myVertex=-1;
	//---- loop over all vertices ----------------------------
	for(unsigned ivtx = 0;ivtx <  VertexHandle->size();ivtx++) {
	  //---- loop over the tracks associated with the vertex ---
	  if (!((* VertexHandle)[ivtx].isFake()) && (* VertexHandle)[ivtx].ndof() >= DEF_GOODVTX_NDOF && fabs((* VertexHandle)[ivtx].z()) <= DEF_GOODVTX_Z) {
	    for(reco::Vertex::trackRef_iterator i_vtxTrk = (* VertexHandle)[ivtx].tracks_begin(); i_vtxTrk != (* VertexHandle)[ivtx].tracks_end(); ++i_vtxTrk) {
	      //---- match the jet track to the track from the vertex ----
	      reco::TrackRef trkRef(i_vtxTrk->castTo<reco::TrackRef>());
	      //---- check if the tracks match -------------------------
	      if (trkRef == jt->trackRef()) 
		{
		  myVertex=ivtx;
		  if(ivtx == 0){
		    nChargedHadronsgoodvtx += 1;
		    p4ChargedHadronsgoodvtx_uncl += p4;
		    sumptChargedHadronsgoodvtx_uncl += p4.Pt();
		    nChargedHadronsnoothervtx += 1;
		    p4ChargedHadronsnoothervtx_uncl += p4;
		    sumptChargedHadronsnoothervtx_uncl += p4.Pt();
		  }
		  break;
		}
	    }
	  }	  
	}
	if (myVertex==-1){
	  nChargedHadronsnoothervtx += 1;
	  p4ChargedHadronsnoothervtx_uncl += p4;	  
	  sumptChargedHadronsnoothervtx_uncl += p4.Pt();
	}
	
      }
      else if (id==PFCandidate::e) { // electrons
	nElectrons += 1;
	p4Electrons_uncl += p4;
	sumptElectrons_uncl += p4.Pt();
      }
      else if (id==PFCandidate::mu) { // muons
	nMuons += 1;
	p4Muons_uncl += p4;
 	sumptMuons_uncl += p4.Pt();
     }
      else if (id==PFCandidate::gamma) { // photons
	nPhotons += 1;
	p4Photons_uncl += p4;
	sumptPhotons_uncl += p4.Pt();
      }
      else if (id==PFCandidate::h0) { // neutral hadrons
	nNeutralHadrons += 1;
	p4NeutralHadrons_uncl += p4;
	sumptNeutralHadrons_uncl += p4.Pt();
      }
      else if (id==PFCandidate::h_HF) { // HF hadrons
	nHFHadrons += 1;
	p4HFHadrons_uncl += p4;
	sumptHFHadrons_uncl += p4.Pt();
      }
      else if (id==PFCandidate::egamma_HF) { // HF EM clusters
	nHFEM += 1;
	p4HFEM_uncl += p4;
	sumptHFEM_uncl += p4.Pt();
      }
      else cout << "SCREAMMMMMM" << endl;
    }
  }
  
  npfcand_all =  nPfCand_all;
  const TLorentzVector *p = 0;
  p = &p4PfCand_all;
  epfcand_all = p->E();
  if(p->E()){
    ptpfcand_all = p->Pt();
    etapfcand_all = p->Eta();
    phipfcand_all = p->Phi();
  }else{	    
    ptpfcand_all = 0.;
    etapfcand_all = -999.;
    phipfcand_all = -999.;
  }	

  nChargedHadrons_uncl =  nChargedHadrons;
  p = &p4ChargedHadrons_uncl;
  eChargedHadrons_uncl = p->E();
  if(p->E()){
    ptChargedHadrons_uncl = p->Pt();
    etaChargedHadrons_uncl = p->Eta();
    phiChargedHadrons_uncl = p->Phi();
  }else{	    
    ptChargedHadrons_uncl = 0.;
    etaChargedHadrons_uncl = -999.;
    phiChargedHadrons_uncl = -999.;
  }	
  
  nChargedHadronsgoodvtx_uncl =  nChargedHadronsgoodvtx;
  p = &p4ChargedHadronsgoodvtx_uncl;
  eChargedHadronsgoodvtx_uncl = p->E();
  if(p->E()){
    ptChargedHadronsgoodvtx_uncl = p->Pt();
    etaChargedHadronsgoodvtx_uncl = p->Eta();
    phiChargedHadronsgoodvtx_uncl = p->Phi();
  }else{	    
    ptChargedHadronsgoodvtx_uncl = 0.;
    etaChargedHadronsgoodvtx_uncl = -999.;
    phiChargedHadronsgoodvtx_uncl = -999.;
  }	
  
  nChargedHadronsnoothervtx_uncl =  nChargedHadronsnoothervtx;
  p = &p4ChargedHadronsnoothervtx_uncl;
  eChargedHadronsnoothervtx_uncl = p->E();
  if(p->E()){
    ptChargedHadronsnoothervtx_uncl = p->Pt();
    etaChargedHadronsnoothervtx_uncl = p->Eta();
    phiChargedHadronsnoothervtx_uncl = p->Phi();
  }else{	    
    ptChargedHadronsnoothervtx_uncl = 0.;
    etaChargedHadronsnoothervtx_uncl = -999.;
    phiChargedHadronsnoothervtx_uncl = -999.;
  }	
  
  nElectrons_uncl =  nElectrons;
  p = &p4Electrons_uncl;
  eElectrons_uncl = p->E();
  if(p->E()){
    ptElectrons_uncl = p->Pt();
    etaElectrons_uncl = p->Eta();
    phiElectrons_uncl = p->Phi();
  }else{
    ptElectrons_uncl = 0.;
    etaElectrons_uncl = -999.;
    phiElectrons_uncl = -999.;
  }
  
  nMuons_uncl =  nMuons;
  p = &p4Muons_uncl;
  eMuons_uncl = p->E();
  if(p->E()){
    ptMuons_uncl = p->Pt();
    etaMuons_uncl = p->Eta();
    phiMuons_uncl = p->Phi();
  }else{
    ptMuons_uncl = 0.;
    etaMuons_uncl = -999.;
    phiMuons_uncl = -999.;
  }
  
  nPhotons_uncl =  nPhotons;
  p = &p4Photons_uncl;
  ePhotons_uncl = p->E();
  if(p->E()){	  
    ptPhotons_uncl = p->Pt();
    etaPhotons_uncl = p->Eta();
    phiPhotons_uncl = p->Phi();
  }else{
    ptPhotons_uncl = 0.;
    etaPhotons_uncl = -999.;
    phiPhotons_uncl = -999.;
  }
  
  nNeutralHadrons_uncl =  nNeutralHadrons;
  p = &p4NeutralHadrons_uncl;
  eNeutralHadrons_uncl = p->E();
  if(p->E()){	  
    ptNeutralHadrons_uncl = p->Pt();
    etaNeutralHadrons_uncl = p->Eta();
    phiNeutralHadrons_uncl = p->Phi();
  }else{
    ptNeutralHadrons_uncl = 0.;
    etaNeutralHadrons_uncl = -999.;
    phiNeutralHadrons_uncl = -999.;
  }
  
  nHFHadrons_uncl =  nHFHadrons;
  p = &p4HFHadrons_uncl;
  eHFHadrons_uncl = p->E();
  if(p->E()){	  
    ptHFHadrons_uncl = p->Pt();
    etaHFHadrons_uncl = p->Eta();
    phiHFHadrons_uncl = p->Phi();
  }else{
    ptHFHadrons_uncl = 0.;
    etaHFHadrons_uncl = -999.;
    phiHFHadrons_uncl = -999.;
  }	    
  
  nHFEM_uncl =  nHFEM;
  p = &p4HFEM_uncl;
  eHFEM_uncl = p->E();
  if(p->E()){	  
    ptHFEM_uncl = p->Pt();
    etaHFEM_uncl = p->Eta();
    phiHFEM_uncl = p->Phi();
  }else{
    ptHFEM_uncl = 0.;
    etaHFEM_uncl = -999.;
    phiHFEM_uncl = -999.;
  }	    

  Int_t nCandinjet = 0;

  if (dumpAKT5Jets_)
    {
      for (PFJetCollection::const_iterator it = pfjetsakt5->begin(); 
	   it != pfjetsakt5->end(); ++it) {
	
	if (nJet_pfakt5>=100) {cout << "number of reco jets pfakt5 is larger than 100. Skipping" << endl; continue;}
	if (nJet_pfakt5 < pfjetnmin_ || it->pt() > pfjetptthr_) {
	  
	  ptJet_pfakt5[nJet_pfakt5] = it->pt();
	  eJet_pfakt5[nJet_pfakt5] = it->energy();	 
	  etaJet_pfakt5[nJet_pfakt5] = it->eta();	 
	  phiJet_pfakt5[nJet_pfakt5] = it->phi();	      

	  // Jet Energy Scale Corrections on-the-fly     
	  PFJet  correctedJet = *it;
	  edm::RefToBase<reco::Jet> jetRef(edm::Ref<PFJetCollection>(pfjetsakt5,nJet_pfakt5));
	  double scale = corrector_pfakt5->correction(*it,iEvent,iSetup);
	  correctedJet.scaleEnergy(scale);
	  ptCorrJet_pfakt5[nJet_pfakt5] = correctedJet.pt();
	  
	  // Extra variables for PFlow studies
	  Int_t nChargedHadrons = 0;
	  Int_t nChargedHadronsgoodvtx = 0;
	  Int_t nChargedHadronsnoothervtx = 0;
	  Int_t nPhotons = 0;
	  Int_t nNeutralHadrons = 0;
	  Int_t nElectrons = 0;
	  Int_t nMuons = 0;
	  Int_t nHFHadrons = 0;
	  Int_t nHFEM = 0;
	  
	  TLorentzVector p4ChargedHadrons;
	  TLorentzVector p4ChargedHadronsgoodvtx;
	  TLorentzVector p4ChargedHadronsnoothervtx;
	  TLorentzVector p4Photons;
	  TLorentzVector p4NeutralHadrons;
	  TLorentzVector p4Electrons;
	  TLorentzVector p4Muons;
	  TLorentzVector p4HFHadrons;
	  TLorentzVector p4HFEM;
	  
	  sumptChargedHadrons_pfakt5[nJet_pfakt5] = 0;
	  sumptChargedHadronsgoodvtx_pfakt5[nJet_pfakt5] = 0;
	  sumptChargedHadronsnoothervtx_pfakt5[nJet_pfakt5] = 0;
	  sumptPhotons_pfakt5[nJet_pfakt5] = 0;
	  sumptNeutralHadrons_pfakt5[nJet_pfakt5] = 0;
	  sumptElectrons_pfakt5[nJet_pfakt5] = 0;
	  sumptMuons_pfakt5[nJet_pfakt5] = 0;
	  sumptHFHadrons_pfakt5[nJet_pfakt5] = 0;
	  sumptHFEM_pfakt5[nJet_pfakt5] = 0;

	  vector<PFCandidatePtr> pfCandidates = it->getPFConstituents();
	  
	  float rms_cands_wrong=0.;

	  float SumW=0;
	  float SumW2=0;
	  float SumDeta=0;
	  float SumDeta2=0;
	  float SumDphi=0;
	  float SumDphi2=0;
	  float SumDetaDphi=0;
	  
	  float SumW_QC=0;
	  float SumW2_QC=0;
	  float SumDeta_QC=0;
	  float SumDeta2_QC=0;
	  float SumDphi_QC=0;
	  float SumDphi2_QC=0;
	  float SumDetaDphi_QC=0;
	  
	  float Eta0=it->eta();
	  float Phi0=it->phi();

        // initialize:
        rmsCandJet_pfakt5[nJet_pfakt5] =  -999.;
        ptDJet_pfakt5[nJet_pfakt5] =      -999.;
        axis1Jet_pfakt5[nJet_pfakt5] =    -999.;
        axis2Jet_pfakt5[nJet_pfakt5] =    -999.;
        pullJet_pfakt5[nJet_pfakt5] =     -999.;
        tanaJet_pfakt5[nJet_pfakt5]  =    -999.;

        rmsCandTrue_QCJet_pfakt5[nJet_pfakt5] =  -999.;
        ptD_QCJet_pfakt5[nJet_pfakt5] =      -999.;
        axis1_QCJet_pfakt5[nJet_pfakt5] =    -999.;
        axis2_QCJet_pfakt5[nJet_pfakt5] =    -999.;
        pull_QCJet_pfakt5[nJet_pfakt5] =     -999.;
        tana_QCJet_pfakt5[nJet_pfakt5]  =    -999.;

        RchgJet_pfakt5[nJet_pfakt5] = 0.;
        RneutralJet_pfakt5[nJet_pfakt5] = 0.;
        RJet_pfakt5[nJet_pfakt5] = 0.;
        Rchg_QCJet_pfakt5[nJet_pfakt5] = 0.;

	  pTMaxJet_pfakt5[nJet_pfakt5] = 0.;
	  pTMaxChgJet_pfakt5[nJet_pfakt5] = 0.;
	  pTMaxNeutralJet_pfakt5[nJet_pfakt5] = 0.;
        pTMaxChg_QCJet_pfakt5[nJet_pfakt5] = 0.;

        nChg_ptCutJet_pfakt5[nJet_pfakt5] = 0;
        nChg_QCJet_pfakt5[nJet_pfakt5] = 0;
        nChg_ptCut_QCJet_pfakt5[nJet_pfakt5] = 0;
        nNeutral_ptCutJet_pfakt5[nJet_pfakt5] = 0;

        std::vector<bool> jetPart_forMult,jetPart_forAxis;

	  
	  for(int i=0;i<it->nConstituents();++i)
		{

		reco::TrackRef itrk ;
		reco::PFCandidatePtr  part = it->getPFConstituent(i);

		double pt=part->pt();
		double eta=part->eta();
		double phi=part->phi();

		PFCandidate::ParticleType id = part->particleId();
		// Convert particle momentum to normal TLorentzVector, wrong type :(
		math::XYZTLorentzVectorD const& p4t = part->p4();
		TLorentzVector p4(p4t.px(), p4t.py(), p4t.pz(), p4t.energy());
		TLorentzVector jetp4;
		jetp4.SetPtEtaPhiE(it->pt(), it->eta(), it->phi(), it->energy());

		if (id==PFCandidate::h) { // charged hadrons
		  nChargedHadrons += 1;
		  p4ChargedHadrons += p4;
		  sumptChargedHadrons_pfakt5[nJet_pfakt5] += p4.Pt();
		}
		if (id==PFCandidate::e) { // electrons
		  nElectrons += 1;
		  p4Electrons += p4;
		  sumptElectrons_pfakt5[nJet_pfakt5] += p4.Pt();
		}
		if (id==PFCandidate::mu) { // muons
		  nMuons += 1;
		  p4Muons += p4;
		  sumptMuons_pfakt5[nJet_pfakt5] += p4.Pt();
		}
		if (id==PFCandidate::gamma) { // photons
		  nPhotons += 1;
		  p4Photons += p4;
		  sumptPhotons_pfakt5[nJet_pfakt5] += p4.Pt();
		}
		if (id==PFCandidate::h0) { // neutral hadrons
		  nNeutralHadrons += 1;
		  p4NeutralHadrons += p4;
		  sumptNeutralHadrons_pfakt5[nJet_pfakt5] += p4.Pt();
		}
		if (id==PFCandidate::h_HF) { // HF hadrons
		  nHFHadrons += 1;
		  p4HFHadrons += p4;
		  sumptHFHadrons_pfakt5[nJet_pfakt5] += p4.Pt();
		}
		if (id==PFCandidate::egamma_HF) { // HF EM clusters
		  nHFEM += 1;
		  p4HFEM += p4;
		  sumptHFEM_pfakt5[nJet_pfakt5] += p4.Pt();
		}



		if (part.isNonnull())
		  itrk = (*part).trackRef();
		if (pt > pTMaxJet_pfakt5[nJet_pfakt5]) 
		  pTMaxJet_pfakt5[nJet_pfakt5] = pt;
		if (itrk.isNonnull() && pt > pTMaxChgJet_pfakt5[nJet_pfakt5]) 
		  pTMaxChgJet_pfakt5[nJet_pfakt5] = pt;
		if (!itrk.isNonnull() && pt > pTMaxNeutralJet_pfakt5[nJet_pfakt5]) 
		  pTMaxNeutralJet_pfakt5[nJet_pfakt5] = pt;
		if (!itrk.isNonnull() && pt > 1.0) 
		  nNeutral_ptCutJet_pfakt5[nJet_pfakt5]++;
		
		bool trkForAxis = false;
		bool trkForMult = false;
		
		//-----matching with vertex tracks-------
		if (!itrk.isNonnull()) { 
		  trkForMult = true;
		  trkForAxis = true;
		}
		else {
		  if (pt > 1.0)
		    nChg_ptCutJet_pfakt5[nJet_pfakt5]++;
		  float dZmin = 999;
		  int index_min = 999;
		  reco::VertexCollection::const_iterator vtxClose;
		  for(unsigned ivtx = 0;ivtx < VertexHandle->size();ivtx++) {
		    float dZ_cut = fabs(itrk->dz((*VertexHandle)[ivtx].position()));
		    float sumpT = 0;
		    for(reco::Vertex::trackRef_iterator itk = (*VertexHandle)[ivtx].tracks_begin();itk!=(*VertexHandle)[ivtx].tracks_end(); ++itk) {
		      sumpT = sumpT + ((*itk)->pt())*((*itk)->pt());
		    }
		    if (dZ_cut < dZmin) {
		      dZmin = dZ_cut;
		      index_min = ivtx;
		        //  std::cout<<"dz=="<<dZ_cut<<std::endl;
		    }
		  }//Loop over vertices 
		  if (index_min == 0) {
		    float dz = itrk->dz((*VertexHandle)[0].position());
		    float d0 = itrk->dxy((*VertexHandle)[0].position());
		    float vtx_xError = (*VertexHandle)[0].xError();
		    float vtx_yError = (*VertexHandle)[0].yError();
		    float vtx_zError = (*VertexHandle)[0].zError();
		    float d0_sigma=sqrt(pow(itrk->d0Error(),2) + pow(vtx_xError,2) + pow(vtx_yError,2));
		    float dz_sigma=sqrt(pow(itrk->dzError(),2) + pow(vtx_zError,2));
		    if (itrk->quality(reco::TrackBase::qualityByName("highPurity")) && fabs(dz/dz_sigma) < 5.) {
		      trkForAxis = true;
		      if (fabs(d0/d0_sigma) < 5.)
		        trkForMult = true;
		    }//
		  }
		  if (trkForMult)
		    nChg_QCJet_pfakt5[nJet_pfakt5]++;
		  if (itrk.isNonnull() && trkForMult && pt > 1.0)
		    nChg_ptCut_QCJet_pfakt5[nJet_pfakt5]++;
		  if (pt > pTMaxChg_QCJet_pfakt5[nJet_pfakt5] && trkForAxis) 
		    pTMaxChg_QCJet_pfakt5[nJet_pfakt5] = pt;
		}// for charged particles only


            jetPart_forMult.push_back(trkForMult);
            jetPart_forAxis.push_back(trkForAxis);

		double dphi = 2*atan(tan((phi-Phi0)/2));      
		double deta = eta-Eta0;
		SumW+=pt;
		SumW2+=pt*pt;
		SumDeta+=pt*pt*deta;
		SumDeta2+=pt*pt*deta*deta;
		SumDphi+=pt*pt*dphi;
		SumDphi2+=pt*pt*dphi*dphi;
		SumDetaDphi+=pt*pt*deta*dphi;
		float deltaR = jetp4.DeltaR(p4);
		rms_cands_wrong += (p4.Pt()*p4.Pt()*deltaR*deltaR);
		if (trkForAxis) {
		  SumW_QC+=pt;
		  SumW2_QC+=pt*pt;
		  SumDeta_QC+=pt*pt*deta;
		  SumDeta2_QC+=pt*pt*deta*deta;
		  SumDphi_QC+=pt*pt*dphi;
		  SumDphi2_QC+=pt*pt*dphi*dphi;
		  SumDetaDphi_QC+=pt*pt*deta*dphi;
		}

	}
	float ave_deta = SumDeta/SumW2;
	float ave_dphi = SumDphi/SumW2;
	float  ave_deta2 = SumDeta2/SumW2;
	float  ave_dphi2 = SumDphi2/SumW2;
      float a = ave_deta2-ave_deta*ave_deta;
      float b = ave_dphi2-ave_dphi*ave_dphi;
      float c = -(SumDetaDphi/SumW2-ave_deta*ave_dphi);
      float delta = sqrt(fabs((a-b)*(a-b)+4*c*c));
      if (a+b+delta > 0) {
        axis1Jet_pfakt5[nJet_pfakt5] = sqrt(0.5*(a+b+delta));
      }
      if (a+b-delta > 0) {  
        axis2Jet_pfakt5[nJet_pfakt5] = sqrt(0.5*(a+b-delta));
      }
      if (c != 0) {
        tanaJet_pfakt5[nJet_pfakt5] = 0.5*(b-a+delta)/c;
      }	

	float ave_deta_QC = SumDeta_QC/SumW2_QC;
	float ave_dphi_QC = SumDphi_QC/SumW2_QC;
	float  ave_deta2_QC = SumDeta2_QC/SumW2_QC;
	float  ave_dphi2_QC = SumDphi2_QC/SumW2_QC;
      float a_QC = ave_deta2_QC-ave_deta_QC*ave_deta_QC;
      float b_QC = ave_dphi2_QC-ave_dphi_QC*ave_dphi_QC;
      float c_QC = -(SumDetaDphi_QC/SumW2_QC-ave_deta_QC*ave_dphi_QC);
      float delta_QC = sqrt(fabs((a_QC-b_QC)*(a_QC-b_QC)+4*c_QC*c_QC));
      if (a_QC+b_QC+delta_QC > 0) {
        axis1_QCJet_pfakt5[nJet_pfakt5] = sqrt(0.5*(a_QC+b_QC+delta_QC));
      }
      if (a_QC+b_QC-delta_QC > 0) {  
        axis2_QCJet_pfakt5[nJet_pfakt5] = sqrt(0.5*(a_QC+b_QC-delta_QC));
      }
      if (c_QC != 0) {
        tana_QCJet_pfakt5[nJet_pfakt5] = 0.5*(b_QC-a_QC+delta_QC)/c_QC;
      }	

      ptDJet_pfakt5[nJet_pfakt5] =sqrt( SumW2/ (SumW*SumW));
      ptD_QCJet_pfakt5[nJet_pfakt5] =sqrt( SumW2_QC/ (SumW_QC*SumW_QC));

      // this is thw rong definition of rms, kept only for backwards compatibility:
	rmsCandJet_pfakt5[nJet_pfakt5] = rms_cands_wrong/SumW2;

	//-------calculate pull------
    	float ddetaR_sum(0.0), ddphiR_sum(0.0),ddetaR_sum_QC(0.0), ddphiR_sum_QC(0.0);
      float sum_ddR = 0.;
      float sum_ddR_QC = 0.;
    	for(int i=0; i<it->nConstituents(); ++i) {
			double pt=it->getJetConstituentsQuick()[i]->pt();
			double eta=it->getJetConstituentsQuick()[i]->eta();
			double phi=it->getJetConstituentsQuick()[i]->phi();
			double dphi = 2*atan(tan((phi-Phi0)/2));      
			double deta = eta-Eta0;
  		    float weight = pt*pt;
  		    float ddeta, ddphi,ddR;
  		    ddeta = deta - ave_deta ;//jetPart_deta[i] - ave_deta ; 
  		    ddphi = 2*atan(tan(( dphi - ave_dphi)/2.)) ;
  		    ddR = sqrt(ddeta*ddeta + ddphi*ddphi);
		    sum_ddR += ddR *ddR* weight;
  		    ddetaR_sum += ddR*ddeta*weight;
  		    ddphiR_sum += ddR*ddphi*weight;
                if (jetPart_forAxis[i]) { // this should be ave_deta_QC
  		      float ddeta_QC = deta - ave_deta_QC ;//jetPart_deta[i] - ave_deta ; 
  		      float ddphi_QC = 2*atan(tan(( dphi - ave_dphi_QC)/2.)) ;
  		      float ddR_QC = sqrt(ddeta_QC*ddeta_QC + ddphi_QC*ddphi_QC);
		      sum_ddR_QC += ddR_QC *ddR_QC* weight;
  		      ddetaR_sum_QC += ddR_QC*ddeta_QC*weight;
  		      ddphiR_sum_QC += ddR_QC*ddphi_QC*weight;
                }
  		  }//second loop over constituents  
  if (SumW2 > 0) {
    float ddetaR_ave = ddetaR_sum/SumW2;
    float ddphiR_ave = ddphiR_sum/SumW2;
    pullJet_pfakt5[nJet_pfakt5] = sqrt(ddetaR_ave*ddetaR_ave+ddphiR_ave*ddphiR_ave);
  }

  if (SumW2_QC > 0) {
    float ddetaR_ave_QC = ddetaR_sum_QC/SumW2_QC;
    float ddphiR_ave_QC = ddphiR_sum_QC/SumW2_QC;
    pull_QCJet_pfakt5[nJet_pfakt5] = sqrt(ddetaR_ave_QC*ddetaR_ave_QC+ddphiR_ave_QC*ddphiR_ave_QC);
  }

  rmsCandTrueJet_pfakt5[nJet_pfakt5] = sqrt( sum_ddR / SumW2);
  rmsCandTrue_QCJet_pfakt5[nJet_pfakt5] = sqrt( sum_ddR_QC / SumW2_QC);

  RchgJet_pfakt5[nJet_pfakt5] = pTMaxChgJet_pfakt5[nJet_pfakt5]/SumW;
  RneutralJet_pfakt5[nJet_pfakt5] = pTMaxNeutralJet_pfakt5[nJet_pfakt5]/SumW;
  RJet_pfakt5[nJet_pfakt5] = pTMaxJet_pfakt5[nJet_pfakt5]/SumW;
  Rchg_QCJet_pfakt5[nJet_pfakt5] = pTMaxChg_QCJet_pfakt5[nJet_pfakt5]/SumW_QC;




//	  float sumPt_cands=0.;
//	  float sumPt2_cands=0.;
//	  
//	  for (vector<PFCandidatePtr>::const_iterator jt = pfCandidates.begin();
//	       jt != pfCandidates.end(); ++jt) {
//	    
//	    PFCandidate::ParticleType id = (*jt)->particleId();
//	    // Convert particle momentum to normal TLorentzVector, wrong type :(
//	    math::XYZTLorentzVectorD const& p4t = (*jt)->p4();
//	    TLorentzVector p4(p4t.px(), p4t.py(), p4t.pz(), p4t.energy());
//	    TLorentzVector jetp4;
//	    jetp4.SetPtEtaPhiE(it->pt(), it->eta(), it->phi(), it->energy());
//	    if(p4.Pt()!=0){
//	      sumPt_cands += p4.Pt();
//	      sumPt2_cands += (p4.Pt()*p4.Pt());
//	      //float deltaR = it->p4().DeltaR(p4);
//	      float deltaR = jetp4.DeltaR(p4);
//	      rms_cands += (p4.Pt()*p4.Pt()*deltaR*deltaR);
//	    }

	    
	  //} //for PFCandidates

	  //ptDJet_pfakt5[nJet_pfakt5] = sqrt( sumPt2_cands )/sumPt_cands;
	  //rmsCandJet_pfakt5[nJet_pfakt5] = rms_cands/sumPt2_cands;

	  PileupJetIdentifier jetIdentifer_vars = jetMVACalculator.computeIdVariables(&correctedJet, scale, selectedVtx, vertexCollection);
	  jetId_dRMean_pfakt5[nJet_pfakt5]=jetIdentifer_vars.dRMean();
	  jetId_frac01_pfakt5[nJet_pfakt5]=jetIdentifer_vars.frac01();
	  jetId_frac02_pfakt5[nJet_pfakt5]=jetIdentifer_vars.frac02();
	  jetId_frac03_pfakt5[nJet_pfakt5]=jetIdentifer_vars.frac03();
	  jetId_frac04_pfakt5[nJet_pfakt5]=jetIdentifer_vars.frac04();
	  jetId_frac05_pfakt5[nJet_pfakt5]=jetIdentifer_vars.frac05();
	  jetId_nNeutrals_pfakt5[nJet_pfakt5]=jetIdentifer_vars.nNeutrals();
	  jetId_beta_pfakt5[nJet_pfakt5]=jetIdentifer_vars.beta();
	  jetId_betaStar_pfakt5[nJet_pfakt5]=jetIdentifer_vars.betaStar();
	  jetId_dZ_pfakt5[nJet_pfakt5]=jetIdentifer_vars.dZ();
	  jetId_nCharged_pfakt5[nJet_pfakt5]=jetIdentifer_vars.nCharged();
	  jetId_dR2Mean_pfakt5[nJet_pfakt5]=jetIdentifer_vars.dR2Mean();
	  jetId_betaStarClassic_pfakt5[nJet_pfakt5]=jetIdentifer_vars.betaStarClassic();
	  
	  for(unsigned int imva=0; imva<jetMVAAlgos.size(); imva++){
	    PileupJetIdAlgo* ialgo = (jetId_algos[imva]);
	    ialgo->set(jetIdentifer_vars);
	    PileupJetIdentifier id = ialgo->computeMva();
	    // if (jetMVAAlgos.size() != 4) cout << "problem with jet mva" << jetMVAAlgos.size() << endl;
	    // mva values
	    if (imva==0) jetIdSimple_mva_pfakt5[nJet_pfakt5]   = id.mva() ;
	    if (imva==1) jetIdFull_mva_pfakt5[nJet_pfakt5]     = id.mva() ;
	    if (imva==2) jetIdCutBased_mva_pfakt5[nJet_pfakt5] = id.mva() ;
	    // WP
	    if (imva==0) jetIdSimple_wp_pfakt5[nJet_pfakt5]   = id.idFlag() ;
	    if (imva==1) jetIdFull_wp_pfakt5[nJet_pfakt5]     = id.idFlag() ;
	    if (imva==2) jetIdCutBased_wp_pfakt5[nJet_pfakt5] = id.idFlag() ;
	  }	    

	  const TLorentzVector *p = 0;
	  
	  nChargedHadrons_pfakt5[nJet_pfakt5] =  nChargedHadrons;
	  p = &p4ChargedHadrons;
	  eChargedHadrons_pfakt5[nJet_pfakt5] = p->E() / it->energy();
	  if(p->E()){
	    ptChargedHadrons_pfakt5[nJet_pfakt5] = p->Pt();
	    etaChargedHadrons_pfakt5[nJet_pfakt5] = p->Eta();
	    phiChargedHadrons_pfakt5[nJet_pfakt5] = p->Phi();
	  }else{	    
	    ptChargedHadrons_pfakt5[nJet_pfakt5] = 0.;
	    etaChargedHadrons_pfakt5[nJet_pfakt5] = -999.;
	    phiChargedHadrons_pfakt5[nJet_pfakt5] = -999.;
	  }	
	  
	  nElectrons_pfakt5[nJet_pfakt5] =  nElectrons;
	  p = &p4Electrons;
	  eElectrons_pfakt5[nJet_pfakt5] = p->E() / it->energy();
	  if(p->E()){
	    ptElectrons_pfakt5[nJet_pfakt5] = p->Pt();
	    etaElectrons_pfakt5[nJet_pfakt5] = p->Eta();
	    phiElectrons_pfakt5[nJet_pfakt5] = p->Phi();
	  }else{
	    ptElectrons_pfakt5[nJet_pfakt5] = 0.;
	    etaElectrons_pfakt5[nJet_pfakt5] = -999.;
	    phiElectrons_pfakt5[nJet_pfakt5] = -999.;
	  }
	  
	  nMuons_pfakt5[nJet_pfakt5] =  nMuons;
	  p = &p4Muons;
	  eMuons_pfakt5[nJet_pfakt5] = p->E() / it->energy();
	  if(p->E()){
	    ptMuons_pfakt5[nJet_pfakt5] = p->Pt();
	    etaMuons_pfakt5[nJet_pfakt5] = p->Eta();
	    phiMuons_pfakt5[nJet_pfakt5] = p->Phi();
	  }else{
	    ptMuons_pfakt5[nJet_pfakt5] = 0.;
	    etaMuons_pfakt5[nJet_pfakt5] = -999.;
	    phiMuons_pfakt5[nJet_pfakt5] = -999.;
	  }
	  
	  nPhotons_pfakt5[nJet_pfakt5] =  nPhotons;
	  p = &p4Photons;
	  ePhotons_pfakt5[nJet_pfakt5] = p->E() / it->energy();
	  if(p->E()){	  
	    ptPhotons_pfakt5[nJet_pfakt5] = p->Pt();
	    etaPhotons_pfakt5[nJet_pfakt5] = p->Eta();
	    phiPhotons_pfakt5[nJet_pfakt5] = p->Phi();
	  }else{
	    ptPhotons_pfakt5[nJet_pfakt5] = 0.;
	    etaPhotons_pfakt5[nJet_pfakt5] = -999.;
	    phiPhotons_pfakt5[nJet_pfakt5] = -999.;
	  }
	  
	  nNeutralHadrons_pfakt5[nJet_pfakt5] =  nNeutralHadrons;
	  p = &p4NeutralHadrons;
	  eNeutralHadrons_pfakt5[nJet_pfakt5] = p->E() / it->energy();
	  if(p->E()){	  
	    ptNeutralHadrons_pfakt5[nJet_pfakt5] = p->Pt();
	    etaNeutralHadrons_pfakt5[nJet_pfakt5] = p->Eta();
	    phiNeutralHadrons_pfakt5[nJet_pfakt5] = p->Phi();
	  }else{
	    ptNeutralHadrons_pfakt5[nJet_pfakt5] = 0.;
	    etaNeutralHadrons_pfakt5[nJet_pfakt5] = -999.;
	    phiNeutralHadrons_pfakt5[nJet_pfakt5] = -999.;
	  }
	  
	  nHFHadrons_pfakt5[nJet_pfakt5] =  nHFHadrons;
	  p = &p4HFHadrons;
	  eHFHadrons_pfakt5[nJet_pfakt5] = p->E() / it->energy();
	  if(p->E()){	  
	    ptHFHadrons_pfakt5[nJet_pfakt5] = p->Pt();
	    etaHFHadrons_pfakt5[nJet_pfakt5] = p->Eta();
	    phiHFHadrons_pfakt5[nJet_pfakt5] = p->Phi();
	  }else{
	    ptHFHadrons_pfakt5[nJet_pfakt5] = 0.;
	    etaHFHadrons_pfakt5[nJet_pfakt5] = -999.;
	    phiHFHadrons_pfakt5[nJet_pfakt5] = -999.;
	  }	    
	  
	  nHFEM_pfakt5[nJet_pfakt5] =  nHFEM;
	  p = &p4HFEM;
	  eHFEM_pfakt5[nJet_pfakt5] = p->E() / it->energy();
	  if(p->E()){	  
	    ptHFEM_pfakt5[nJet_pfakt5] = p->Pt();
	    etaHFEM_pfakt5[nJet_pfakt5] = p->Eta();
	    phiHFEM_pfakt5[nJet_pfakt5] = p->Phi();
	  }else{
	    ptHFEM_pfakt5[nJet_pfakt5] = 0.;
	    etaHFEM_pfakt5[nJet_pfakt5] = -999.;
	    phiHFEM_pfakt5[nJet_pfakt5] = -999.;
	  }	    	  
	  
	  int index = nJet_pfakt5;
	  combinedSecondaryVertexBJetTags_pfakt5[index] = -999.;
	  combinedSecondaryVertexMVABJetTags_pfakt5[index] = -999.;
	  jetBProbabilityBJetTags_pfakt5[index] = -999.;
	  jetProbabilityBJetTags_pfakt5[index] =  -999.;
	  simpleSecondaryVertexHighEffBJetTags_pfakt5[index] =  -999.;
	  simpleSecondaryVertexHighPurBJetTags_pfakt5[index] =  -999.;
	  softMuonBJetTags_pfakt5[index] =  -999.;
	  softMuonByIP3dBJetTags_pfakt5[index] =  -999.;
	  softMuonByPtBJetTags_pfakt5[index] =  -999.;
	  //softElectronBJetTags_pfakt5[index] =  (*softElectronBJetTags)[index].second ;
	  softElectronByIP3dBJetTags_pfakt5[index] =   -999.;
	  softElectronByPtBJetTags_pfakt5[index]     = -999.;    
	  trackCountingHighPurBJetTags_pfakt5[index] = -999.;
	  trackCountingHighEffBJetTags_pfakt5[index] = -999.;
	  
	  //PU id
	  reco::TrackRefVector vTrks(it->getTrackRefs());
	  float sumTrkPt(0.0);
	  float sumTrkPtBeta[100],sumTrkPtBetaStar[100];
	  for (int ivtx=0;ivtx<100;++ivtx)
	    {
	      sumTrkPtBeta[ivtx]=0.;
	      sumTrkPtBetaStar[ivtx]=0.;
	    }

	  //---- loop over the tracks of the jet ----
	  for(reco::TrackRefVector::const_iterator i_trk = vTrks.begin(); i_trk != vTrks.end(); i_trk++) {
	    if ( VertexHandle->size() == 0) break;
	    sumTrkPt += (*i_trk)->pt();
	    int myVertex=-1;
	    //---- loop over all vertices ----------------------------
	    for(unsigned ivtx = 0;ivtx <  VertexHandle->size();ivtx++) {
	      //---- loop over the tracks associated with the vertex ---
	      if (!((* VertexHandle)[ivtx].isFake()) && (* VertexHandle)[ivtx].ndof() >= DEF_GOODVTX_NDOF && fabs((* VertexHandle)[ivtx].z()) <= DEF_GOODVTX_Z) {
		for(reco::Vertex::trackRef_iterator i_vtxTrk = (* VertexHandle)[ivtx].tracks_begin(); i_vtxTrk != (* VertexHandle)[ivtx].tracks_end(); ++i_vtxTrk) {
		  //---- match the jet track to the track from the vertex ----
		  reco::TrackRef trkRef(i_vtxTrk->castTo<reco::TrackRef>());
		  //---- check if the tracks match -------------------------
		  if (trkRef == (*i_trk)) 
		    {
		      myVertex=ivtx;
		      break;
		    }
		}
	      }
	    }
	    if (myVertex>-1)
	      {
		for(unsigned ivtx = 0;ivtx <  VertexHandle->size();ivtx++) {
		  if (ivtx== (unsigned int) myVertex)
		    sumTrkPtBeta[ivtx] += (*i_trk)->pt();
		  else
		    sumTrkPtBetaStar[ivtx] += (*i_trk)->pt();
		  if (ivtx == 0) {
		    nChargedHadronsgoodvtx += 1;
		    nChargedHadronsnoothervtx += 1;
		    TLorentzVector p4((*i_trk)->px(), (*i_trk)->py(), (*i_trk)->pz(), (*i_trk)->p());
		    p4ChargedHadronsgoodvtx += p4;
		    p4ChargedHadronsnoothervtx += p4;
		    sumptChargedHadronsgoodvtx_pfakt5[nJet_pfakt5] += p4.Pt();
		  }
		}
	      } else {
	      TLorentzVector p4((*i_trk)->px(), (*i_trk)->py(), (*i_trk)->pz(), (*i_trk)->p());
	      nChargedHadronsnoothervtx += 1;
	      p4ChargedHadronsnoothervtx += p4;
	      sumptChargedHadronsnoothervtx_pfakt5[nJet_pfakt5] += p4.Pt();
	    }

	  }
	  
	  nChargedHadronsgoodvtx_pfakt5[nJet_pfakt5] =  nChargedHadronsgoodvtx;
	  p = &p4ChargedHadronsgoodvtx;
	  eChargedHadronsgoodvtx_pfakt5[nJet_pfakt5] = p->E() / it->energy();
	  if(p->E()){
	    ptChargedHadronsgoodvtx_pfakt5[nJet_pfakt5] = p->Pt();
	    etaChargedHadronsgoodvtx_pfakt5[nJet_pfakt5] = p->Eta();
	    phiChargedHadronsgoodvtx_pfakt5[nJet_pfakt5] = p->Phi();
	  }else{	    
	    ptChargedHadronsgoodvtx_pfakt5[nJet_pfakt5] = 0.;
	    etaChargedHadronsgoodvtx_pfakt5[nJet_pfakt5] = -999.;
	    phiChargedHadronsgoodvtx_pfakt5[nJet_pfakt5] = -999.;
	  }	
	  
	  nChargedHadronsnoothervtx_pfakt5[nJet_pfakt5] =  nChargedHadronsnoothervtx;
	  p = &p4ChargedHadronsnoothervtx;
	  eChargedHadronsnoothervtx_pfakt5[nJet_pfakt5] = p->E() / it->energy();
	  if(p->E()){
	    ptChargedHadronsnoothervtx_pfakt5[nJet_pfakt5] = p->Pt();
	    etaChargedHadronsnoothervtx_pfakt5[nJet_pfakt5] = p->Eta();
	    phiChargedHadronsnoothervtx_pfakt5[nJet_pfakt5] = p->Phi();
	  }else{	    
	    ptChargedHadronsnoothervtx_pfakt5[nJet_pfakt5] = 0.;
	    etaChargedHadronsnoothervtx_pfakt5[nJet_pfakt5] = -999.;
	    phiChargedHadronsnoothervtx_pfakt5[nJet_pfakt5] = -999.;
	  }	
	
	  if (sumTrkPt > 0) {
	    for(unsigned ivtx = 0; ivtx <  VertexHandle->size();ivtx++) {
	      beta_pfakt5[index][ivtx]     = sumTrkPtBeta[ivtx]/sumTrkPt;
	      betaStar_pfakt5[index][ivtx] = sumTrkPtBetaStar[ivtx]/sumTrkPt;
	    }
	  }

	  if(it->pt() > 10. ){
	    // B tagging
	    combinedSecondaryVertexBJetTags_pfakt5[index] =  (*combinedSecondaryVertexBJetTags)[index].second ;
	    combinedSecondaryVertexMVABJetTags_pfakt5[index] =  (*combinedSecondaryVertexMVABJetTags)[index].second ;
	    jetBProbabilityBJetTags_pfakt5[index] =  (*jetBProbabilityBJetTags)[index].second ;
	    jetProbabilityBJetTags_pfakt5[index] =  (*jetProbabilityBJetTags)[index].second ;
	    simpleSecondaryVertexHighEffBJetTags_pfakt5[index] =  (*simpleSecondaryVertexHighEffBJetTags)[index].second ;
	    simpleSecondaryVertexHighPurBJetTags_pfakt5[index] =  (*simpleSecondaryVertexHighPurBJetTags)[index].second ;
	    
	    if( !softMuonBJetTags.isValid() || softMuonBJetTags->size() != pfjetsakt5->size() ) {
	      if(!softMuonBJetTags.isValid()) cout << "softMuonBJetTags not valid" << endl;
	      else   cout << "softMuonBJetTags: " << softMuonBJetTags->size() << " pfjetsakt5: " << pfjetsakt5->size() << endl;
	      cout << "jet index: " << index << endl;
	    } else {
	      softMuonBJetTags_pfakt5[index] =  (*softMuonBJetTags)[index].second ;
	    }
	    softMuonByIP3dBJetTags_pfakt5[index] =  (*softMuonByIP3dBJetTags)[index].second ;
	    softMuonByPtBJetTags_pfakt5[index] =  (*softMuonByPtBJetTags)[index].second ;
	    //softElectronBJetTags_pfakt5[index] =  (*softElectronBJetTags)[index].second ;
	    softElectronByIP3dBJetTags_pfakt5[index] =  (*softElectronByIP3dBJetTags)[index].second ;
	    softElectronByPtBJetTags_pfakt5[index] =  (*softElectronByPtBJetTags)[index].second;
	    trackCountingHighPurBJetTags_pfakt5[index] =  (*trackCountingHighPurBJetTags)[index].second ;
	    trackCountingHighEffBJetTags_pfakt5[index] =  (*trackCountingHighEffBJetTags)[index].second ;
	  }
	  ++nJet_pfakt5;
	  
	} // if >pfjetptthr     
      } // pfakt5
    }
  
  if (dumpPFAKT7Jets_)
    {
      for (PFJetCollection::const_iterator it = pfjetsakt7->begin(); 
	   it != pfjetsakt7->end(); ++it) {
	
	if (nJet_pfakt7>=100) {cout << "number of reco jets pfakt7 is larger than 100. Skipping" << endl; continue;}
	if (nJet_pfakt7 < pfjetnmin_ || it->pt() > pfjetptthr_) {
	  
	  ptJet_pfakt7[nJet_pfakt7] = it->pt();
	  eJet_pfakt7[nJet_pfakt7] = it->energy();	 
	  etaJet_pfakt7[nJet_pfakt7] = it->eta();	 
	  phiJet_pfakt7[nJet_pfakt7] = it->phi();	      
	  
	  // Jet Energy Scale Corrections on-the-fly     
	  PFJet  correctedJet = *it;
	  edm::RefToBase<reco::Jet> jetRef(edm::Ref<PFJetCollection>(pfjetsakt7,nJet_pfakt7));
	  double scale = corrector_pfakt5->correction(*it,iEvent,iSetup);
	  correctedJet.scaleEnergy(scale);
	  ptCorrJet_pfakt7[nJet_pfakt7] = correctedJet.pt();
	  
	  ++nJet_pfakt7;
	  
	} 
      }
    }

  if (dumpKT6Jets_)
    {
      for (PFJetCollection::const_iterator it = pfjetskt6->begin(); 
	   it != pfjetskt6->end(); ++it) {
	
	if (nJet_pfkt6>=100) {cout << "number of reco jets pfkt6 is larger than 100. Skipping" << endl; continue;}
	if (nJet_pfkt6 < pfjetnmin_ || it->pt() > pfjetptthr_) {
	  
	  ptJet_pfkt6[nJet_pfkt6] = it->pt();
	  eJet_pfkt6[nJet_pfkt6] = it->energy();	 
	  etaJet_pfkt6[nJet_pfkt6] = it->eta();	 
	  phiJet_pfkt6[nJet_pfkt6] = it->phi();	      
	  
	  nJet_pfkt6++;
	}
      }
    }

   // Fill caloMET
   const CaloMETCollection *calometcol = calomethandle.product();
   CaloMET const& calomet = calometcol->front();
   sMet = calomet.sumEt();
   eMet = calomet.pt();	 
   phiMet = calomet.phi();	      
   signifMet = calomet.mEtSig();	      
   
   // Fill muJEScaloMET
   // const CaloMETCollection *muJEScalometcol = muJESCorrMEThandle.product();
   // CaloMET const& muJEScalomet = muJEScalometcol->front();
   // smuCorrMet = muJEScalomet.sumEt();
   // emuCorrMet = muJEScalomet.pt();	 
   // phimuCorrMet = muJEScalomet.phi();	      
   // signifmuCorrMet = muJEScalomet.mEtSig();	      

   // Fill corrcaloMET
   const CaloMETCollection *mucalometcol = muCorrMEThandle.product();
   CaloMET const& mucalomet = mucalometcol->front();
   sCorrMet = mucalomet.sumEt();
   eCorrMet = mucalomet.pt();	 
   phiCorrMet = mucalomet.phi();	      
   signifCorrMet = mucalomet.mEtSig();	      

   // Fill caloMETNoHF
   const CaloMETCollection *calometcolNoHF = caloMEThandleNoHF.product();
   CaloMET const& calometNoHF = calometcolNoHF->front();
   sNoHFMet = calometNoHF.sumEt();
   eNoHFMet = calometNoHF.pt();	 
   phiNoHFMet = calometNoHF.phi();	      
   signifNoHFMet = calometNoHF.mEtSig();	      

   // Fill tcMET
   const View<MET> *tcmetcol = tcmethandle.product();
   MET const& tcmet = tcmetcol->front();
   stcMet = tcmet.sumEt();
   etcMet = tcmet.pt();
   phitcMet = tcmet.phi();
   signiftcMet = tcmet.mEtSig();	      

   // Fill pfMET

   const PFMETCollection *pfmetcol = pfmethandle.product();
   PFMET const& pfmet = pfmetcol->front();

   spfMet = pfmet.sumEt();
   epfMet = pfmet.pt();
   phipfMet = pfmet.phi();
   signifpfMet = pfmet.mEtSig();	      

   const PFMETCollection *pfMetType1Col = pfMetType1Handle.product();
   PFMET const& pfMetType1 = pfMetType1Col->front();

   spfMetType1 = pfMetType1.sumEt();
   epfMetType1 = pfMetType1.pt();
   phipfMetType1 = pfMetType1.phi();
   signifpfMetType1 = pfMetType1.mEtSig();	      

   edm::Handle<reco::PFMET> globalPfMet_h;
   iEvent.getByLabel(edm::InputTag("ClusteredPFMetProducerStd:globalPfMet"), globalPfMet_h) ;
   sglobalPfMet = globalPfMet_h->sumEt();
   eglobalPfMet = globalPfMet_h->pt();
   phiglobalPfMet = globalPfMet_h->phi();
   signifglobalPfMet = globalPfMet_h->mEtSig();	      

   edm::Handle<reco::PFMET> centralPfMet_h;
   iEvent.getByLabel(edm::InputTag("ClusteredPFMetProducerStd:centralPfMet"), centralPfMet_h) ;
   scentralPfMet = centralPfMet_h->sumEt();
   ecentralPfMet = centralPfMet_h->pt();
   phicentralPfMet = centralPfMet_h->phi();
   signifcentralPfMet = centralPfMet_h->mEtSig();	      

   edm::Handle< std::vector<reco::PFMET> > assocPfMet_h;
   iEvent.getByLabel(edm::InputTag("ClusteredPFMetProducerStd:assocPfMet"), assocPfMet_h) ;
   edm::Handle< std::vector<reco::PFMET> > assocOtherVtxPfMet_h;
   iEvent.getByLabel(edm::InputTag("ClusteredPFMetProducerStd:assocOtherVtxPfMet"), assocOtherVtxPfMet_h) ;
   edm::Handle< std::vector<reco::PFMET> > trkPfMet_h;
   iEvent.getByLabel(edm::InputTag("ClusteredPFMetProducerStd:trkPfMet"), trkPfMet_h) ;
   edm::Handle< std::vector<reco::PFMET> > cleanPfMet_h;
   iEvent.getByLabel(edm::InputTag("ClusteredPFMetProducerStd:cleanPfMet"), cleanPfMet_h) ;

   // new for 52X
   // edm::Handle< std::vector<reco::MET> > cleanedSaclayPfMet_h;
   // iEvent.getByLabel(edm::InputTag("cleanMETProducer:cleanPUMET"), cleanedSaclayPfMet_h) ;
   // edm::Handle< std::vector<reco::MET> > minTypeICleanSaclayPfMet_h;
   // iEvent.getByLabel(edm::InputTag("cleanMETProducer:minCPUMET"), minTypeICleanSaclayPfMet_h) ;

   for (unsigned int iVtx=0;iVtx<VertexHandle->size();++iVtx)
     {
       eassocPfMet[iVtx] = (*assocPfMet_h)[iVtx].pt();
       phiassocPfMet[iVtx] = (*assocPfMet_h)[iVtx].phi();
       signifassocPfMet[iVtx] = (*assocPfMet_h)[iVtx].mEtSig();	      

       eassocOtherVtxPfMet[iVtx] = (*assocOtherVtxPfMet_h)[iVtx].pt();
       phiassocOtherVtxPfMet[iVtx] = (*assocOtherVtxPfMet_h)[iVtx].phi();
       signifassocOtherVtxPfMet[iVtx] = (*assocOtherVtxPfMet_h)[iVtx].mEtSig();	      

       etrkPfMet[iVtx] = (*trkPfMet_h)[iVtx].pt();
       phitrkPfMet[iVtx] = (*trkPfMet_h)[iVtx].phi();
       signiftrkPfMet[iVtx] = (*trkPfMet_h)[iVtx].mEtSig();	      

       ecleanPfMet[iVtx] = (*cleanPfMet_h)[iVtx].pt();
       phicleanPfMet[iVtx] = (*cleanPfMet_h)[iVtx].phi();
       signifcleanPfMet[iVtx] = (*cleanPfMet_h)[iVtx].mEtSig();	      

       /*
       ecleanedSaclayPfMet[iVtx] = (*cleanedSaclayPfMet_h)[iVtx].pt();
       phicleanedSaclayPfMet[iVtx] = (*cleanedSaclayPfMet_h)[iVtx].phi();
       signifcleanedSaclayPfMet[iVtx] = (*cleanedSaclayPfMet_h)[iVtx].mEtSig();	      

       eminTypeICleanSaclayPfMet[iVtx] = (*minTypeICleanSaclayPfMet_h)[iVtx].pt();
       phiminTypeICleanSaclayPfMet[iVtx] = (*minTypeICleanSaclayPfMet_h)[iVtx].phi();
       signifminTypeICleanSaclayPfMet[iVtx] = (*minTypeICleanSaclayPfMet_h)[iVtx].mEtSig();	      
       */
     }

   edm::Handle< std::vector<double> >globalPfMetSums_h; 
   iEvent.getByLabel(edm::InputTag("ClusteredPFMetProducerStd:globalPfMetSums"), globalPfMetSums_h) ;

   for (unsigned int iVtx=0;iVtx<12;++iVtx)
       globalPfSums[iVtx] = (*globalPfMetSums_h)[iVtx];


   sMetGen = 0.;
   eMetGen = 0.;
   phiMetGen = 0.;
   signifMetGen = 0.;

   sMetGen2 = 0.;
   eMetGen2 = 0.;
   phiMetGen2 = 0.;

   if( isMC ) {

     // Fill gen MET

     const GenMETCollection *genmetcol = genmethandle.product();
     GenMET const& genmet = genmetcol->front();

     sMetGen = genmet.sumEt();
     eMetGen = genmet.energy();
     phiMetGen = genmet.phi();
     signifMetGen = genmet.mEtSig();

     const GenMETCollection *genmetcol2 = genmethandle2.product();
     GenMET const& genmet2 = genmetcol2->front();

     sMetGen2 = genmet2.sumEt();
     eMetGen2 = genmet2.energy();
     phiMetGen2 = genmet2.phi();
     
   } //if is MC

       edm::Handle<reco::MuonCollection> muonHandle;
       iEvent.getByLabel("muons",muonHandle);

       nMuonsReco = 0;
       //       std::cout <<  muonHandle->size() << std::endl;
       for(int x=0;x < min((int) muonHandle->size(),200); x++)
	 {
	   //	   std::cout << "+++" << x << std::endl;
	   ++nMuonsReco;
	   muon_pt[x]  = (*muonHandle)[x].pt();
	   muon_energy[x]  = (*muonHandle)[x].energy();
	   muon_px[x]  = (*muonHandle)[x].px();
	   muon_py[x]  = (*muonHandle)[x].py();
	   muon_pz[x]  = (*muonHandle)[x].pz();
	   muon_phi[x] = (*muonHandle)[x].phi();
	   muon_eta[x] = (*muonHandle)[x].eta();
	   muon_charge[x] = (*muonHandle)[x].charge();
	   muon_vx[x]  = (*muonHandle)[x].vx();
	   muon_vy[x]  = (*muonHandle)[x].vy();
	   muon_vz[x]  = (*muonHandle)[x].vz();

	   // PF isolation (for 2012 data) - chiara
	   if (IS2012) {
	     muon_pfiso04_chHad[x] = (*muonHandle)[x].pfIsolationR04().sumChargedHadronPt;
	     muon_pfiso04_chPar[x] = (*muonHandle)[x].pfIsolationR04().sumChargedParticlePt;
	     muon_pfiso04_nHad[x]  = (*muonHandle)[x].pfIsolationR04().sumNeutralHadronEt;
	     muon_pfiso04_Phot[x]  = (*muonHandle)[x].pfIsolationR04().sumPhotonEt;
	     muon_pfiso04_PUPt[x]  = (*muonHandle)[x].pfIsolationR04().sumPUPt;
	     muon_pfiso03_chHad[x] = (*muonHandle)[x].pfIsolationR03().sumChargedHadronPt;
	     muon_pfiso03_chPar[x] = (*muonHandle)[x].pfIsolationR03().sumChargedParticlePt;
	     muon_pfiso03_nHad[x]  = (*muonHandle)[x].pfIsolationR03().sumNeutralHadronEt;
	     muon_pfiso03_Phot[x]  = (*muonHandle)[x].pfIsolationR03().sumPhotonEt;
	     muon_pfiso03_PUPt[x]  = (*muonHandle)[x].pfIsolationR03().sumPUPt;
	     muon_isPFMuon[x]      = (*muonHandle)[x].isPFMuon();   
	   } else {
	     muon_pfiso04_chHad[x] = -999.;
	     muon_pfiso04_chPar[x] = -999.;
	     muon_pfiso04_nHad[x]  = -999.;
	     muon_pfiso04_Phot[x]  = -999.;
	     muon_pfiso04_PUPt[x]  = -999.;
	     muon_pfiso03_chHad[x] = -999.;
	     muon_pfiso03_chPar[x] = -999.;
	     muon_pfiso03_nHad[x]  = -999.;
	     muon_pfiso03_Phot[x]  = -999.;
	     muon_pfiso03_PUPt[x]  = -999.;
	     muon_isPFMuon[x]      = -999;
	   }

	   //tia's stuff
	   muon_isGlobalMuon[x] =     (*muonHandle)[x].isGlobalMuon();
	   muon_isTrackerMuon[x] =    (*muonHandle)[x].isTrackerMuon();
	   muon_isStandAloneMuon[x] = (*muonHandle)[x].isStandAloneMuon();
	   
	   //FIXME: Seems needed  top and bottom referene point, but not sure how to do that !!!
	   reco::TrackRef moTrkref;
	   if((muon_isGlobalMuon[x]) || (muon_isTrackerMuon[x])){
	     moTrkref = (*muonHandle)[x].innerTrack();
	     muon_InnerTrack_isNonnull[x] = (*muonHandle)[x].innerTrack().isNonnull();
	   }
	   else{
	     moTrkref = (*muonHandle)[x].outerTrack();
	     muon_OuterTrack_isNonnull[x] =   (*muonHandle)[x].outerTrack().isNonnull();
	   }
	   //For StandAlone
	   muon_OuterPoint_x[x]=0.0;
	   muon_OuterPoint_y[x]=0.0;
	   muon_OuterPoint_z[x]=0.0;
	   //For Global,Tracker
	   muon_InnerPoint_x[x]=0.0;
	   muon_InnerPoint_y[x]=0.0;
	   muon_InnerPoint_z[x]=0.0;
	   
	   if((muon_OuterTrack_isNonnull[x])){//stand-alone
	     muon_OuterPoint_x[x]= moTrkref->referencePoint().x();
	     muon_OuterPoint_y[x]= moTrkref->referencePoint().y();
	     muon_OuterPoint_z[x]= moTrkref->referencePoint().z();
	   }
	   if(muon_InnerTrack_isNonnull[x]){//global,tracker
	     muon_InnerPoint_x[x]= moTrkref->referencePoint().x();
	     muon_InnerPoint_y[x]= moTrkref->referencePoint().y();
	     muon_InnerPoint_z[x]= moTrkref->referencePoint().z();
	   }
	   // default isolation variables 0.3
	   MuonIsolation Iso03  = (*muonHandle)[x].isolationR03(); // Default  variables 0.3
	   muon_trackIso[x] = Iso03.sumPt;
	   muon_ecalIso[x] = Iso03.emEt;
	   muon_hcalIso[x] = Iso03.hadEt;
	   muon_relIso[x] =(muon_trackIso[x] + muon_ecalIso[x] + muon_hcalIso[x])/muon_pt[x];

	   muon_normChi2[x]= -99;
	   muon_validHits[x]= -99;
	   
	   if((*muonHandle)[x].globalTrack().isNonnull() ){
	     muon_normChi2[x] = (*muonHandle)[x].globalTrack()->chi2()/(*muonHandle)[x].globalTrack()->ndof();
	     muon_validHits[x] = (*muonHandle)[x].globalTrack()->hitPattern().numberOfValidMuonHits();
	   }
	   
	   muon_numberOfMatches[x] = (*muonHandle)[x].numberOfMatches();

	   muon_tkHits[x] =-99;
	   muon_pixHits[x]=-99;

	   if((*muonHandle)[x].track().isNonnull() ){
	     muon_tkHits[x] = (*muonHandle)[x].track()->numberOfValidHits();
	     muon_pixHits[x] = (*muonHandle)[x].track()->hitPattern().numberOfValidPixelHits();
	     muon_trkLayerWithMeas[x] = (*muonHandle)[x].track()->hitPattern().trackerLayersWithMeasurement();  // chiara
	   }
	   

	 }

   if (dumpBeamHaloInformations_)
     {
       edm::Handle<BeamHaloSummary> TheBeamHaloSummary ;
       iEvent.getByLabel("BeamHaloSummary", TheBeamHaloSummary) ;


       isBeamHaloIDTightPass    = false;
       isBeamHaloIDLoosePass    = false;

       isBeamHaloEcalLoosePass   = false;
       isBeamHaloHcalLoosePass   = false;
       isBeamHaloCSCLoosePass    = false;
       isBeamHaloGlobalLoosePass = false;

       isBeamHaloEcalTightPass   = false;
       isBeamHaloHcalTightPass   = false;
       isBeamHaloCSCTightPass    = false;
       isBeamHaloGlobalTightPass = false;


       isSmellsLikeHalo_Tag  = false;
       isLooseHalo_Tag       = false;
       isTightHalo_Tag       = false;
       isExtremeTightHalo_Tag = false;


       if(TheBeamHaloSummary.isValid()) {
	 const BeamHaloSummary TheSummary = (*TheBeamHaloSummary.product() );

	 if(TheSummary.EcalLooseHaloId())isBeamHaloEcalLoosePass   = true;
	 if(TheSummary.HcalLooseHaloId())isBeamHaloHcalLoosePass   = true;

	 if(TheSummary.EcalTightHaloId())isBeamHaloEcalTightPass   = true;
	 if(TheSummary.HcalTightHaloId())isBeamHaloHcalTightPass   = true;

	 if(TheSummary.CSCLooseHaloId())isBeamHaloCSCLoosePass    = true;
	 if(TheSummary.CSCTightHaloId())isBeamHaloCSCTightPass    = true;

	 if(TheSummary.GlobalLooseHaloId())isBeamHaloGlobalLoosePass = true;
	 if(TheSummary.GlobalTightHaloId())isBeamHaloGlobalTightPass = true;

	 if(TheSummary.EventSmellsLikeHalo())isSmellsLikeHalo_Tag = true;
	 if(TheSummary.LooseId())isLooseHalo_Tag = true;
	 if(TheSummary.TightId())isTightHalo_Tag = true;
	 if(TheSummary.ExtremeTightId())isExtremeTightHalo_Tag = true;


	 if( TheSummary.EcalLooseHaloId()  && TheSummary.HcalLooseHaloId() &&
	     TheSummary.CSCLooseHaloId()   && TheSummary.GlobalLooseHaloId() )
	   isBeamHaloIDLoosePass = true;

	 if( TheSummary.EcalTightHaloId()  && TheSummary.HcalTightHaloId() &&
	     TheSummary.CSCTightHaloId()   && TheSummary.GlobalTightHaloId() )
	   isBeamHaloIDTightPass = true;

       }//not empty




       //cosmic muon
       edm::Handle<reco::MuonCollection>cosmicMuonHandle;
       iEvent.getByLabel("muonsFromCosmics",cosmicMuonHandle);
       const reco::MuonCollection & cosmicMuons = *cosmicMuonHandle;

       nCosmicMuons=0;
       for(int x=0;x < min((int)cosmicMuons.size(),200);x++){

	 cosmicmuon_pt[x]  = cosmicMuons[x].pt();
	 cosmicmuon_energy[x]  = cosmicMuons[x].energy();
	 cosmicmuon_px[x]  = cosmicMuons[x].px();
	 cosmicmuon_py[x]  = cosmicMuons[x].py();
	 cosmicmuon_pz[x]  = cosmicMuons[x].pz();
	 cosmicmuon_phi[x] = cosmicMuons[x].phi();
	 cosmicmuon_eta[x] = cosmicMuons[x].eta();
	 cosmicmuon_charge[x] = cosmicMuons[x].charge();

	 //tia's stuff
	 cosmicmuon_isGlobalMuon[x] =     cosmicMuons[x].isGlobalMuon();
	 cosmicmuon_isTrackerMuon[x] =    cosmicMuons[x].isTrackerMuon();
	 cosmicmuon_isStandAloneMuon[x] = cosmicMuons[x].isStandAloneMuon();


	 reco::TrackRef cmoTrkref = cosmicMuons[x].outerTrack();
	 cosmicmuon_OuterTrack_isNonnull[x] = cosmicMuons[x].outerTrack().isNonnull();
	 
	 cosmicmuon_OuterPoint_x[x]=0.0;
	 cosmicmuon_OuterPoint_y[x]=0.0;
	 cosmicmuon_OuterPoint_z[x]=0.0;
	 // standalone muon variables
	 if(cosmicmuon_OuterTrack_isNonnull[x]){
	   cosmicmuon_OuterPoint_x[x]= cmoTrkref->referencePoint().x();
	   cosmicmuon_OuterPoint_y[x]= cmoTrkref->referencePoint().y();
	   cosmicmuon_OuterPoint_z[x]= cmoTrkref->referencePoint().z();
	 }


	 nCosmicMuons++;
       }//end of for loop

     }
   //   std::cout << "===" << nMuonsReco << std::endl;
   //event++;  
   m_tree->Fill();
}

std::map<std::string,float> GammaJetAnalyzer::particleLevelIsolation(const GenParticleCollection* genParticles, const GenParticle* particle, float DRsize, int status, float DRVetoSize, float chargedPtMin, float neutralEMPtMin, float neutralHadPtMin)
{
  std::map<std::string,float> sumPt;

  sumPt["Charged"]=0.;
  sumPt["EMNeutral"]=0.;
  sumPt["HADNeutral"]=0.;

  for (GenParticleCollection::const_iterator p = genParticles->begin();
       p != genParticles->end(); ++p) {
       if (p->status() != status) 
	 continue;

       if (reco::deltaR(*p,*particle)>DRsize)
	 continue;
       
       if ((reco::deltaR(*p,*particle)<DRVetoSize && p->pdgId()==particle->pdgId()) || ( &(*p) == particle ) || (p->mother() == particle) ) //remove self- veto
	 continue;
       
       if (abs(p->charge())!=0 && p->pt()>chargedPtMin)
	 sumPt["Charged"]+=p->pt();
       else if (abs(p->pdgId())==22 && p->pt()>neutralEMPtMin)
	 sumPt["EMNeutral"]+=p->pt();
       else if (p->pt()>neutralHadPtMin)
       	 sumPt["HADNeutral"]+=p->pt();
  }
  
  return sumPt;
}
// ------------ method called once each job just before starting event loop  ------------
void 
GammaJetAnalyzer::beginJob()
{

  //m_tree = fs_->make<TTree>("pippo","Analysis tree");
  outfile = TFile::Open(outFileName.c_str(), "RECREATE");
  outfile->mkdir("myanalysis");
  outfile->cd("myanalysis");

  m_tree = new TTree ("pippo","Analysis tree") ;
  //  m_tree->SetAutoSave (10000000) ;
  m_tree->Branch("genpt",&genpt,"genpt/F");
  m_tree->Branch("genProcessId",&genProcessId,"genProcessId/I");
  m_tree->Branch("genQScale",&genQScale,"genQScale/F");

  m_tree->Branch("isMC",&isMC,"isMC/O");

  //MET Filters
  m_tree->Branch("passEcalLaserFilter",&passEcalLaserFilter,"passEcalLaserFilter/O");
  m_tree->Branch("passHBHENoiseFilter",&passHBHENoiseFilter,"passHBHENoiseFilter/O");
  m_tree->Branch("passCSCTightHaloFilter",&passCSCTightHaloFilter,"passCSCTightHaloFilter/O");
  m_tree->Branch("passhcalLaserEventFilter",&passhcalLaserEventFilter,"passhcalLaserEventFilter/O");
  m_tree->Branch("passEcalDeadCellTriggerPrimitiveFilter",&passEcalDeadCellTriggerPrimitiveFilter,"passEcalDeadCellTriggerPrimitiveFilter/O");
  m_tree->Branch("passtrackingFailureFilter",&passtrackingFailureFilter,"passtrackingFailureFilter/O");
  m_tree->Branch("passeeBadScFilter",&passeeBadScFilter,"passeeBadScFilter/O");

  m_tree->Branch("store",&store,"store/I");
  m_tree->Branch("lbn",&lbn,"lbn/I");
  m_tree->Branch("bx",&bx,"bx/I");
  m_tree->Branch("orbit",&orbit,"orbit/I");
  m_tree->Branch("run",&run,"run/I");
  m_tree->Branch("event",&event,"event/I");

  m_tree->Branch("rhoPF",&rho,"rhoPF/F");
  m_tree->Branch("rhoCalo",&rhoCalo,"rhoCalo/F");
  m_tree->Branch("rhoAllJets",&rhoAllJets,"rhoAllJets/F");

  // Problem: nMC==100 always, and sometimes last particle has very high pT
  // => could be losing interesting particles, even quarks/gluons (status==2)
  // hmmm, status==1 particles have pT less than ~4.5 GeV for pThat>500 => ok
  m_tree->Branch("nMC",&nMC,"nMC/I");
  m_tree->Branch("pdgIdMC",&pdgIdMC,"pdgIdMC[nMC]/I");
  m_tree->Branch("statusMC",&statusMC,"statusMC[nMC]/I");
  m_tree->Branch("motherIDMC",&motherIDMC,"motherIDMC[nMC]/I");
  // Most MC particles have mass, but why do photons (status=1 and 3) have mass?
  //m_tree->Branch("massMC ",&massMC ,"massMC[nMC]/F");
  //to a good approximation, m = sqrt(e^2 - (pt*cosh(eta))^2), when m>1e-6 GeV
  m_tree->Branch("ptMC ",&ptMC ,"ptMC[nMC]/F");
  m_tree->Branch("eMC  ",&eMC  ,"eMC[nMC]/F");
  m_tree->Branch("etaMC",&etaMC,"etaMC[nMC]/F");
  m_tree->Branch("phiMC",&phiMC,"phiMC[nMC]/F");

  m_tree->Branch("isoParticleChargedDR01MC",&isoParticleChargedDR01MC,"isoParticleChargedDR01MC[nMC]/F");
  m_tree->Branch("isoParticleChargedDR02MC",&isoParticleChargedDR02MC,"isoParticleChargedDR02MC[nMC]/F");
  m_tree->Branch("isoParticleChargedDR03MC",&isoParticleChargedDR03MC,"isoParticleChargedDR03MC[nMC]/F");
  m_tree->Branch("isoParticleChargedDR04MC",&isoParticleChargedDR04MC,"isoParticleChargedDR04MC[nMC]/F");
  m_tree->Branch("isoParticleChargedDR05MC",&isoParticleChargedDR05MC,"isoParticleChargedDR05MC[nMC]/F");

  m_tree->Branch("isoParticleEMNeutralDR01MC",&isoParticleEMNeutralDR01MC,"isoParticleEMNeutralDR01MC[nMC]/F");
  m_tree->Branch("isoParticleEMNeutralDR02MC",&isoParticleEMNeutralDR02MC,"isoParticleEMNeutralDR02MC[nMC]/F");
  m_tree->Branch("isoParticleEMNeutralDR03MC",&isoParticleEMNeutralDR03MC,"isoParticleEMNeutralDR03MC[nMC]/F");
  m_tree->Branch("isoParticleEMNeutralDR04MC",&isoParticleEMNeutralDR04MC,"isoParticleEMNeutralDR04MC[nMC]/F");
  m_tree->Branch("isoParticleEMNeutralDR05MC",&isoParticleEMNeutralDR05MC,"isoParticleEMNeutralDR05MC[nMC]/F");

  m_tree->Branch("isoParticleHADNeutralDR01MC",&isoParticleHADNeutralDR01MC,"isoParticleHADNeutralDR01MC[nMC]/F");
  m_tree->Branch("isoParticleHADNeutralDR02MC",&isoParticleHADNeutralDR02MC,"isoParticleHADNeutralDR02MC[nMC]/F");
  m_tree->Branch("isoParticleHADNeutralDR03MC",&isoParticleHADNeutralDR03MC,"isoParticleHADNeutralDR03MC[nMC]/F");
  m_tree->Branch("isoParticleHADNeutralDR04MC",&isoParticleHADNeutralDR04MC,"isoParticleHADNeutralDR04MC[nMC]/F");
  m_tree->Branch("isoParticleHADNeutralDR05MC",&isoParticleHADNeutralDR05MC,"isoParticleHADNeutralDR05MC[nMC]/F");

  m_tree->Branch("isoPartonDR01MC",&isoPartonDR01MC,"isoPartonDR01MC[nMC]/F");
  m_tree->Branch("isoPartonDR02MC",&isoPartonDR02MC,"isoPartonDR02MC[nMC]/F");
  m_tree->Branch("isoPartonDR03MC",&isoPartonDR03MC,"isoPartonDR03MC[nMC]/F");
  m_tree->Branch("isoPartonDR04MC",&isoPartonDR04MC,"isoPartonDR04MC[nMC]/F");
  m_tree->Branch("isoPartonDR05MC",&isoPartonDR05MC,"isoPartonDR05MC[nMC]/F");

  m_tree->Branch("pu_n", &pu_n, "pu_n/I");
  m_tree->Branch("pu_true_n", &pu_true_n, "pu_true_n/I");
  //m_tree->Branch("pu_bunchcrossing", &pu_bunchcrossing, "pu_bunchcrossing/I");
  m_tree->Branch("pu_zpos", &pu_zpos, "pu_zpos[pu_n]/F");
  m_tree->Branch("pu_sumpt_lowpt", &pu_sumpt_lowpt, "pu_sumpt_lowpt[pu_n]/F");
  m_tree->Branch("pu_sumpt_highpt", &pu_sumpt_highpt, "pu_sumpt_highpt[pu_n]/F");
  m_tree->Branch("pu_ntrks_lowpt", &pu_ntrks_lowpt, "pu_ntrks_lowpt[pu_n]/F");
  m_tree->Branch("pu_ntrks_highpt", &pu_ntrks_highpt, "pu_ntrks_highpt[pu_n]/F");


  m_tree->Branch("nPhot",&nPhot,"nPhot/I");
  m_tree->Branch("ptPhot ",&ptPhot ,"ptPhot[nPhot]/F");
  m_tree->Branch("ePhot  ",&ePhot  ,"ePhot[nPhot]/F");
  m_tree->Branch("escPhot  ",&escPhot  ,"escPhot[nPhot]/F");
  m_tree->Branch("escRegrPhot  ",&escRegrPhot  ,"escRegrPhot[nPhot]/F");
  m_tree->Branch("escRegrPhotError  ",&escRegrPhotError  ,"escRegrPhotError[nPhot]/F");
  m_tree->Branch("escPhFixPhot  ",&escPhFixPhot  ,"escPhFixPhot[nPhot]/F");
  m_tree->Branch("escPhFixPhotError  ",&escPhFixPhotError  ,"escPhFixPhotError[nPhot]/F");
  m_tree->Branch("escRawPhot  ",&escRawPhot  ,"escRawPhot[nPhot]/F");
  m_tree->Branch("etascPhot  ",&etascPhot  ,"etascPhot[nPhot]/F");
  m_tree->Branch("phiscPhot  ",&phiscPhot  ,"phiscPhot[nPhot]/F");
  m_tree->Branch("xscPhot  ",&xscPhot  ,"xscPhot[nPhot]/F");
  m_tree->Branch("yscPhot  ",&yscPhot  ,"yscPhot[nPhot]/F");
  m_tree->Branch("zscPhot  ",&zscPhot  ,"zscPhot[nPhot]/F");
  m_tree->Branch("xcaloPhot  ",&xcaloPhot  ,"xcaloPhot[nPhot]/F");
  m_tree->Branch("ycaloPhot  ",&ycaloPhot  ,"ycaloPhot[nPhot]/F");
  m_tree->Branch("zcaloPhot  ",&zcaloPhot  ,"zcaloPhot[nPhot]/F");
  m_tree->Branch("eseedPhot  ",&eseedPhot  ,"eseedPhot[nPhot]/F");
  m_tree->Branch("etaPhot",&etaPhot,"etaPhot[nPhot]/F");
  m_tree->Branch("phiPhot",&phiPhot,"phiPhot[nPhot]/F");
  m_tree->Branch("timePhot",&timePhot,"timePhot[nPhot]/F");
  m_tree->Branch("e4SwissCrossPhot",&e4SwissCrossPhot,"e4SwissCrossPhot[nPhot]/F");
  m_tree->Branch("hasPixelSeedPhot",&hasPixelSeedPhot,"hasPixelSeedPhot[nPhot]/I");
  m_tree->Branch("hasMatchedPromptElePhot",&hasMatchedPromptElePhot,"hasMatchedPromptElePhot[nPhot]/I");
  m_tree->Branch("hasMatchedConvPhot",&hasMatchedConvPhot,"hasMatchedConvPhot[nPhot]/I");
  m_tree->Branch("isEBPhot",&isEBPhot, "isEBPhot[nPhot]/O");
  m_tree->Branch("isEEPhot",&isEEPhot, "isEEPhot[nPhot]/O");
  m_tree->Branch("isEBEEGapPhot",&isEBEEGapPhot, "isEBEEGapPhot[nPhot]/O");

  m_tree->Branch("ntracksConvPhot",nTracksConvPhot,"nTracksConvPhot[nPhot]/I");
  m_tree->Branch("isValidVtxConvPhot",isValidVtxConvPhot,"isValidVtxConvPhot[nPhot]/O");
  m_tree->Branch("pairInvmassConvPhot",pairInvariantMassConvPhot,"pairInvariantMassConvPhot[nPhot]/F");
  m_tree->Branch("pairCotThetaSeperationConvPhot",pairCotThetaSeparationConvPhot,"pairCotThetaSeparationConvPhot[nPhot]/F");
  m_tree->Branch("pairmomentumXConvPhot",pairMomentum_xConvPhot,"pairMomentum_xConvPhot[nPhot]/F");
  m_tree->Branch("pairmomentumYConvPhot",pairMomentum_yConvPhot,"pairMomentum_yConvPhot[nPhot]/F");
  m_tree->Branch("pairmomentumZConvPhot",pairMomentum_zConvPhot,"pairMomentum_zConvPhot[nPhot]/F");
  m_tree->Branch("chi2ConvPhot",chi2ConvPhot,"chi2ConvPhot[nPhot]/F");
  m_tree->Branch("nDofConvPhot",nDofConvPhot,"nDofConvPhot[nPhot]/F");
  m_tree->Branch("eOverPConvPhot",eOverPConvPhot,"eOverPConvPhot[nPhot]/F");
  m_tree->Branch("convVxConvPhot",conv_vxConvPhot,"conv_vxConvPhot[nPhot]/F");
  m_tree->Branch("convVyConvPhot",conv_vyConvPhot,"conv_vyConvPhot[nPhot]/F");
  m_tree->Branch("convVzConvPhot",conv_vzConvPhot,"conv_vzConvPhot[nPhot]/F");
  m_tree->Branch("distOfMinimumApproachConvPhot",distOfMinimumApproachConvPhot,"distOfMinimumApproachConvPhot[nPhot]/F");
  m_tree->Branch("dPhiTracksAtVtxConvPhot",dPhiTracksAtVtxConvPhot,"dPhiTracksAtVtxConvPhot[nPhot]/F");
//   m_tree->Branch("dPhiTracksAtEcalConvPhot",dPhiTracksAtEcalConvPhot,"dPhiTracksAtEcalConvPhot[nPhot]/F");
//   m_tree->Branch("dEtaTracksAtEcalConvPhot",dEtaTracksAtEcalConvPhot,"dEtaTracksAtEcalConvPhot[nPhot]/F");

  // Default photon ID
  m_tree->Branch("pid_isEM",&pid_isEM,"pid_isEM[nPhot]/O");
  m_tree->Branch("pid_isLoose",&pid_isLoose,"pid_isLoose[nPhot]/O");
  m_tree->Branch("pid_isTight",&pid_isTight,"pid_isTight[nPhot]/O");
  m_tree->Branch("pid_jurECAL",&pid_jurECAL,"pid_jurECAL[nPhot]/F");
  m_tree->Branch("pid_twrHCAL",&pid_twrHCAL,"pid_twrHCAL[nPhot]/F");
  m_tree->Branch("pid_HoverE",&pid_HoverE,"pid_HoverE[nPhot]/F");
  m_tree->Branch("pid_hlwTrack",&pid_hlwTrack,"pid_hlwTarck[nPhot]/F");
  m_tree->Branch("pid_hlwTrackNoDz",&pid_hlwTrackNoDz,"pid_hlwTrackNoDz[nPhot]/F");
  m_tree->Branch("pid_hlwTrackForCiC",&pid_hlwTrackForCiC,"pid_hlwTrackBestRank[40][100]/F");
  m_tree->Branch("pid_etawid",&pid_etawid,"pid_etawid[nPhot]/F");

  m_tree->Branch("pid_jurECAL03",&pid_jurECAL03,"pid_jurECAL03[nPhot]/F");
  m_tree->Branch("pid_twrHCAL03",&pid_twrHCAL03,"pid_twrHCAL03[nPhot]/F");
  m_tree->Branch("pid_hlwTrack03",&pid_hlwTrack03,"pid_hlwTrack03[nPhot]/F");
  m_tree->Branch("pid_hlwTrack03NoDz",&pid_hlwTrack03NoDz,"pid_hlwTrack03NoDz[nPhot]/F");
  m_tree->Branch("pid_hlwTrack03ForCiC",&pid_hlwTrack03ForCiC,"pid_hlwTrack03ForCiC[40][100]/F");

  m_tree->Branch("pid_pfIsoCharged01ForCiC",&pid_pfIsoCharged01ForCiC,"pid_pfIsoCharged01ForCiC[40][100]/F");
  m_tree->Branch("pid_pfIsoCharged02ForCiC",&pid_pfIsoCharged02ForCiC,"pid_pfIsoCharged02ForCiC[40][100]/F");
  m_tree->Branch("pid_pfIsoCharged03ForCiC",&pid_pfIsoCharged03ForCiC,"pid_pfIsoCharged03ForCiC[40][100]/F");
  m_tree->Branch("pid_pfIsoCharged04ForCiC",&pid_pfIsoCharged04ForCiC,"pid_pfIsoCharged04ForCiC[40][100]/F");
  m_tree->Branch("pid_pfIsoCharged05ForCiC",&pid_pfIsoCharged05ForCiC,"pid_pfIsoCharged05ForCiC[40][100]/F");
  m_tree->Branch("pid_pfIsoCharged06ForCiC",&pid_pfIsoCharged06ForCiC,"pid_pfIsoCharged06ForCiC[40][100]/F");

  m_tree->Branch("pid_pfIsoPhotons01ForCiC",&pid_pfIsoPhotons01ForCiC,"pid_pfIsoPhotons01ForCiC[nPhot]/F");
  m_tree->Branch("pid_pfIsoPhotons02ForCiC",&pid_pfIsoPhotons02ForCiC,"pid_pfIsoPhotons02ForCiC[nPhot]/F");
  m_tree->Branch("pid_pfIsoPhotons03ForCiC",&pid_pfIsoPhotons03ForCiC,"pid_pfIsoPhotons03ForCiC[nPhot]/F");
  m_tree->Branch("pid_pfIsoPhotons04ForCiC",&pid_pfIsoPhotons04ForCiC,"pid_pfIsoPhotons04ForCiC[nPhot]/F");
  m_tree->Branch("pid_pfIsoPhotons05ForCiC",&pid_pfIsoPhotons05ForCiC,"pid_pfIsoPhotons05ForCiC[nPhot]/F");
  m_tree->Branch("pid_pfIsoPhotons06ForCiC",&pid_pfIsoPhotons06ForCiC,"pid_pfIsoPhotons06ForCiC[nPhot]/F");

  m_tree->Branch("pid_pfIsoNeutrals01ForCiC",&pid_pfIsoNeutrals01ForCiC,"pid_pfIsoNeutrals01ForCiC[nPhot]/F");
  m_tree->Branch("pid_pfIsoNeutrals02ForCiC",&pid_pfIsoNeutrals02ForCiC,"pid_pfIsoNeutrals02ForCiC[nPhot]/F");
  m_tree->Branch("pid_pfIsoNeutrals03ForCiC",&pid_pfIsoNeutrals03ForCiC,"pid_pfIsoNeutrals03ForCiC[nPhot]/F");
  m_tree->Branch("pid_pfIsoNeutrals04ForCiC",&pid_pfIsoNeutrals04ForCiC,"pid_pfIsoNeutrals04ForCiC[nPhot]/F");
  m_tree->Branch("pid_pfIsoNeutrals05ForCiC",&pid_pfIsoNeutrals05ForCiC,"pid_pfIsoNeutrals05ForCiC[nPhot]/F");
  m_tree->Branch("pid_pfIsoNeutrals06ForCiC",&pid_pfIsoNeutrals06ForCiC,"pid_pfIsoNeutrals06ForCiC[nPhot]/F");

  m_tree->Branch("pid_pfIsoFPRCharged02",&pid_pfIsoFPRCharged02,"pid_pfIsoFPRCharged02[nPhot]/F");
  m_tree->Branch("pid_pfIsoFPRNeutral02",&pid_pfIsoFPRNeutral02,"pid_pfIsoFPRNeutral02[nPhot]/F");
  m_tree->Branch("pid_pfIsoFPRPhoton02",&pid_pfIsoFPRPhoton02,"pid_pfIsoFPRPhoton02[nPhot]/F");
  m_tree->Branch("pid_pfIsoFPRRandomConeCharged02",&pid_pfIsoFPRRandomConeCharged02,"pid_pfIsoFPRRandomConeCharged02[nPhot]/F");
  m_tree->Branch("pid_pfIsoFPRRandomConeNeutral02",&pid_pfIsoFPRRandomConeNeutral02,"pid_pfIsoFPRRandomConeNeutral02[nPhot]/F");
  m_tree->Branch("pid_pfIsoFPRRandomConePhoton02",&pid_pfIsoFPRRandomConePhoton02,"pid_pfIsoFPRRandomConePhoton02[nPhot]/F");
  m_tree->Branch("pid_pfIsoFPRRandomConeEta02",&pid_pfIsoFPRRandomConeEta02,"pid_pfIsoFPRRandomConeEta02[nPhot]/F");
  m_tree->Branch("pid_pfIsoFPRRandomConePhi02",&pid_pfIsoFPRRandomConePhi02,"pid_pfIsoFPRRandomConePhi02[nPhot]/F");

  m_tree->Branch("pid_pfIsoFPRCharged03",&pid_pfIsoFPRCharged03,"pid_pfIsoFPRCharged03[nPhot]/F");
  m_tree->Branch("pid_pfIsoFPRNeutral03",&pid_pfIsoFPRNeutral03,"pid_pfIsoFPRNeutral03[nPhot]/F");
  m_tree->Branch("pid_pfIsoFPRPhoton03",&pid_pfIsoFPRPhoton03,"pid_pfIsoFPRPhoton03[nPhot]/F");
  m_tree->Branch("pid_pfIsoFPRRandomConeCharged03",&pid_pfIsoFPRRandomConeCharged03,"pid_pfIsoFPRRandomConeCharged03[nPhot]/F");
  m_tree->Branch("pid_pfIsoFPRRandomConeNeutral03",&pid_pfIsoFPRRandomConeNeutral03,"pid_pfIsoFPRRandomConeNeutral03[nPhot]/F");
  m_tree->Branch("pid_pfIsoFPRRandomConePhoton03",&pid_pfIsoFPRRandomConePhoton03,"pid_pfIsoFPRRandomConePhoton03[nPhot]/F");
  m_tree->Branch("pid_pfIsoFPRRandomConeEta03",&pid_pfIsoFPRRandomConeEta03,"pid_pfIsoFPRRandomConeEta03[nPhot]/F");
  m_tree->Branch("pid_pfIsoFPRRandomConePhi03",&pid_pfIsoFPRRandomConePhi03,"pid_pfIsoFPRRandomConePhi03[nPhot]/F");

  m_tree->Branch("pid_pfIsoFPRCharged04",&pid_pfIsoFPRCharged04,"pid_pfIsoFPRCharged04[nPhot]/F");
  m_tree->Branch("pid_pfIsoFPRNeutral04",&pid_pfIsoFPRNeutral04,"pid_pfIsoFPRNeutral04[nPhot]/F");
  m_tree->Branch("pid_pfIsoFPRPhoton04",&pid_pfIsoFPRPhoton04,"pid_pfIsoFPRPhoton04[nPhot]/F");
  m_tree->Branch("pid_pfIsoFPRRandomConeCharged04",&pid_pfIsoFPRRandomConeCharged04,"pid_pfIsoFPRRandomConeCharged04[nPhot]/F");
  m_tree->Branch("pid_pfIsoFPRRandomConeNeutral04",&pid_pfIsoFPRRandomConeNeutral04,"pid_pfIsoFPRRandomConeNeutral04[nPhot]/F");
  m_tree->Branch("pid_pfIsoFPRRandomConePhoton04",&pid_pfIsoFPRRandomConePhoton04,"pid_pfIsoFPRRandomConePhoton04[nPhot]/F");
  m_tree->Branch("pid_pfIsoFPRRandomConeEta04",&pid_pfIsoFPRRandomConeEta04,"pid_pfIsoFPRRandomConeEta04[nPhot]/F");
  m_tree->Branch("pid_pfIsoFPRRandomConePhi04",&pid_pfIsoFPRRandomConePhi04,"pid_pfIsoFPRRandomConePhi04[nPhot]/F");

  m_tree->Branch("ptiso004Phot",&ptiso004Phot,"ptiso004Phot[nPhot]/F");
  m_tree->Branch("ntrkiso004Phot",&ntrkiso004Phot,"ntrkiso004Phot[nPhot]/I");
  m_tree->Branch("ptiso035Phot",&ptiso035Phot,"ptiso035Phot[nPhot]/F");
  m_tree->Branch("ntrkiso035Phot",&ntrkiso035Phot,"ntrkiso035Phot[nPhot]/I");
  m_tree->Branch("ptiso04Phot",&ptiso04Phot,"ptiso04Phot[nPhot]/F");
  m_tree->Branch("ntrkiso04Phot",&ntrkiso04Phot,"ntrkiso04Phot[nPhot]/I");

  m_tree->Branch("hcalovecal04Phot",&hcalovecal04Phot,"hcalovecal04Phot[nPhot]/F"); 
  m_tree->Branch("ecaliso04Phot",&ecaliso04Phot,"ecaliso04Phot[nPhot]/F");  

  m_tree->Branch("pid_scetawid",&pid_scetawid,"pid_scetawid[nPhot]/F");
  m_tree->Branch("pid_scphiwid",&pid_scphiwid,"pid_scphiwid[nPhot]/F");
  m_tree->Branch("pid_lambdaRatio",&pid_lambdaRatio,"pid_lambdaRatio[nPhot]/F");
  m_tree->Branch("pid_esXwidth",&pid_esXwidth,"pid_esXwidth[nPhot]/F");
  m_tree->Branch("pid_esYwidth",&pid_esYwidth,"pid_esYwidth[nPhot]/F");
  m_tree->Branch("pid_esXshape",&pid_esXShape,"pid_esXshape[40][21]/F"); 
  m_tree->Branch("pid_esYshape",&pid_esYShape,"pid_esYshape[40][21]/F"); 

  m_tree->Branch("sMajMajPhot",&sMajMajPhot,"sMajMaj2Phot[nPhot]/F");
  m_tree->Branch("sMinMinPhot",&sMinMinPhot,"sMinMin2Phot[nPhot]/F");
  m_tree->Branch("alphaPhot",&alphaPhot,"alphaPhot[nPhot]/F");
  m_tree->Branch("sEtaEtaPhot",&sEtaEtaPhot,"sEtaEtaPhot[nPhot]/F");
  m_tree->Branch("sEtaPhiPhot",&sEtaPhiPhot,"sEtaPhiPhot[nPhot]/F");
  m_tree->Branch("sPhiPhiPhot",&sPhiPhiPhot,"sPhiPhiPhot[nPhot]/F");
  m_tree->Branch("E1Phot",&E1Phot,"E1Phot[nPhot]/F");
  m_tree->Branch("E2OverE9Phot",&E2OverE9Phot,"E2OverE9Phot[nPhot]/F");
  m_tree->Branch("E4Phot",&E4Phot,"E4Phot[nPhot]/F");
  m_tree->Branch("E9Phot",&E9Phot,"E9Phot[nPhot]/F");
  m_tree->Branch("E25Phot",&E25Phot,"E25Phot[nPhot]/F");
  m_tree->Branch("ieleassocPhot",&ieleassocPhot,"ieleassocPhot[nPhot]/I");
  m_tree->Branch("pid_deltaRToTrackPhot",&pid_deltaRToTrackPhot,"pid_deltaRToTrackPhot[nPhot]/F");

  m_tree->Branch("nElePhot",&nElePhot,"nElePhot/I");
  m_tree->Branch("pid_jurECALElePhot ",&pid_jurECALElePhot ,"pid_jurECALElePhot[nElePhot]/F");
  m_tree->Branch("pid_twrHCALElePhot ",&pid_twrHCALElePhot ,"pid_twrHCALElePhot[nElePhot]/F");
  m_tree->Branch("pid_HoverEElePhot ",&pid_HoverEElePhot ,"pid_HoverEElePhot[nElePhot]/F");
  m_tree->Branch("pid_hlwTrackElePhot ",&pid_hlwTrackElePhot ,"pid_hlwTrackElePhot[nElePhot]/F");
  m_tree->Branch("pid_etawidElePhot ",&pid_etawidElePhot ,"pid_etawidElePhot[nElePhot]/F");
  m_tree->Branch("pid_dphivtxElePhot ",&pid_dphivtxElePhot ,"pid_dphivtxElePhot[nElePhot]/F");
  m_tree->Branch("pid_detavtxElePhot ",&pid_detavtxElePhot ,"pid_detavtxElePhot[nElePhot]/F");
  m_tree->Branch("pid_mishitsElePhot ",&pid_mishitsElePhot ,"pid_mishitsElePhot[nElePhot]/I");
  m_tree->Branch("pid_distElePhot ",&pid_distElePhot ,"pid_distElePhot[nElePhot]/F");
  m_tree->Branch("pid_dcotElePhot ",&pid_dcotElePhot ,"pid_dcotElePhot[nElePhot]/F");
  m_tree->Branch("pid_ptElePhot ",&pid_ptElePhot ,"pid_ptElePhot[nElePhot]/F");

  if (dumpAKT5Jets_)
    {
      m_tree->Branch("nJet_akt5",&nJet_akt5,"nJet_akt5/I");
      m_tree->Branch("ptJet_akt5 ",&ptJet_akt5 ,"ptJet_akt5[nJet_akt5]/F");
      m_tree->Branch("ptCorrJet_akt5 ",&ptCorrJet_akt5 ,"ptCorrJet_akt5[nJet_akt5]/F");
      m_tree->Branch("eJet_akt5  ",&eJet_akt5  ,"eJet_akt5[nJet_akt5]/F");
      m_tree->Branch("etaJet_akt5",&etaJet_akt5,"etaJet_akt5[nJet_akt5]/F");
      m_tree->Branch("phiJet_akt5",&phiJet_akt5,"phiJet_akt5[nJet_akt5]/F");
      m_tree->Branch("emfJet_akt5",&emfJet_akt5,"emfJet_akt5[nJet_akt5]/F");
      m_tree->Branch("n90Jet_akt5",&n90Jet_akt5,"n90Jet_akt5[nJet_akt5]/F");
      m_tree->Branch("n90HitsJet_akt5",&n90HitsJet_akt5,"n90HitsJet_akt5[nJet_akt5]/F");
      m_tree->Branch("fHPDJet_akt5",&fHPDJet_akt5,"fHPDJet_akt5[nJet_akt5]/F");
      m_tree->Branch("fRBXJet_akt5",&fRBXJet_akt5,"fRBXJet_akt5[nJet_akt5]/F");
    }

  if (dumpAKT7Jets_)
    {
      m_tree->Branch("nJet_akt7",&nJet_akt7,"nJet_akt7/I");
      m_tree->Branch("ptJet_akt7 ",&ptJet_akt7 ,"ptJet_akt7[nJet_akt7]/F");
      m_tree->Branch("ptCorrJet_akt7 ",&ptCorrJet_akt7 ,"ptCorrJet_akt7[nJet_akt5]/F");
      m_tree->Branch("eJet_akt7  ",&eJet_akt7  ,"eJet_akt7[nJet_akt7]/F");
      m_tree->Branch("etaJet_akt7",&etaJet_akt7,"etaJet_akt7[nJet_akt7]/F");
      m_tree->Branch("phiJet_akt7",&phiJet_akt7,"phiJet_akt7[nJet_akt7]/F");
      m_tree->Branch("emfJet_akt7",&emfJet_akt7,"emfJet_akt7[nJet_akt7]/F");
      m_tree->Branch("n90Jet_akt7",&n90Jet_akt7,"n90Jet_akt7[nJet_akt7]/F");
      m_tree->Branch("n90HitsJet_akt7",&n90HitsJet_akt7,"n90HitsJet_akt7[nJet_akt7]/F");
      m_tree->Branch("fHPDJet_akt7",&fHPDJet_akt7,"fHPDJet_akt7[nJet_akt7]/F");
      m_tree->Branch("fRBXJet_akt7",&fRBXJet_akt7,"fRBXJet_akt7[nJet_akt7]/F");
    }

  if (dumpJPTAKT5Jets_)
    {
      m_tree->Branch("nJet_jptak5",&nJet_jptak5,"nJet_jptak5/I");
      m_tree->Branch("ptJet_jptak5 ",&ptJet_jptak5 ,"ptJet_jptak5[nJet_jptak5]/F");
      m_tree->Branch("ptCorrJet_jptak5 ",&ptCorrJet_jptak5 ,"ptCorrJet_jptak5[nJet_jptak5]/F");
      m_tree->Branch("eJet_jptak5  ",&eJet_jptak5  ,"eJet_jptak5[nJet_jptak5]/F");
      m_tree->Branch("etaJet_jptak5",&etaJet_jptak5,"etaJet_jptak5[nJet_jptak5]/F");
      m_tree->Branch("phiJet_jptak5",&phiJet_jptak5,"phiJet_jptak5[nJet_jptak5]/F");
      m_tree->Branch("emfJet_jptak5",&emfJet_jptak5,"emfJet_jptak5[nJet_jptak5]/F");
    }

  if (dumpKT4Jets_)
    {
      m_tree->Branch("nJet_pfkt4",&nJet_pfkt4,"nJet_pfkt4/I");
      m_tree->Branch("ptJet_pfkt4 ",&ptJet_pfkt4 ,"ptJet_pfkt4[nJet_pfkt4]/F");
      m_tree->Branch("eJet_pfkt4  ",&eJet_pfkt4  ,"eJet_pfkt4[nJet_pfkt4]/F");
      m_tree->Branch("etaJet_pfkt4",&etaJet_pfkt4,"etaJet_pfkt4[nJet_pfkt4]/F");
      m_tree->Branch("phiJet_pfkt4",&phiJet_pfkt4,"phiJet_pfkt4[nJet_pfkt4]/F");
    }

  if (dumpPFAKT5Jets_)
    {
      m_tree->Branch("nJet_pfakt5",&nJet_pfakt5,"nJet_pfakt5/I");
      m_tree->Branch("ptJet_pfakt5 ",&ptJet_pfakt5 ,"ptJet_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("ptCorrJet_pfakt5 ",&ptCorrJet_pfakt5 ,"ptCorrJet_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("eJet_pfakt5  ",&eJet_pfakt5  ,"eJet_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("etaJet_pfakt5",&etaJet_pfakt5,"etaJet_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("phiJet_pfakt5",&phiJet_pfakt5,"phiJet_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("ptDJet_pfakt5",&ptDJet_pfakt5,"ptDJet_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("rmsCandJet_pfakt5",&rmsCandJet_pfakt5,"rmsCandJet_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("rmsCandTrueJet_pfakt5",&rmsCandTrueJet_pfakt5,"rmsCandTrueJet_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("rmsCandTrue_QCJet_pfakt5",&rmsCandTrue_QCJet_pfakt5,"rmsCandTrue_QCJet_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("axis1Jet_pfakt5",&axis1Jet_pfakt5,"axis1Jet_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("axis2Jet_pfakt5",&axis2Jet_pfakt5,"axis2Jet_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("pullJet_pfakt5",&pullJet_pfakt5,"pullJet_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("tanaJet_pfakt5",&tanaJet_pfakt5,"tanaJet_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("ptD_QCJet_pfakt5",&ptD_QCJet_pfakt5,"ptD_QCJet_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("rmsCandTrue_QCJet_pfakt5",&rmsCandTrue_QCJet_pfakt5,"rmsCandTrue_QCJet_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("axis1_QCJet_pfakt5",&axis1_QCJet_pfakt5,"axis1_QCJet_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("axis2_QCJet_pfakt5",&axis2_QCJet_pfakt5,"axis2_QCJet_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("pull_QCJet_pfakt5",&pull_QCJet_pfakt5,"pull_QCJet_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("tana_QCJet_pfakt5",&tana_QCJet_pfakt5,"tana_QCJet_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("RchgJet_pfakt5",&RchgJet_pfakt5,"RchgJet_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("RneutralJet_pfakt5",&RneutralJet_pfakt5,"RneutralJet_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("RJet_pfakt5",&RJet_pfakt5,"RJet_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("Rchg_QCJet_pfakt5",&Rchg_QCJet_pfakt5,"Rchg_QCJet_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("nChg_ptCutJet_pfakt5",&nChg_ptCutJet_pfakt5,"nChg_ptCutJet_pfakt5[nJet_pfakt5]/I");
      m_tree->Branch("nChg_QCJet_pfakt5",&nChg_QCJet_pfakt5,"nChg_QCJet_pfakt5[nJet_pfakt5]/I");
      m_tree->Branch("nChg_ptCut_QCJet_pfakt5",&nChg_ptCut_QCJet_pfakt5,"nChg_ptCut_QCJet_pfakt5[nJet_pfakt5]/I");
      m_tree->Branch("nNeutral_ptCutJet_pfakt5",&nNeutral_ptCutJet_pfakt5,"nNeutral_ptCutJet_pfakt5[nJet_pfakt5]/I");
      m_tree->Branch("pTMaxJet_pfakt5",&pTMaxJet_pfakt5,"pTMaxJet_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("pTMaxChgJet_pfakt5",&pTMaxChgJet_pfakt5,"pTMaxChgJet_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("pTMaxNeutralJet_pfakt5",&pTMaxNeutralJet_pfakt5,"pTMaxNeutralJet_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("pTMaxChg_QCJet_pfakt5",&pTMaxChg_QCJet_pfakt5,"pTMaxChg_QCJet_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("jetId_dRMean_pfakt5",&jetId_dRMean_pfakt5,"jetId_dRMean_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("jetId_frac01_pfakt5",&jetId_frac01_pfakt5,"jetId_frac01_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("jetId_frac02_pfakt5",&jetId_frac02_pfakt5,"jetId_frac02_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("jetId_frac03_pfakt5",&jetId_frac03_pfakt5,"jetId_frac03_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("jetId_frac04_pfakt5",&jetId_frac04_pfakt5,"jetId_frac04_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("jetId_frac05_pfakt5",&jetId_frac05_pfakt5,"jetId_frac05_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("jetId_nNeutrals_pfakt5",&jetId_nNeutrals_pfakt5,"jetId_nNeutrals_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("jetId_beta_pfakt5",&jetId_beta_pfakt5,"jetId_beta_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("jetId_betaStar_pfakt5",&jetId_betaStar_pfakt5,"jetId_betaStar_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("jetId_dZ_pfakt5",&jetId_dZ_pfakt5,"jetId_dZ_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("jetId_nCharged_pfakt5",&jetId_nCharged_pfakt5,"jetId_nCharged_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("jetId_dR2Mean_pfakt5",&jetId_dR2Mean_pfakt5,"jetId_dR2Mean_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("jetId_betaStarClassic_pfakt5",&jetId_betaStarClassic_pfakt5,"jetId_betaStarClassic_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("jetIdSimple_mva_pfakt5",&jetIdSimple_mva_pfakt5,"jetIdSimple_mva_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("jetIdSimple_wp_pfakt5",&jetIdSimple_wp_pfakt5,"jetIdSimple_wp_pfakt5[nJet_pfakt5]/I");
      m_tree->Branch("jetIdFull_mva_pfakt5",&jetIdFull_mva_pfakt5,"jetIdFull_mva_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("jetIdFull_wp_pfakt5",&jetIdFull_wp_pfakt5,"jetIdFull_wp_pfakt5[nJet_pfakt5]/I");
      m_tree->Branch("jetIdCutBased_mva_pfakt5",&jetIdCutBased_mva_pfakt5,"jetIdCutBased_mva_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("jetIdCutBased_wp_pfakt5",&jetIdCutBased_wp_pfakt5,"jetIdCutBased_wp_pfakt5[nJet_pfakt5]/I");

      //Calculated for each vertex
      m_tree->Branch("beta_pfakt5",&beta_pfakt5,"beta_pfakt5[100][100]/F");
      m_tree->Branch("betaStar_pfakt5",&betaStar_pfakt5,"betaStar_pfakt5[100][100]/F");

      m_tree->Branch("combinedSecondaryVertexBJetTags", &combinedSecondaryVertexBJetTags_pfakt5, "combinedSecondaryVertexBJetTags[nJet_pfakt5]/F");
      m_tree->Branch("combinedSecondaryVertexMVABJetTags", &combinedSecondaryVertexMVABJetTags_pfakt5, "combinedSecondaryVertexMVABJetTags[nJet_pfakt5]/F");
      m_tree->Branch("jetBProbabilityBJetTags", &jetBProbabilityBJetTags_pfakt5, "jetBProbabilityBJetTags[nJet_pfakt5]/F");
      m_tree->Branch("jetProbabilityBJetTags", &jetProbabilityBJetTags_pfakt5, "jetProbabilityBJetTags[nJet_pfakt5]/F");
      m_tree->Branch("simpleSecondaryVertexHighEffBJetTags", &simpleSecondaryVertexHighEffBJetTags_pfakt5, "simpleSecondaryVertexHighEffBJetTags[nJet_pfakt5]/F");
      m_tree->Branch("simpleSecondaryVertexHighPurBJetTags", &simpleSecondaryVertexHighPurBJetTags_pfakt5, "simpleSecondaryVertexHighPurBJetTags[nJet_pfakt5]/F");
      m_tree->Branch("softMuonBJetTags", &softMuonBJetTags_pfakt5, "softMuonBJetTags[nJet_pfakt5]/F");
      m_tree->Branch("softMuonByIP3dBJetTags", &softMuonByIP3dBJetTags_pfakt5, "softMuonByIP3dBJetTags[nJet_pfakt5]/F");
      m_tree->Branch("softMuonByPtBJetTags", &softMuonByPtBJetTags_pfakt5, "softMuonByPtBJetTags[nJet_pfakt5]/F");
      m_tree->Branch("softElectronByIP3dBJetTags", &softElectronByIP3dBJetTags_pfakt5, "softElectronByIP3dBJetTags[nJet_pfakt5]/F");
      m_tree->Branch("softElectronByPtBJetTags", &softElectronByPtBJetTags_pfakt5, "softElectronByPtBJetTags[nJet_pfakt5]/F");
      m_tree->Branch("trackCountingHighPurBJetTags", &trackCountingHighPurBJetTags_pfakt5, "trackCountingHighPurBJetTags[nJet_pfakt5]/F");
      m_tree->Branch("trackCountingHighEffBJetTags", &trackCountingHighEffBJetTags_pfakt5, "trackCountingHighEffBJetTags[nJet_pfakt5]/F");
      
      m_tree->Branch("npfcand_all",&npfcand_all,"npfcand_all/I");
      m_tree->Branch("nChargedHadrons_uncl",&nChargedHadrons_uncl,"nChargedHadrons_uncl/I");
      m_tree->Branch("nChargedHadronsgoodvtx_uncl",&nChargedHadronsgoodvtx_uncl,"nChargedHadronsgoodvtx_uncl/I");
      m_tree->Branch("nChargedHadronsnoothervtx_uncl",&nChargedHadronsnoothervtx_uncl,"nChargedHadronsothervtx_uncl/I");
      m_tree->Branch("nPhotons_uncl",       &nPhotons_uncl,       "nPhotons_uncl/I");
      m_tree->Branch("nMuons_uncl",         &nMuons_uncl,         "nMuons_uncl/I");
      m_tree->Branch("nElectrons_uncl",     &nElectrons_uncl,     "nElectrons_uncl/I");
      m_tree->Branch("nNeutralHadrons_uncl",&nNeutralHadrons_uncl,"nNeutralHadrons_uncl/I");
      m_tree->Branch("nHFHadrons_uncl",     &nHFHadrons_uncl,     "nHFHadrons_uncl/I");
      m_tree->Branch("nHFEM_uncl",     &nHFEM_uncl,     "nHFEM_uncl/I");
      
      m_tree->Branch("epfcand_all",&epfcand_all,"epfcand_all/F");
      m_tree->Branch("eChargedHadrons_uncl",&eChargedHadrons_uncl,"eChargedHadrons_uncl/F");
      m_tree->Branch("eChargedHadronsgoodvtx_uncl",&eChargedHadronsgoodvtx_uncl,"eChargedHadronsgoodvtx_uncl/F");
      m_tree->Branch("eChargedHadronsnoothervtx_uncl",&eChargedHadronsnoothervtx_uncl,"eChargedHadronsnoothervtx_uncl/F");
      m_tree->Branch("ePhotons_uncl",&ePhotons_uncl,"ePhotons_uncl/F");
      m_tree->Branch("eMuons_uncl",&eMuons_uncl,"eMuons_uncl/F");
      m_tree->Branch("eElectrons_uncl",&eElectrons_uncl,"eElectrons_uncl/F");
      m_tree->Branch("eNeutralHadrons_uncl",&eNeutralHadrons_uncl,"eNeutralHadrons_uncl/F");
      m_tree->Branch("eHFHadrons_uncl",&eHFHadrons_uncl,"eHFHadrons_uncl/F");
      m_tree->Branch("eHFEM_uncl",&eHFEM_uncl,"eHFEM_uncl/F");
      m_tree->Branch("ptpfcand_all",&ptpfcand_all,"ptpfcand_all/F");
      m_tree->Branch("ptChargedHadrons_uncl",&ptChargedHadrons_uncl,"ptChargedHadrons_uncl/F");
      m_tree->Branch("ptChargedHadronsgoodvtx_uncl",&ptChargedHadronsgoodvtx_uncl,"ptChargedHadronsgoodvtx_uncl/F");
      m_tree->Branch("ptChargedHadronsnoothervtx_uncl",&ptChargedHadronsnoothervtx_uncl,"ptChargedHadronsnoothervtx_uncl/F");
      m_tree->Branch("ptPhotons_uncl",&ptPhotons_uncl,"ptPhotons_uncl/F");
      m_tree->Branch("ptMuons_uncl",&ptMuons_uncl,"ptMuons_uncl/F");
      m_tree->Branch("ptElectrons_uncl",&ptElectrons_uncl,"ptElectrons_uncl/F");
      m_tree->Branch("ptNeutralHadrons_uncl",&ptNeutralHadrons_uncl,"ptNeutralHadrons_uncl/F");
      m_tree->Branch("ptHFHadrons_uncl",&ptHFHadrons_uncl,"ptHFHadrons_uncl/F");
      m_tree->Branch("ptHFEM_uncl",&ptHFEM_uncl,"ptHFEM_uncl/F");
      m_tree->Branch("ptpfcand_all",&ptpfcand_all,"ptpfcand_all/F");
      m_tree->Branch("etaChargedHadrons_uncl",&etaChargedHadrons_uncl,"etaChargedHadrons_uncl/F");
      m_tree->Branch("etaChargedHadronsgoodvtx_uncl",&etaChargedHadronsgoodvtx_uncl,"etaChargedHadronsgoodvtx_uncl/F");
      m_tree->Branch("etaChargedHadronsnoothervtx_uncl",&etaChargedHadronsnoothervtx_uncl,"etaChargedHadronsnoothervtx_uncl/F");
      m_tree->Branch("etaPhotons_uncl",&etaPhotons_uncl,"etaPhotons_uncl/F");
      m_tree->Branch("etaMuons_uncl",&etaMuons_uncl,"etaMuons_uncl/F");
      m_tree->Branch("etaElectrons_uncl",&etaElectrons_uncl,"etaElectrons_uncl/F");
      m_tree->Branch("etaNeutralHadrons_uncl",&etaNeutralHadrons_uncl,"etaNeutralHadrons_uncl/F");
      m_tree->Branch("etaHFHadrons_uncl",&etaHFHadrons_uncl,"etaHFHadrons_uncl/F");
      m_tree->Branch("etaHFEM_uncl",&etaHFEM_uncl,"etaHFEM_uncl/F");
      m_tree->Branch("ptpfcand_all",&ptpfcand_all,"ptpfcand_all/F");
      m_tree->Branch("phiChargedHadrons_uncl",&phiChargedHadrons_uncl,"phiChargedHadrons_uncl/F");
      m_tree->Branch("phiChargedHadronsgoodvtx_uncl",&phiChargedHadronsgoodvtx_uncl,"phiChargedHadronsgoodvtx_uncl/F");
      m_tree->Branch("phiChargedHadronsnoothervtx_uncl",&phiChargedHadronsnoothervtx_uncl,"phiChargedHadronsnoothervtx_uncl/F");
      m_tree->Branch("phiPhotons_uncl",&phiPhotons_uncl,"phiPhotons_uncl/F");
      m_tree->Branch("phiMuons_uncl",&phiMuons_uncl,"phiMuons_uncl/F");
      m_tree->Branch("phiElectrons_uncl",&phiElectrons_uncl,"phiElectrons_uncl/F");
      m_tree->Branch("phiNeutralHadrons_uncl",&phiNeutralHadrons_uncl,"phiNeutralHadrons_uncl/F");
      m_tree->Branch("phiHFHadrons_uncl",&phiHFHadrons_uncl,"phiHFHadrons_uncl/F");
      m_tree->Branch("phiHFEM_uncl",&phiHFEM_uncl,"phiHFEM_uncl/F");
      m_tree->Branch("sumptpfcand_all",&sumptpfcand_all,"sumptpfcand_all/F");
      m_tree->Branch("sumptChargedHadrons_uncl",&sumptChargedHadrons_uncl,"sumptChargedHadrons_uncl/F");
      m_tree->Branch("sumptChargedHadronsgoodvtx_uncl",&sumptChargedHadronsgoodvtx_uncl,"sumptChargedHadronsgoodvtx_uncl/F");
      m_tree->Branch("sumptChargedHadronsnoothervtx_uncl",&sumptChargedHadronsnoothervtx_uncl,"sumptChargedHadronsnoothervtx_uncl/F");
      m_tree->Branch("sumptPhotons_uncl",&sumptPhotons_uncl,"sumptPhotons_uncl/F");
      m_tree->Branch("sumptMuons_uncl",&sumptMuons_uncl,"sumptMuons_uncl/F");
      m_tree->Branch("sumptElectrons_uncl",&sumptElectrons_uncl,"sumptElectrons_uncl/F");
      m_tree->Branch("sumptNeutralHadrons_uncl",&sumptNeutralHadrons_uncl,"sumptNeutralHadrons_uncl/F");
      m_tree->Branch("sumptHFHadrons_uncl",&sumptHFHadrons_uncl,"sumptHFHadrons_uncl/F");
      m_tree->Branch("sumptHFEM_uncl",&sumptHFEM_uncl,"sumptHFEM_uncl/F");
      m_tree->Branch("sumptpfcand_all",&sumptpfcand_all,"sumptpfcand_all/F");

      //   // Extra variables for PFlow studies
      m_tree->Branch("nChargedHadrons_pfakt5",nChargedHadrons_pfakt5,"nChargedHadrons_pfakt5[nJet_pfakt5]/I");
      m_tree->Branch("nChargedHadronsgoodvtx_pfakt5",nChargedHadronsgoodvtx_pfakt5,"nChargedHadronsgoodvtx_pfakt5[nJet_pfakt5]/I");
      m_tree->Branch("nChargedHadronsnoothervtx_pfakt5",nChargedHadronsnoothervtx_pfakt5,"nChargedHadronsnoothervtx_pfakt5[nJet_pfakt5]/I");
      m_tree->Branch("nPhotons_pfakt5",       nPhotons_pfakt5,       "nPhotons_pfakt5[nJet_pfakt5]/I");
      m_tree->Branch("nMuons_pfakt5",         nMuons_pfakt5,         "nMuons_pfakt5[nJet_pfakt5]/I");
      m_tree->Branch("nElectrons_pfakt5",     nElectrons_pfakt5,     "nElectrons_pfakt5[nJet_pfakt5]/I");
      m_tree->Branch("nNeutralHadrons_pfakt5",nNeutralHadrons_pfakt5,"nNeutralHadrons_pfakt5[nJet_pfakt5]/I");
      m_tree->Branch("nHFHadrons_pfakt5",     nHFHadrons_pfakt5,     "nHFHadrons_pfakt5[nJet_pfakt5]/I");
      m_tree->Branch("nHFEM_pfakt5",     nHFEM_pfakt5,     "nHFEM_pfakt5[nJet_pfakt5]/I");
      
      m_tree->Branch("eChargedHadrons_pfakt5",eChargedHadrons_pfakt5,"eChargedHadrons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("eChargedHadronsgoodvtx_pfakt5",eChargedHadronsgoodvtx_pfakt5,"eChargedHadronsgoodvtx_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("eChargedHadronsnoothervtx_pfakt5",eChargedHadronsnoothervtx_pfakt5,"eChargedHadronsnoothervtx_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("ePhotons_pfakt5",ePhotons_pfakt5,"ePhotons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("eMuons_pfakt5",eMuons_pfakt5,"eMuons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("eElectrons_pfakt5",eElectrons_pfakt5,"eElectrons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("eNeutralHadrons_pfakt5",eNeutralHadrons_pfakt5,"eNeutralHadrons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("eHFHadrons_pfakt5",eHFHadrons_pfakt5,"eHFHadrons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("eHFEM_pfakt5",eHFEM_pfakt5,"eHFEM_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("ptChargedHadrons_pfakt5",ptChargedHadrons_pfakt5,"ptChargedHadrons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("ptChargedHadronsgoodvtx_pfakt5",ptChargedHadronsgoodvtx_pfakt5,"ptChargedHadronsgoodvtx_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("ptChargedHadronsnoothervtx_pfakt5",ptChargedHadronsnoothervtx_pfakt5,"ptChargedHadronsnoothervtx_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("ptPhotons_pfakt5",ptPhotons_pfakt5,"ptPhotons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("ptMuons_pfakt5",ptMuons_pfakt5,"ptMuons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("ptElectrons_pfakt5",ptElectrons_pfakt5,"ptElectrons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("ptNeutralHadrons_pfakt5",ptNeutralHadrons_pfakt5,"ptNeutralHadrons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("ptHFHadrons_pfakt5",ptHFHadrons_pfakt5,"ptHFHadrons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("ptHFEM_pfakt5",ptHFEM_pfakt5,"ptHFEM_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("etaChargedHadrons_pfakt5",etaChargedHadrons_pfakt5,"etaChargedHadrons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("etaChargedHadronsgoodvtx_pfakt5",etaChargedHadronsgoodvtx_pfakt5,"etaChargedHadronsgoodvtx_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("etaChargedHadronsnoothervtx_pfakt5",etaChargedHadronsnoothervtx_pfakt5,"etaChargedHadronsnoothervtx_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("etaPhotons_pfakt5",etaPhotons_pfakt5,"etaPhotons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("etaMuons_pfakt5",etaMuons_pfakt5,"etaMuons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("etaElectrons_pfakt5",etaElectrons_pfakt5,"etaElectrons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("etaNeutralHadrons_pfakt5",etaNeutralHadrons_pfakt5,"etaNeutralHadrons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("etaHFHadrons_pfakt5",etaHFHadrons_pfakt5,"etaHFHadrons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("etaHFEM_pfakt5",etaHFEM_pfakt5,"etaHFEM_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("phiChargedHadrons_pfakt5",phiChargedHadrons_pfakt5,"phiChargedHadrons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("phiChargedHadronsgoodvtx_pfakt5",phiChargedHadronsgoodvtx_pfakt5,"phiChargedHadronsgoodvtx_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("phiChargedHadronsnoothervtx_pfakt5",phiChargedHadronsnoothervtx_pfakt5,"phiChargedHadronsnoothervtx_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("phiPhotons_pfakt5",phiPhotons_pfakt5,"phiPhotons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("phiMuons_pfakt5",phiMuons_pfakt5,"phiMuons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("phiElectrons_pfakt5",phiElectrons_pfakt5,"phiElectrons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("phiNeutralHadrons_pfakt5",phiNeutralHadrons_pfakt5,"phiNeutralHadrons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("phiHFHadrons_pfakt5",phiHFHadrons_pfakt5,"phiHFHadrons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("phiHFEM_pfakt5",phiHFEM_pfakt5,"phiHFEM_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("sumptChargedHadrons_pfakt5",sumptChargedHadrons_pfakt5,"sumptChargedHadrons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("sumptChargedHadronsgoodvtx_pfakt5",sumptChargedHadronsgoodvtx_pfakt5,"sumptChargedHadronsgoodvtx_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("sumptChargedHadronsnoothervtx_pfakt5",sumptChargedHadronsnoothervtx_pfakt5,"sumptChargedHadronsnoothervtx_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("sumptPhotons_pfakt5",sumptPhotons_pfakt5,"sumptPhotons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("sumptMuons_pfakt5",sumptMuons_pfakt5,"sumptMuons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("sumptElectrons_pfakt5",sumptElectrons_pfakt5,"sumptElectrons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("sumptNeutralHadrons_pfakt5",sumptNeutralHadrons_pfakt5,"sumptNeutralHadrons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("sumptHFHadrons_pfakt5",sumptHFHadrons_pfakt5,"sumptHFHadrons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("sumptHFEM_pfakt5",sumptHFEM_pfakt5,"sumptHFEM_pfakt5[nJet_pfakt5]/F");
    }

  if (dumpPFAKT7Jets_)
    {
      m_tree->Branch("nJet_pfakt7",&nJet_pfakt7,"nJet_pfakt7/I");
      m_tree->Branch("ptJet_pfakt7 ",&ptJet_pfakt7 ,"ptJet_pfakt7[nJet_pfakt7]/F");
      m_tree->Branch("ptCorrJet_pfakt7 ",&ptCorrJet_pfakt7 ,"ptCorrJet_pfakt7[nJet_pfakt7]/F");
      m_tree->Branch("eJet_pfakt7  ",&eJet_pfakt7  ,"eJet_pfakt7[nJet_pfakt7]/F");
      m_tree->Branch("etaJet_pfakt7",&etaJet_pfakt7,"etaJet_pfakt7[nJet_pfakt7]/F");
      m_tree->Branch("phiJet_pfakt7",&phiJet_pfakt7,"phiJet_pfakt7[nJet_pfakt7]/F");
    }

  if (dumpKT6Jets_)
    {
      m_tree->Branch("nJet_pfkt6",&nJet_pfkt6,"nJet_pfkt6/I");
      m_tree->Branch("ptJet_pfkt6 ",&ptJet_pfkt6 ,"ptJet_pfkt6[nJet_pfkt6]/F");
      m_tree->Branch("eJet_pfkt6  ",&eJet_pfkt6  ,"eJet_pfkt6[nJet_pfkt6]/F");
      m_tree->Branch("etaJet_pfkt6",&etaJet_pfkt6,"etaJet_pfkt6[nJet_pfkt6]/F");
      m_tree->Branch("phiJet_pfkt6",&phiJet_pfkt6,"phiJet_pfkt6[nJet_pfkt6]/F");
    }

  if (dumpAKT5Jets_ || dumpJPTAKT5Jets_ || dumpPFAKT5Jets_)
    {
      m_tree->Branch("nJetGen_akt5",&nJetGen_akt5,"nJetGen_akt5/I");
      m_tree->Branch("ptJetGen_akt5 ",&ptJetGen_akt5,
		     "ptJetGen_akt5[nJetGen_akt5]/F");
      m_tree->Branch("eJetGen_akt5  ",&eJetGen_akt5,
		     "eJetGen_akt5[nJetGen_akt5]/F");
      m_tree->Branch("etaJetGen_akt5",&etaJetGen_akt5,
		     "etaJetGen_akt5[nJetGen_akt5]/F");
      m_tree->Branch("phiJetGen_akt5",&phiJetGen_akt5,
		     "phiJetGen_akt5[nJetGen_akt5]/F");

      //   // Extra variables for PFlow studies
      m_tree->Branch("nMuonsGen_akt5", nMuonsGen_akt5,
		     "nMuonsGen_akt5[nJetGen_akt5]/I");
      m_tree->Branch("nElectronsGen_akt5", nElectronsGen_akt5,
		     "nElectronsGen_akt5[nJetGen_akt5]/I");
      m_tree->Branch("nPhotonsGen_akt5", nPhotonsGen_akt5,
		     "nPhotonsGen_akt5[nJetGen_akt5]/I");
      m_tree->Branch("nTracksGen_akt5", nTracksGen_akt5,
		     "nTracksGen_akt5[nJetGen_akt5]/I");
      m_tree->Branch("nNeutralHadronsGen_akt5", nNeutralHadronsGen_akt5,
		     "nNeutralHadronsGen_akt5[nJetGen_akt5]/I");
      m_tree->Branch("nHFHadronsGen_akt5", nHFHadronsGen_akt5,
		     "nHFHadronsGen_akt5[nJetGen_akt5]/I");
      m_tree->Branch("nHFEMGen_akt5", nHFEMGen_akt5,
		     "nHFEMGen_akt5[nJetGen_akt5]/I");
      
      m_tree->Branch("nNeutronsGen_akt5", nNeutronsGen_akt5,
		     "nNeutronsGen_akt5[nJetGen_akt5]/I");
      m_tree->Branch("nK0LGen_akt5", nK0LGen_akt5,
		     "nK0LGen_akt5[nJetGen_akt5]/I");
      m_tree->Branch("nK0SGen_akt5", nK0SGen_akt5,
		     "nK0SGen_akt5[nJetGen_akt5]/I");
      m_tree->Branch("nLambdasGen_akt5", nLambdasGen_akt5,
		     "nLambdasGen_akt5[nJetGen_akt5]/I");
      m_tree->Branch("nCsiGen_akt5", nCsiGen_akt5,
		     "nCsiGen_akt5[nJetGen_akt5]/I");
      m_tree->Branch("nOtherNeutralHadronsGen_akt5", nOtherNeutralHadronsGen_akt5,
		     "nOtherNeutralHadronsGen_akt5[nJetGen_akt5]/I");
      
      m_tree->Branch("eMuonsGen_akt5", eMuonsGen_akt5,
		     "eMuonsGen_akt5[nJetGen_akt5]/F");
      m_tree->Branch("eElectronsGen_akt5", eElectronsGen_akt5,
		     "eElectronsGen_akt5[nJetGen_akt5]/F");
      m_tree->Branch("ePhotonsGen_akt5", ePhotonsGen_akt5,
		     "ePhotonsGen_akt5[nJetGen_akt5]/F");
      m_tree->Branch("eTracksGen_akt5", eTracksGen_akt5,
		     "eTracksGen_akt5[nJetGen_akt5]/F");
      m_tree->Branch("eNeutralHadronsGen_akt5", eNeutralHadronsGen_akt5,
		     "eNeutralHadronsGen_akt5[nJetGen_akt5]/F");
      m_tree->Branch("eHFHadronsGen_akt5", eHFHadronsGen_akt5,
		     "eHFHadronsGen_akt5[nJetGen_akt5]/F");
      m_tree->Branch("eHFEMGen_akt5", eHFEMGen_akt5,
		     "eHFEMGen_akt5[nJetGen_akt5]/F");
      
      m_tree->Branch("eNeutronsGen_akt5", eNeutronsGen_akt5,
		     "eNeutronsGen_akt5[nJetGen_akt5]/F");
      m_tree->Branch("eK0LGen_akt5", eK0LGen_akt5,
		     "eK0LGen_akt5[nJetGen_akt5]/F");
      m_tree->Branch("eK0SGen_akt5", eK0SGen_akt5,
		     "eK0SGen_akt5[nJetGen_akt5]/F");
      m_tree->Branch("eLambdasGen_akt5", eLambdasGen_akt5,
		     "eLambdasGen_akt5[nJetGen_akt5]/F");
      m_tree->Branch("eCsiGen_akt5", eCsiGen_akt5,
		     "eCsiGen_akt5[nJetGen_akt5]/F");
      m_tree->Branch("eOtherNeutralHadronsGen_akt5", eOtherNeutralHadronsGen_akt5,
		     "eOtherNeutralHadronsGen_akt5[nJetGen_akt5]/F");
    }

  if (dumpAKT7Jets_ || dumpJPTAKT7Jets_ || dumpPFAKT7Jets_)
    {
      m_tree->Branch("nJetGen_akt7",&nJetGen_akt7,"nJetGen_akt7/I");
      m_tree->Branch("ptJetGen_akt7 ",&ptJetGen_akt7 ,"ptJetGen_akt7[nJetGen_akt7]/F");
      m_tree->Branch("eJetGen_akt7  ",&eJetGen_akt7  ,"eJetGen_akt7[nJetGen_akt7]/F");
      m_tree->Branch("etaJetGen_akt7",&etaJetGen_akt7,"etaJetGen_akt7[nJetGen_akt7]/F");
      m_tree->Branch("phiJetGen_akt7",&phiJetGen_akt7,"phiJetGen_akt7[nJetGen_akt7]/F");
    }
  if (dumpKT4Jets_)
    {
      m_tree->Branch("nJetGen_kt4",&nJetGen_kt4,"nJetGen_kt4/I");
      m_tree->Branch("ptJetGen_kt4 ",&ptJetGen_kt4 ,"ptJetGen_kt4[nJetGen_kt4]/F");
      m_tree->Branch("eJetGen_kt4  ",&eJetGen_kt4  ,"eJetGen_kt4[nJetGen_kt4]/F");
      m_tree->Branch("etaJetGen_kt4",&etaJetGen_kt4,"etaJetGen_kt4[nJetGen_kt4]/F");
      m_tree->Branch("phiJetGen_kt4",&phiJetGen_kt4,"phiJetGen_kt4[nJetGen_kt4]/F");
    }

  if (dumpKT6Jets_)
    {
      m_tree->Branch("nJetGen_kt6",&nJetGen_kt6,"nJetGen_kt6/I");
      m_tree->Branch("ptJetGen_kt6 ",&ptJetGen_kt6 ,"ptJetGen_kt6[nJetGen_kt6]/F");
      m_tree->Branch("eJetGen_kt6  ",&eJetGen_kt6  ,"eJetGen_kt6[nJetGen_kt6]/F");
      m_tree->Branch("etaJetGen_kt6",&etaJetGen_kt6,"etaJetGen_kt6[nJetGen_kt6]/F");
      m_tree->Branch("phiJetGen_kt6",&phiJetGen_kt6,"phiJetGen_kt6[nJetGen_kt6]/F");
    }


  //vertex info
  m_tree->Branch("nvertex",&nvertex,"nvertex/I");

  // pz,eta always zero, px,py redundant
  m_tree->Branch("sMet  ",&sMet  ,"sMet/F");
  m_tree->Branch("eMet  ",&eMet  ,"eMet/F");
  m_tree->Branch("phiMet",&phiMet,"phiMet/F");
  m_tree->Branch("signifMet",&signifMet,"signifMet/F");

  m_tree->Branch("sCorrMet  ",&sCorrMet  ,"sCorrMet/F");
  m_tree->Branch("eCorrMet  ",&eCorrMet  ,"eCorrMet/F");
  m_tree->Branch("phiCorrMet",&phiCorrMet,"phiCorrMet/F");
  m_tree->Branch("signifCorrMet",&signifCorrMet,"signifCorrMet/F");

  // new for 52X
  // m_tree->Branch("smuCorrMet  ",&smuCorrMet  ,"smuCorrMet/F");
  // m_tree->Branch("emuCorrMet  ",&emuCorrMet  ,"emuCorrMet/F");
  // m_tree->Branch("phimuCorrMet",&phimuCorrMet,"phimuCorrMet/F");
  // m_tree->Branch("signifmuCorrMet",&signifmuCorrMet,"signifmuCorrMet/F");

  m_tree->Branch("sNoHFMet  ",&sNoHFMet  ,"sNoHFMet/F");
  m_tree->Branch("eNoHFMet  ",&eNoHFMet  ,"eNoHFMet/F");
  m_tree->Branch("phiNoHFMet",&phiNoHFMet,"phiNoHFMet/F");
  m_tree->Branch("signifNoHFMet",&signifNoHFMet,"signifNoHFMet/F");

  m_tree->Branch("stcMet",&stcMet  ,"stcMet/F");
  m_tree->Branch("etcMet",&etcMet  ,"etcMet/F");
  m_tree->Branch("phitcMet",&phitcMet,"phitcMet/F");
  m_tree->Branch("signiftcMet",&signiftcMet,"signiftcMet/F");

  m_tree->Branch("sglobalPfMet",&sglobalPfMet  ,"sglobalPfMet/F");
  m_tree->Branch("eglobalPfMet",&eglobalPfMet  ,"eglobalPfMet/F");
  m_tree->Branch("phiglobalPfMet",&phiglobalPfMet,"phiglobalPfMet/F");
  m_tree->Branch("signifglobalPfMet",&signifglobalPfMet,"signifglobalPfMet/F");

  m_tree->Branch("scentralPfMet",&scentralPfMet  ,"scentralPfMet/F");
  m_tree->Branch("ecentralPfMet",&ecentralPfMet  ,"ecentralPfMet/F");
  m_tree->Branch("phicentralPfMet",&phicentralPfMet,"phicentralPfMet/F");
  m_tree->Branch("signifcentralPfMet",&signifcentralPfMet,"signifcentralPfMet/F");

  m_tree->Branch("eassocPfMet",&eassocPfMet  ,"eassocPfMet[nvertex]/F");
  m_tree->Branch("phiassocPfMet",&phiassocPfMet,"phiassocPfMet[nvertex]/F");
  m_tree->Branch("signifassocPfMet",&signifassocPfMet,"signifassocPfMet[nvertex]/F");

  m_tree->Branch("eassocOtherVtxPfMet",&eassocOtherVtxPfMet  ,"eassocOtherVtxPfMet[nvertex]/F");
  m_tree->Branch("phiassocOtherVtxPfMet",&phiassocOtherVtxPfMet,"phiassocOtherVtxPfMet[nvertex]/F");
  m_tree->Branch("signifassocOtherVtxPfMet",&signifassocOtherVtxPfMet,"signifassocOtherVtxPfMet[nvertex]/F");

  m_tree->Branch("etrkPfMet",&etrkPfMet  ,"etrkPfMet[nvertex]/F");
  m_tree->Branch("phitrkPfMet",&phitrkPfMet,"phitrkPfMet[nvertex]/F");
  m_tree->Branch("signiftrkPfMet",&signiftrkPfMet,"signiftrkPfMet[nvertex]/F");

  m_tree->Branch("ecleanPfMet",&ecleanPfMet  ,"ecleanPfMet[nvertex]/F");
  m_tree->Branch("phicleanPfMet",&phicleanPfMet,"phicleanPfMet[nvertex]/F");
  m_tree->Branch("signifcleanPfMet",&signifcleanPfMet,"signifcleanPfMet[nvertex]/F");

  // new for 52X
  // m_tree->Branch("ecleanedSaclayPfMet",&ecleanedSaclayPfMet  ,"ecleanedSaclayPfMet[nvertex]/F");
  // m_tree->Branch("phicleanedSaclayPfMet",&phicleanedSaclayPfMet,"phicleanedSaclayPfMet[nvertex]/F");
  // m_tree->Branch("signifcleanedSaclayPfMet",&signifcleanedSaclayPfMet,"signifcleanedSaclayPfMet[nvertex]/F");
  // m_tree->Branch("eminTypeICleanSaclayPfMet",&eminTypeICleanSaclayPfMet  ,"eminTypeICleanSaclayPfMet[nvertex]/F");
  // m_tree->Branch("phiminTypeICleanSaclayPfMet",&phiminTypeICleanSaclayPfMet,"phiminTypeICleanSaclayPfMet[nvertex]/F");
  // m_tree->Branch("signifminTypeICleanSaclayPfMet",&signifminTypeICleanSaclayPfMet,"signifminTypeICleanSaclayPfMet[nvertex]/F");

  m_tree->Branch("globalPfSums",&globalPfSums  ,"globalPfSums[12]/F");

  m_tree->Branch("spfMet",&spfMet  ,"spfMet/F");
  m_tree->Branch("epfMet",&epfMet  ,"epfMet/F");
  m_tree->Branch("phipfMet",&phipfMet,"phipfMet/F");
  m_tree->Branch("signifpfMet",&signifpfMet,"signifpfMet/F");

  m_tree->Branch("spfMetType1",&spfMetType1  ,"spfMetType1/F");
  m_tree->Branch("epfMetType1",&epfMetType1  ,"epfMetType1/F");
  m_tree->Branch("phipfMetType1",&phipfMetType1,"phipfMetType1/F");
  m_tree->Branch("signifpfMetType1",&signifpfMetType1,"signifpfMetType1/F");

  m_tree->Branch("sMetGen",&sMetGen  ,"sMetGen/F");
  m_tree->Branch("eMetGen",&eMetGen  ,"eMetGen/F");
  m_tree->Branch("phiMetGen",&phiMetGen,"phiMetGen/F");
  m_tree->Branch("signifMetGen",&signifMetGen,"signifMetGen/F");

  m_tree->Branch("sMetGen2",&sMetGen2  ,"sMetGen2/F");
  m_tree->Branch("eMetGen2",&eMetGen2  ,"eMetGen2/F");
  m_tree->Branch("phiMetGen2",&phiMetGen2,"phiMetGen2/F");



  m_tree->Branch("vxMC",&vxMC,"vxMC/F");
  m_tree->Branch("vyMC",&vyMC,"vyMC/F");
  m_tree->Branch("vzMC",&vzMC,"vzMC/F");

  m_tree->Branch("vx",&vx,"vx[nvertex]/F");
  m_tree->Branch("vy",&vy,"vy[nvertex]/F");
  m_tree->Branch("vz",&vz,"vz[nvertex]/F");
  m_tree->Branch("vntracks",&vntracks,"vntracks[nvertex]/F");
  m_tree->Branch("vchi2",&vchi2,"vchi2[nvertex]/F");
  m_tree->Branch("vndof",&vndof,"vndof[nvertex]/F");
  m_tree->Branch("vlogsumpt2",&vlogsumpt2,"vlogsumpt2[nvertex]/F");

  m_tree->Branch("nPreselPhotonPairs",&nPreselPhotonPairs,"nPreselPhotonPairs/I");
  m_tree->Branch("indexPreselPhot1",&indexPreselPhot1,"indexPreselPhot1[nPreselPhotonPairs]/I");
  m_tree->Branch("indexPreselPhot2",&indexPreselPhot2,"indexPreselPhot2[nPreselPhotonPairs]/I");
  m_tree->Branch("vrankPhotonPairs",&vrankPhotonPairs,"vrankPhotonPairs[nPreselPhotonPairs]/I");
  m_tree->Branch("vevtMvaPhotonPairs",&vevtMvaPhotonPairs,"vevtMvaPhotonPairs[nPreselPhotonPairs]/F");
  m_tree->Branch("vevtProbPhotonPairs",&vevtProbPhotonPairs,"vevtProbPhotonPairs[nPreselPhotonPairs]/F");
  m_tree->Branch("vptbalPhotonPairs",&vptbalPhotonPairs,"vptbalPhotonPairs[nPreselPhotonPairs]/F");
  m_tree->Branch("vptasymPhotonPairs",&vptasymPhotonPairs,"vptasymPhotonPairs[nPreselPhotonPairs]/F");

  // Set trigger bits of interest
  nHLT = 13;
  hltTriggers["HLT_Photon10_L1R"] = 0;
  hltTriggers["HLT_Photon10_LooseEcalIso_TrackIso_L1R"] = 1;
  hltTriggers["HLT_Photon15_L1R"] = 2;
  hltTriggers["HLT_Photon20_L1R"] = 3;
  hltTriggers["HLT_Photon20_LooseEcalIso_TrackIso_L1R"] = 4;
  hltTriggers["HLT_Photon25_L1R"] = 5;
  hltTriggers["HLT_Photon25_LooseEcalIso_TrackIso_L1R"] = 6;
  hltTriggers["HLT_Photon30_L1R_1E31"] = 7;
  hltTriggers["HLT_Photon70_L1R"] = 8;

  hltTriggers["HLT_Photon10_Cleaned_L1R"] = 9;
  hltTriggers["HLT_Photon15_Cleaned_L1R"] = 10;
  hltTriggers["HLT_Photon20_Cleaned_L1R"] = 11;
  hltTriggers["HLT_Photon30_Cleaned_L1R"] = 12;


  //m_tree->Branch("hltPass",&hltPass,"hltPass/O");
  //m_tree->Branch("hltCount",&hltCount,"hltCount/I");
  m_tree->Branch("nHLT",&nHLT,"nHLT/I");
  m_tree->Branch("hltNamesLen",&hltNamesLen,"hltNamesLen/I");
  //m_tree->Branch("HLTNames",&aHLTNames,"HLTNames[hltNamesLen]/C,6000");
  m_tree->Branch("HLTNames",&aHLTNames);
  m_tree->Branch("HLTResults",&aHLTResults);
  //m_tree->Branch("HLTResults",&aHLTResults,"HLTResults[nHLT]/O");



  //Trigger objects for trigger Matching
  m_tree->Branch("trg17_SC_ele_n", &ElectronRefs0_n,"ElectronRefs0_n/I");
  m_tree->Branch("trg17_SC_ele_eta", &ElectronRefs0_eta,"ElectronRefs0_eta[ElectronRefs0_n]/F");
  m_tree->Branch("trg17_SC_ele_et", &ElectronRefs0_et,"ElectronRefs0_et[ElectronRefs0_n]/F");
  m_tree->Branch("trg17_SC_ele_phi", &ElectronRefs0_phi,"ElectronRefs0_phi[ElectronRefs0_n]/F");
  m_tree->Branch("trg32_ele_n", &ElectronRefs1_n,"ElectronRefs1_n/I");
  m_tree->Branch("trg32_ele_eta", &ElectronRefs1_eta,"ElectronRefs1_eta[ElectronRefs1_n]/F");
  m_tree->Branch("trg32_ele_et", &ElectronRefs1_et,"ElectronRefs1_et[ElectronRefs1_n]/F");
  m_tree->Branch("trg32_ele_phi", &ElectronRefs1_phi,"ElectronRefs1_phi[ElectronRefs1_n]/F");


  //Ele17_Ele8 trigger objects
  m_tree->Branch("trg8_ele_n", &ElectronRefs2_n,"ElectronRefs2_n/I");
  m_tree->Branch("trg8_ele_eta", &ElectronRefs2_eta,"ElectronRefs2_eta[ElectronRefs2_n]/F");
  m_tree->Branch("trg8_ele_et", &ElectronRefs2_et,"ElectronRefs2_et[ElectronRefs2_n]/F");
  m_tree->Branch("trg8_ele_phi", &ElectronRefs2_phi,"ElectronRefs2_phi[ElectronRefs2_n]/F");
  m_tree->Branch("trg17_ele_n", &ElectronRefs3_n,"ElectronRefs3_n/I");
  m_tree->Branch("trg17_ele_eta", &ElectronRefs3_eta,"ElectronRefs3_eta[ElectronRefs3_n]/F");
  m_tree->Branch("trg17_ele_et", &ElectronRefs3_et,"ElectronRefs3_et[ElectronRefs3_n]/F");
  m_tree->Branch("trg17_ele_phi", &ElectronRefs3_phi,"ElectronRefs3_phi[ElectronRefs3_n]/F");

  //Ele17_Ele8_mass50 trigger objects
  m_tree->Branch("trg8_mass50_ele_n", &ElectronRefs4_n,"ElectronRefs4_n/I");
  m_tree->Branch("trg8_mass50_ele_eta", &ElectronRefs4_eta,"ElectronRefs4_eta[ElectronRefs4_n]/F");
  m_tree->Branch("trg8_mass50_ele_et", &ElectronRefs4_et,"ElectronRefs4_et[ElectronRefs4_n]/F");
  m_tree->Branch("trg8_mass50_ele_phi", &ElectronRefs4_phi,"ElectronRefs4_phi[ElectronRefs4_n]/F");
  m_tree->Branch("trg17_mass50_ele_n", &ElectronRefs5_n,"ElectronRefs5_n/I");
  m_tree->Branch("trg17_mass50_ele_eta", &ElectronRefs5_eta,"ElectronRefs5_eta[ElectronRefs5_n]/F");
  m_tree->Branch("trg17_mass50_ele_et", &ElectronRefs5_et,"ElectronRefs5_et[ElectronRefs5_n]/F");
  m_tree->Branch("trg17_mass50_ele_phi", &ElectronRefs5_phi,"ElectronRefs5_phi[ElectronRefs5_n]/F");


  //Ele20_SC4_mass50 trigger objects
  m_tree->Branch("trg4_mass50_SC_n", &ElectronRefs6_n,"ElectronRefs6_n/I");
  m_tree->Branch("trg4_mass50_SC_eta", &ElectronRefs6_eta,"ElectronRefs6_eta[ElectronRefs6_n]/F");
  m_tree->Branch("trg4_mass50_SC_et", &ElectronRefs6_et,"ElectronRefs6_et[ElectronRefs6_n]/F");
  m_tree->Branch("trg4_mass50_SC_phi", &ElectronRefs6_phi,"ElectronRefs6_phi[ElectronRefs6_n]/F");
  m_tree->Branch("trg20_mass50_ele_n", &ElectronRefs7_n,"ElectronRefs7_n/I");
  m_tree->Branch("trg20_mass50_ele_eta", &ElectronRefs7_eta,"ElectronRefs7_eta[ElectronRefs7_n]/F");
  m_tree->Branch("trg20_mass50_ele_et", &ElectronRefs7_et,"ElectronRefs7_et[ElectronRefs7_n]/F");
  m_tree->Branch("trg20_mass50_ele_phi", &ElectronRefs7_phi,"ElectronRefs7_phi[ElectronRefs7_n]/F");


  //PhotonID CaloVL trigger objects
  m_tree->Branch("trg20_phoIDCaloVL_n", &PhotonRefs0_n,"PhotonRefs0_n/I");
  m_tree->Branch("trg20_phoIDCaloVL_eta", &PhotonRefs0_eta,"PhotonRefs0_eta[PhotonRefs0_n]/F");
  m_tree->Branch("trg20_phoIDCaloVL_et", &PhotonRefs0_et,"PhotonRefs0_et[PhotonRefs0_n]/F");
  m_tree->Branch("trg20_phoIDCaloVL_phi", &PhotonRefs0_phi,"PhotonRefs0_phi[PhotonRefs0_n]/F");
  m_tree->Branch("trg30_phoIDCaloVL_n", &PhotonRefs1_n,"PhotonRefs1_n/I");
  m_tree->Branch("trg30_phoIDCaloVL_eta", &PhotonRefs1_eta,"PhotonRefs1_eta[PhotonRefs1_n]/F");
  m_tree->Branch("trg30_phoIDCaloVL_et", &PhotonRefs1_et,"PhotonRefs1_et[PhotonRefs1_n]/F");
  m_tree->Branch("trg30_phoIDCaloVL_phi", &PhotonRefs1_phi,"PhotonRefs1_phi[PhotonRefs1_n]/F");
  m_tree->Branch("trg50_phoIDCaloVL_n", &PhotonRefs2_n,"PhotonRefs2_n/I");
  m_tree->Branch("trg50_phoIDCaloVL_eta", &PhotonRefs2_eta,"PhotonRefs2_eta[PhotonRefs2_n]/F");
  m_tree->Branch("trg50_phoIDCaloVL_et", &PhotonRefs2_et,"PhotonRefs2_et[PhotonRefs2_n]/F");
  m_tree->Branch("trg50_phoIDCaloVL_phi", &PhotonRefs2_phi,"PhotonRefs2_phi[PhotonRefs2_n]/F");
  m_tree->Branch("trg75_phoIDCaloVL_n", &PhotonRefs3_n,"PhotonRefs3_n/I");
  m_tree->Branch("trg75_phoIDCaloVL_eta", &PhotonRefs3_eta,"PhotonRefs3_eta[PhotonRefs3_n]/F");
  m_tree->Branch("trg75_phoIDCaloVL_et", &PhotonRefs3_et,"PhotonRefs3_et[PhotonRefs3_n]/F");
  m_tree->Branch("trg75_phoIDCaloVL_phi", &PhotonRefs3_phi,"PhotonRefs3_phi[PhotonRefs3_n]/F");
  m_tree->Branch("trg90_phoIDCaloVL_n", &PhotonRefs4_n,"PhotonRefs4_n/I");
  m_tree->Branch("trg90_phoIDCaloVL_eta", &PhotonRefs4_eta,"PhotonRefs4_eta[PhotonRefs4_n]/F");
  m_tree->Branch("trg90_phoIDCaloVL_et", &PhotonRefs4_et,"PhotonRefs4_et[PhotonRefs4_n]/F");
  m_tree->Branch("trg90_phoIDCaloVL_phi", &PhotonRefs4_phi,"PhotonRefs4_phi[PhotonRefs4_n]/F");

  m_tree->Branch("nEle",&nEle,"nEle/I");
  m_tree->Branch("electron_px",&electron_px,"electron_px[nEle]/F");
  m_tree->Branch("electron_py",&electron_py,"electron_py[nEle]/F");
  m_tree->Branch("electron_pz",&electron_pz,"electron_pz[nEle]/F");
  m_tree->Branch("electron_vx",&electron_vx,"electron_vx[nEle]/F");
  m_tree->Branch("electron_vy",&electron_vy,"electron_vy[nEle]/F");
  m_tree->Branch("electron_vz",&electron_vz,"electron_vz[nEle]/F");
  m_tree->Branch("electron_pt",&electron_pt,"electron_pt[nEle]/F");
  m_tree->Branch("electron_eta",&electron_eta,"electron_eta[nEle]/F");
  m_tree->Branch("electron_phi",&electron_phi,"electron_phi[nEle]/F");
  m_tree->Branch("electron_energy",&electron_energy,"electron_energy[nEle]/F");
  m_tree->Branch("electron_ecalEnergy",&electron_ecalEnergy,"electron_ecalEnergy[nEle]/F");     // chiara
  m_tree->Branch("electron_trackPatVtx",&electron_trackPatVtx,"electron_trackPatVtx[nEle]/F");  // chiara
  m_tree->Branch("electron_charge",&electron_charge,"electron_charge[nEle]/F");
  m_tree->Branch("electron_fBrem",&electron_fBrem,"electron_fBrem[nEle]/F");
  m_tree->Branch("electron_dist",&electron_dist,"electron_dist[nEle]/F");
  m_tree->Branch("electron_dcot",&electron_dcot,"electron_dcot[nEle]/F");
  m_tree->Branch("electron_misHits",&electron_misHits,"electron_misHits[nEle]/I");
  m_tree->Branch("electron_matchedConv",&electron_matchedConv,"electron_matchedConv[nEle]/I");  // chiara
  m_tree->Branch("electron_seedType",&electron_seedType,"electron_seedType[nEle]/I");
  m_tree->Branch("electron_EoP",&electron_EoP,"electron_EoP[nEle]/F");
  m_tree->Branch("electron_OneOverEMinusOneOverP",&electron_OneOverEMinusOneOverP,"electron_OneOverEMinusOneOverP[nEle]/F");
  m_tree->Branch("electron_r9",&electron_r9,"electron_r9[nEle]/F");
  m_tree->Branch("electron_nSubClusters",&electron_nSubClusters,"electron_nSubClusters[nEle]/I");
  m_tree->Branch("electron_trkIso",&electron_trkIso,"electron_trkIso[nEle]/F");
  m_tree->Branch("electron_ecalIso",&electron_ecalIso,"electron_ecalIso[nEle]/F");
  m_tree->Branch("electron_hcalIso",&electron_hcalIso,"electron_hcalIso[nEle]/F");
  m_tree->Branch("electron_trkIso03",&electron_trkIso03,"electron_trkIso03[nEle]/F");
  m_tree->Branch("electron_ecalIso03",&electron_ecalIso03,"electron_ecalIso03[nEle]/F");
  m_tree->Branch("electron_hcalIso03",&electron_hcalIso03,"electron_hcalIso03[nEle]/F");
  m_tree->Branch("electron_SigmaIetaIeta",&electron_SigmaIetaIeta,"electron_SigmaIetaIeta[nEle]/F");
  m_tree->Branch("electron_SigmaIphiIphi",&electron_SigmaIphiIphi,"electron_SigmaIphiIphi[nEle]/F");
  m_tree->Branch("electron_dEtaIn",&electron_dEtaIn,"electron_dEtaIn[nEle]/F");
  m_tree->Branch("electron_dPhiIn",&electron_dPhiIn,"electron_dPhiIn[nEle]/F");
  m_tree->Branch("electron_HoE",&electron_HoE,"electron_HoE[nEle]/F");
  m_tree->Branch("electron_pFlowMVA",&electron_pFlowMVA,"electron_pFlowMVA[nEle]/F");
  m_tree->Branch("electron_sc_energy",&electron_sc_energy,"electron_sc_energy[nEle]/F");
  m_tree->Branch("electron_sc_eta",&electron_sc_eta,"electron_sc_eta[nEle]/F");
  m_tree->Branch("electron_sc_phi",&electron_sc_phi,"electron_sc_phi[nEle]/F");
  // chiara
  m_tree->Branch("electron_mvaNonTrig",&electron_mvaNonTrig,"electron_mvaNonTrig[nEle]/F");
  m_tree->Branch("electron_mvaTrig",   &electron_mvaTrig,   "electron_mvaTrig[nEle]/F");
  m_tree->Branch("electron_chHad03Iso",&electron_chHad03Iso,"electron_chHad03Iso[nEle]/F");
  m_tree->Branch("electron_nHad03Iso", &electron_nHad03Iso, "electron_nHad03Iso[nEle]/F");
  m_tree->Branch("electron_phot03Iso", &electron_phot03Iso, "electron_phot03Iso[nEle]/F");
  m_tree->Branch("electron_chHad04Iso",&electron_chHad04Iso,"electron_chHad04Iso[nEle]/F");
  m_tree->Branch("electron_nHad04Iso", &electron_nHad04Iso, "electron_nHad04Iso[nEle]/F");
  m_tree->Branch("electron_phot04Iso", &electron_phot04Iso, "electron_phot04Iso[nEle]/F");
  m_tree->Branch("electron_chHad05Iso",&electron_chHad05Iso,"electron_chHad05Iso[nEle]/F");
  m_tree->Branch("electron_nHad05Iso", &electron_nHad05Iso, "electron_nHad05Iso[nEle]/F");
  m_tree->Branch("electron_phot05Iso", &electron_phot05Iso, "electron_phot05Iso[nEle]/F");
  // chiara

  if (dumpBeamHaloInformations_)
    {
      m_tree->Branch("isBeamHaloGlobalLoosePass",&isBeamHaloGlobalLoosePass,"isBeamHaloGlobalLoosePass/O");
      m_tree->Branch("isBeamHaloGlobalTightPass",&isBeamHaloGlobalTightPass,"isBeamHaloGloablTightPass/O");
      m_tree->Branch("isBeamHaloHcalLoosePass",&isBeamHaloHcalLoosePass,"isBeamHaloHcalLoosePass/O");
      m_tree->Branch("isBeamHaloHcalTightPass",&isBeamHaloHcalTightPass,"isBeamHaloHcalTightPass/O");
      m_tree->Branch("isBeamHaloCSCLoosePass",&isBeamHaloCSCLoosePass,"isBeamHaloCSCLoosePass/O");
      m_tree->Branch("isBeamHaloCSCTightPass",&isBeamHaloCSCTightPass,"isBeamHaloCSCTightPass/O");
      m_tree->Branch("isBeamHaloEcalLoosePass",&isBeamHaloEcalLoosePass,"isBeamHaloEcalLoosePass/O");
      m_tree->Branch("isBeamHaloEcalTightPass",&isBeamHaloEcalTightPass,"isBeamHaloEcalTightPass/O");
      m_tree->Branch("isBeamHaloIDTightPass",&isBeamHaloIDTightPass,"isBeamHaloIDTightPass/O");
      m_tree->Branch("isBeamHaloIDLoosePass",&isBeamHaloIDLoosePass,"isBeamHaloIDLoosePass/O");
      m_tree->Branch("isSmellsLikeHalo_Tag",&isSmellsLikeHalo_Tag, "isSmellsLikeHalo_Tag/O");
      m_tree->Branch("isLooseHalo_Tag",&isLooseHalo_Tag, "isLooseHalo_Tag/O");
      m_tree->Branch("isTightHalo_Tag",&isTightHalo_Tag, "isTightHalo_Tag/O");
      m_tree->Branch("isExtremeTightHalo_Tag",&isExtremeTightHalo_Tag, "isExtremeTightHalo_Tag/O");

      m_tree->Branch("nMuons",&nMuonsReco,"nMuons/I");
      m_tree->Branch("Muon_px",muon_px,"muon_px[nMuons]/F");
      m_tree->Branch("Muon_py",muon_py,"muon_py[nMuons]/F");
      m_tree->Branch("Muon_pz",muon_pz,"muon_pz[nMuons]/F");
      m_tree->Branch("Muon_vx",muon_vx,"muon_vx[nMuons]/F");
      m_tree->Branch("Muon_vy",muon_vy,"muon_vy[nMuons]/F");
      m_tree->Branch("Muon_vz",muon_vz,"muon_vz[nMuons]/F");
      m_tree->Branch("Muon_pt",muon_pt,"muon_pt[nMuons]/F");
      m_tree->Branch("Muon_eta",muon_eta,"muon_eta[nMuons]/F");
      m_tree->Branch("Muon_phi",muon_phi,"muon_phi[nMuons]/F");
      m_tree->Branch("Muon_energy",muon_energy,"muon_energy[nMuons]/F");
      m_tree->Branch("Muon_charge",muon_charge,"muon_charge[nMuons]/F");
      m_tree->Branch("Muon_isGlobalMuon",muon_isGlobalMuon,"muon_isGlobalMuon[nMuons]/O");
      m_tree->Branch("Muon_isTrackerMuon",muon_isTrackerMuon,"muon_isTrackerMuon[nMuons]/O");
      m_tree->Branch("Muon_isStandAloneMuon",muon_isStandAloneMuon,"muon_isStandAloneMuon[nMuons]/O");
      m_tree->Branch("Muon_trkLayerWithMeas",muon_trkLayerWithMeas,"Muon_trkLayerWithMeas[nMuons]/I"); // chiara
      m_tree->Branch("Muon_InnerTrack_isNonnull",muon_InnerTrack_isNonnull,"muon_InnerTrack_isNonnull[nMuons]/O");
      m_tree->Branch("Muon_OuterTrack_isNonnull",muon_OuterTrack_isNonnull,"muon_OuterTrack_isNonnull[nMuons]/O");
      m_tree->Branch("Muon_OuterPoint_x",muon_OuterPoint_x,"muon_OuterPoint_x[nMuons]/F");
      m_tree->Branch("Muon_OuterPoint_y",muon_OuterPoint_y,"muon_OuterPoint_y[nMuons]/F");
      m_tree->Branch("Muon_OuterPoint_z",muon_OuterPoint_z,"muon_OuterPoint_z[nMuons]/F");
      // chiara
      m_tree->Branch("Muon_isPFMuon",muon_isPFMuon,"muon_isPFMuon[nMuons]/O");  
      m_tree->Branch("Muon_pfiso04_chHad",muon_pfiso04_chHad,"muon_pfiso04_chHad[nMuons]/F");
      m_tree->Branch("Muon_pfiso04_chPar",muon_pfiso04_chPar,"muon_pfiso04_chPar[nMuons]/F");
      m_tree->Branch("Muon_pfiso04_nHad", muon_pfiso04_nHad, "muon_pfiso04_nHad[nMuons]/F");
      m_tree->Branch("Muon_pfiso04_Phot", muon_pfiso04_Phot, "muon_pfiso04_Phot[nMuons]/F");
      m_tree->Branch("Muon_pfiso04_PUPt", muon_pfiso04_PUPt, "muon_pfiso04_PUPt[nMuons]/F");
      m_tree->Branch("Muon_pfiso03_chHad",muon_pfiso03_chHad,"muon_pfiso03_chHad[nMuons]/F");
      m_tree->Branch("Muon_pfiso03_chPar",muon_pfiso03_chPar,"muon_pfiso03_chPar[nMuons]/F");
      m_tree->Branch("Muon_pfiso03_nHad", muon_pfiso03_nHad, "muon_pfiso03_nHad[nMuons]/F");
      m_tree->Branch("Muon_pfiso03_Phot", muon_pfiso03_Phot, "muon_pfiso03_Phot[nMuons]/F");
      m_tree->Branch("Muon_pfiso03_PUPt", muon_pfiso03_PUPt, "muon_pfiso03_PUPt[nMuons]/F");
      // chiara
      //for Global,Tracker muon
      m_tree->Branch("Muon_InnerPoint_x",muon_InnerPoint_x,"muon_InnerPoint_x[nMuons]/F");
      m_tree->Branch("Muon_InnerPoint_y",muon_InnerPoint_y,"muon_InnerPoint_y[nMuons]/F");
      m_tree->Branch("Muon_InnerPoint_z",muon_InnerPoint_z,"muon_InnerPoint_z[nMuons]/F");
      m_tree->Branch("Muon_trackIso",muon_trackIso,"muon_trackIso[nMuons]/F");
      m_tree->Branch("Muon_ecalIso",muon_ecalIso,"muon_ecalIso[nMuons]/F");
      m_tree->Branch("Muon_hcalIso",muon_hcalIso,"muon_hcalIso[nMuons]/F");
      m_tree->Branch("Muon_relIso",muon_relIso,"muon_relIso[nMuons]/F");
      m_tree->Branch("Muon_normChi2",muon_normChi2,"muon_normChi2[nMuons]/I");
      m_tree->Branch("Muon_validHits",muon_validHits,"muon_validHits[nMuons]/I");
      m_tree->Branch("Muon_tkHits",muon_tkHits,"muon_tkHits[nMuons]/I");
      m_tree->Branch("Muon_pixHits",muon_pixHits,"muon_pixHits[nMuons]/I");
      m_tree->Branch("Muon_numberOfMatches",muon_numberOfMatches,"muon_numberOfMatches[nMuons]/I");

      m_tree->Branch("nCosmicMuons",&nCosmicMuons,"nCosmicMuons/I");
      m_tree->Branch("CosmicMuon_px",cosmicmuon_px,"cosmicmuon_px[nCosmicMuons]/F");
      m_tree->Branch("CosmicMuon_py",cosmicmuon_py,"cosmicmuon_py[nCosmicMuons]/F");
      m_tree->Branch("CosmicMuon_pz",cosmicmuon_pz,"cosmicmuon_pz[nCosmicMuons]/F");
      m_tree->Branch("CosmicMuon_pt",cosmicmuon_pt,"cosmicmuon_pt[nCosmicMuons]/F");
      m_tree->Branch("CosmicMuon_eta",cosmicmuon_eta,"cosmicmuon_eta[nCosmicMuons]/F");
      m_tree->Branch("CosmicMuon_phi",cosmicmuon_phi,"cosmicmuon_phi[nCosmicMuons]/F");
      m_tree->Branch("CosmicMuon_energy",cosmicmuon_energy,"cosmicmuon_energy[nCosmicMuons]/F");
      m_tree->Branch("CosmicMuon_charge",cosmicmuon_charge,"cosmicmuon_charge[nCosmicMuons]/F");
      m_tree->Branch("CosmicMuon_isGlobalMuon",cosmicmuon_isGlobalMuon,"cosmicmuon_isGlobalMuon[nCosmicMuons]/O");
      m_tree->Branch("CosmicMuon_isTrackerMuon",cosmicmuon_isTrackerMuon,"cosmicmuon_isTrackerMuon[nCosmicMuons]/O");
      m_tree->Branch("CosmicMuon_isStandAloneMuon",cosmicmuon_isStandAloneMuon,"cosmicmuon_isStandAloneMuon[nCosmicMuons]/O");
      m_tree->Branch("CosmicMuon_InnerTrack_isNonnull",cosmicmuon_InnerTrack_isNonnull,"cosmicmuon_InnerTrack_isNonnull[nCosmicMuons]/O");
      m_tree->Branch("CosmicMuon_OuterTrack_isNonnull",cosmicmuon_OuterTrack_isNonnull,"cosmicmuon_OuterTrack_isNonnull[nCosmicMuons]/O");
      m_tree->Branch("CosmicMuon_OuterPoint_x",cosmicmuon_OuterPoint_x,"cosmicmuon_OuterPoint_x[nCosmicMuons]/F");
      m_tree->Branch("CosmicMuon_OuterPoint_y",cosmicmuon_OuterPoint_y,"cosmicmuon_OuterPoint_y[nCosmicMuons]/F");
      m_tree->Branch("CosmicMuon_OuterPoint_z",cosmicmuon_OuterPoint_z,"cosmicmuon_OuterPoint_z[nCosmicMuons]/F");
    }


  m_tree->Branch("Xsec",  &Xsec_, "Xsec/D");


  event = 0;  
}

// ------------ method called once each job just after ending the event loop  ------------
void GammaJetAnalyzer::endJob() {

  
  outfile->cd("myanalysis");
  m_tree->Write();
  outfile->Close();
  delete jetID_;

  //outfile->Delete();
  
//   //avoid writing the tree second time (automatically)
//  m_tree->Delete();

}

// Method for iterative printing of decay chains
bool GammaJetAnalyzer::printChildren(const SimTrack* p, 
		    map<const SimTrack*, set<const SimTrack*> > const& ptokids,
		    map<const SimTrack*, const SimVertex*> const& ptovtx,
		    int level, bool save) {

  // Print parent
  bool hasvtx = (ptovtx.find(p) != ptovtx.end());
  if (_debug) {
    for (int i = 0; i != 2*level; ++i) cout << " "; // pad with spaces
    cout << Form("* %d (%1.3g GeV, %1.3g",
		 p->type(), p->momentum().pt(),p->momentum().eta());
    if (hasvtx)
      cout << Form(" => r %1.3g cm, z %1.3g cm)",
		   ptovtx.find(p)->second->position().Rho(),
		   ptovtx.find(p)->second->position().z()) << endl;
    else
      cout << ")" << endl;
  }
  
  bool hasKids = (ptokids.find(p) != ptokids.end());

  // Save only SIM tracks not already in GenParticle list
  bool saved = false;
  if (save && level > 0) {
    saved = true;
  }

  // Print children, if any
  if (hasKids) {

    set<const SimTrack*> const& kids = ptokids.find(p)->second;
    for (set<const SimTrack*>::const_iterator iKid = kids.begin();
	 iKid != kids.end(); ++iKid)
      saved |= printChildren(*iKid, ptokids, ptovtx, level+1, save);
  } // if kids

  return saved;
 } // printChildren

// Go down in chain and remove unwanted decays
bool GammaJetAnalyzer::pruneKids(const SimTrack* p,
     map<const SimTrack*, set<const SimTrack*> > & decays,
     map<const SimTrack*, const SimTrack*> & parent,
     map<const SimTrack*, const SimVertex*> & vertex,
     int level) {

  // No children, go one level back
  if (decays.find(p)==decays.end()) return false;
  
  // Prune kids and see if there are any grandchildren left after pruning
  set<const SimTrack*> const& kids = decays.find(p)->second;
  bool hasGrandKids = false;
  bool hasSameType = false;
  unsigned int nPhotons = 0;
  unsigned int nElectrons = 0;
  for (set<const SimTrack*>::const_iterator iKid = kids.begin();
       iKid != kids.end(); ++iKid) {

    bool hasKids = pruneKids(*iKid, decays, parent, vertex, level+1);
    hasGrandKids = hasGrandKids || hasKids;
    hasSameType = hasSameType || (*iKid)->type()==p->type(); 
    if ((*iKid)->type()==22) ++nPhotons;
    if (abs((*iKid)->type())==11) ++nElectrons;
  }
  // if there are grandkids, don't prune kids as we need the whole chain
  if (hasGrandKids) return true;

  // See if we have some reason to prune the kids
  double pt  = p->momentum().pt();
  bool prune = (hasSameType && nPhotons==kids.size()-1) // bremsstrahlung
    || (nElectrons==kids.size() && level!=0) // non-primary photon conversion
    || (abs(p->type())==11 && nPhotons==kids.size()) // weird brem (no e)
    || (abs(p->type())==11 && nPhotons==kids.size()-nElectrons) // ionization
    || (p->type()==111 && pt<2. && nPhotons==kids.size()) // low pT pi0
    || (p->type()==22 && pt<2. && nElectrons==kids.size()); // low pT conv
  // || (kids.size()==1); // not a real decay?
  // (NB: electron brems can produce more than one photon)
  // (NG: electrons can turn into a photon with much less pT?)
  // (NB: photon conversions can produce only one electron)
  // (NB: pizeros can produce only one photon with most of pT)
  // (NB: pi+ decays seem to only produce a muon, no neutrinos) 

  // Prune, i.e. remove the parent decay and all the kids from maps
  if (prune) {

    for (set<const SimTrack*>::const_iterator iKid = kids.begin();
	 iKid != kids.end(); ++iKid) {
      parent.erase(*iKid);
      vertex.erase(*iKid);
    } // for kids
    decays.erase(p);

    return false;
  } // if prune
  else // no pruning done
    return true; 

} // pruneKids

//define this as a plug-in
DEFINE_FWK_MODULE(GammaJetAnalyzer);
