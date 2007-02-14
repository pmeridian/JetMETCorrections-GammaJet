// user include files
#include "JetMETCorrections/GammaJet/interface/GammaJetAnalysis.h" 

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetfwd.h"

#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/ClusterShape.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "DataFormats/Common/interface/Provenance.h"
using namespace std;
namespace cms
{
GammaJetAnalysis::GammaJetAnalysis(const edm::ParameterSet& iConfig)
{
  // get names of modules, producing object collections
  hbheLabel_= iConfig.getParameter<edm::InputTag>("hbheInput");
  hoLabel_=iConfig.getParameter<edm::InputTag>("hoInput");
  hfLabel_=iConfig.getParameter<edm::InputTag>("hfInput");
  mInputCalo = iConfig.getParameter<std::vector<edm::InputTag> >("srcCalo");
  mInputGen = iConfig.getParameter<std::vector<edm::InputTag> >("srcGen");
  ecalLabels_=iConfig.getParameter<std::vector<edm::InputTag> >("ecalInputs");
  myName = iConfig.getParameter<std::string> ("textout");
  
  
// ECAL reconstruction

  islandBarrelBasicClusterCollection_ = iConfig.getParameter<std::string>("islandBarrelBasicClusterCollection");
  islandBarrelBasicClusterProducer_   = iConfig.getParameter<std::string>("islandBarrelBasicClusterProducer");
  islandBarrelBasicClusterShapes_   = iConfig.getParameter<std::string>("islandBarrelBasicClusterShapes");

  islandBarrelSuperClusterCollection_ = iConfig.getParameter<std::string>("islandBarrelSuperClusterCollection");
  islandBarrelSuperClusterProducer_   = iConfig.getParameter<std::string>("islandBarrelSuperClusterProducer");

  correctedIslandBarrelSuperClusterCollection_ = iConfig.getParameter<std::string>("correctedIslandBarrelSuperClusterCollection");
  correctedIslandBarrelSuperClusterProducer_   = iConfig.getParameter<std::string>("correctedIslandBarrelSuperClusterProducer");

  islandEndcapBasicClusterCollection_ = iConfig.getParameter<std::string>("islandEndcapBasicClusterCollection");
  islandEndcapBasicClusterProducer_   = iConfig.getParameter<std::string>("islandEndcapBasicClusterProducer");
  islandEndcapBasicClusterShapes_   = iConfig.getParameter<std::string>("islandEndcapBasicClusterShapes");

  islandEndcapSuperClusterCollection_ = iConfig.getParameter<std::string>("islandEndcapSuperClusterCollection");
  islandEndcapSuperClusterProducer_   = iConfig.getParameter<std::string>("islandEndcapSuperClusterProducer");

  correctedIslandEndcapSuperClusterCollection_ = iConfig.getParameter<std::string>("correctedIslandEndcapSuperClusterCollection");
  correctedIslandEndcapSuperClusterProducer_   = iConfig.getParameter<std::string>("correctedIslandEndcapSuperClusterProducer");

  hybridSuperClusterCollection_ = iConfig.getParameter<std::string>("hybridSuperClusterCollection");
  hybridSuperClusterProducer_   = iConfig.getParameter<std::string>("hybridSuperClusterProducer");

  correctedHybridSuperClusterCollection_ = iConfig.getParameter<std::string>("correctedHybridSuperClusterCollection");
  correctedHybridSuperClusterProducer_   = iConfig.getParameter<std::string>("correctedHybridSuperClusterProducer");

  
  allowMissingInputs_=iConfig.getUntrackedParameter<bool>("AllowMissingInputs",false);
  // get name of output file with histogramms
  fOutputFileName = iConfig.getUntrackedParameter<string>("HistOutFile"); 
  risol[0] = 0.5;
  risol[1] = 0.7;
  risol[2] = 1.0;
  
  ecut[0][0] = 0.09;
  ecut[0][1] = 0.18;
  ecut[0][2] = 0.27;
  
  ecut[1][0] = 0.45;
  ecut[1][1] = 0.9;
  ecut[1][2] = 1.35;
  
  ecut[2][0] = 0.5;
  ecut[2][1] = 1.;
  ecut[2][2] = 1.5;
}


GammaJetAnalysis::~GammaJetAnalysis()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

void GammaJetAnalysis::beginJob( const edm::EventSetup& iSetup)
{
   hOutputFile   = new TFile( fOutputFileName.c_str(), "RECREATE" ) ; 
   myTree = new TTree("GammaJet","GammaJet Tree");
   myTree->Branch("run",  &run, "run/I");
   myTree->Branch("event",  &event, "event/I");
   
// Jet block
   myTree->Branch("NumRecoJets", &NumRecoJets, "NumRecoJets/I");   
   myTree->Branch("JetRecoEt",  JetRecoEt, "JetRecoEt[10]/F");
   myTree->Branch("JetRecoEta",  JetRecoEta, "JetRecoEta[10]/F");
   myTree->Branch("JetRecoPhi",  JetRecoPhi, "JetRecoPhi[10]/F");
   myTree->Branch("JetRecoType",  JetRecoType, "JetRecoType[10]/F");
   
   myTree->Branch("NumGenJets", &NumGenJets, "NumGenJets/I");   
   myTree->Branch("JetGenEt",  JetGenEt, "JetGenEt[10]/F");
   myTree->Branch("JetGenEta",  JetGenEta, "JetGenEta[10]/F");
   myTree->Branch("JetGenPhi",  JetGenPhi, "JetGenPhi[10]/F");
   myTree->Branch("JetGenType",  JetGenType, "JetGenType[10]/F");
   
// Gamma block for ECAL isolated gammas
   myTree->Branch("NumRecoGamma", &NumRecoGamma, "NumRecoGamma/I");
   myTree->Branch("EcalClusDet", &EcalClusDet, "EcalClusDet[20]/I");
   myTree->Branch("GammaRecoEt",  GammaRecoEt, "GammaRecoEt[20]/F");
   myTree->Branch("GammaRecoEta",  GammaRecoEta, "GammaRecoEta[20]/F");
   myTree->Branch("GammaRecoPhi",  GammaRecoPhi, "GammaRecoPhi[20]/F");
   myTree->Branch("GammaIsoEcal",  GammaIsoEcal, "GammaIsoEcal[9][20]/F");

// Tracks block
   myTree->Branch("NumRecoTrack", &NumRecoTrack, "NumRecoTrack/I");
   myTree->Branch("TrackRecoEt",  TrackRecoEt, "TrackRecoEt[200]/F");
   myTree->Branch("TrackRecoEta",  TrackRecoEta, "TrackRecoEta[200]/F");
   myTree->Branch("TrackRecoPhi",  TrackRecoPhi, "TrackRecoPhi[200]/F");

// Particle block
   myTree->Branch("NumPart", &NumPart, "NumPart/I"); 
   myTree->Branch("Status",  Status, "Status[4000]/I");     
   myTree->Branch("Code",  Code, "Code[4000]/I");
   myTree->Branch("Mother1",  Mother1, "Mother1[4000]/I");
   myTree->Branch("partpx",  partpx, "partpx[4000]/F");
   myTree->Branch("partpy",  partpy, "partpy[4000]/F");
   myTree->Branch("partpz",  partpz, "partpz[4000]/F");
   myTree->Branch("parte",  parte, "parte[4000]/F");
   myTree->Branch("partm",  partm, "partm[4000]/F");
   myTree->Branch("partvx",  partvx, "partvx[4000]/F");
   myTree->Branch("partvy",  partvy, "partvy[4000]/F");
   myTree->Branch("partvz",  partvz, "partvz[4000]/F");
   myTree->Branch("partvt",  partvt, "partvt[4000]/F");
// end of tree declaration

   edm::ESHandle<CaloGeometry> pG;
   iSetup.get<IdealGeometryRecord>().get(pG);
   geo = pG.product();

  myout_part = new ofstream((myName+"_part.dat").c_str()); 
  if(!myout_part) cout << " Output file not open!!! "<<endl;
  myout_hcal = new ofstream((myName+"_hcal.dat").c_str()); 
  if(!myout_hcal) cout << " Output file not open!!! "<<endl;
  myout_ecal = new ofstream((myName+"_ecal.dat").c_str()); 
  if(!myout_ecal) cout << " Output file not open!!! "<<endl;
  
  myout_jet = new ofstream((myName+"_jet.dat").c_str()); 
  if(!myout_jet) cout << " Output file not open!!! "<<endl;
  myout_photon = new ofstream((myName+"_photon.dat").c_str()); 
  if(!myout_photon) cout << " Output file not open!!! "<<endl;
   
   
}

void GammaJetAnalysis::endJob()
{

   cout << "===== Start writing user histograms =====" << endl;
   hOutputFile->SetCompressionLevel(2);
   hOutputFile->cd();
   myTree->Write();
   hOutputFile->Close() ;
   cout << "===== End writing user histograms =======" << endl; 
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
GammaJetAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  std::vector<Provenance const*> theProvenance;
  iEvent.getAllProvenance(theProvenance);
  for( std::vector<Provenance const*>::const_iterator ip = theProvenance.begin();
                                                      ip != theProvenance.end(); ip++)
  {
     cout<<" Print all module/label names "<<(**ip).moduleName()<<" "<<(**ip).moduleLabel()<<
     " "<<(**ip).productInstanceName()<<endl;
  }
// Load generator information
// write HEPEVT block into file
   run = iEvent.id().run();
   event = iEvent.id().event();
  (*myout_part)<<"Event "<<iEvent.id().run()<<" "<<iEvent.id().event()<<endl;
  (*myout_jet)<<"Event "<<iEvent.id().run()<<" "<<iEvent.id().event()<<endl;
  (*myout_hcal)<<"Event "<<iEvent.id().run()<<" "<<iEvent.id().event()<<endl;
  (*myout_ecal)<<"Event "<<iEvent.id().run()<<" "<<iEvent.id().event()<<endl;
  (*myout_photon)<<"Event "<<iEvent.id().run()<<" "<<iEvent.id().event()<<endl;
  
  int nevhep = 0;
  vector<HepMC::GenParticle*> theGenPart;
  vector<HepMC::GenParticle*> theGenPartAll;
  map<HepMC::GenParticle*,int> theGenPartPair;
  
  try {
    vector< Handle< HepMCProduct > > EvtHandles ;
    iEvent.getManyByType( EvtHandles ) ;
//    cout<<" Size of HepMCProduct "<<EvtHandles.size()<<endl;
    
   for ( unsigned int i=0; i<EvtHandles.size(); i++ )
   {

      if ( EvtHandles[i].isValid() )
      {

         const HepMC::GenEvent* Evt = EvtHandles[i]->GetEvent() ;

         // take only 1st vertex for now - it's been tested only of PGuns...
         //
         HepMC::GenEvent::vertex_const_iterator Vtx = Evt->vertices_begin() ;

         std::cout<<" ParticleId 7/8= "<<Evt->particle(7)->pdg_id()<<" "<<
                                         Evt->particle(8)->pdg_id()<<std::endl;
         std::cout<<" ParticleId 7/8= "<<Evt->particle(7)->momentum().perp()<<" "<<
                                         Evt->particle(8)->momentum().perp()<<std::endl;

	 
         for (HepMC::GenEvent::particle_const_iterator
            Part = Evt->particles_begin() ; Part!=Evt->particles_end(); Part++ )
         {
            HepLorentzVector Mom = (*Part)->momentum() ;
            double Eta = -log(tan(Mom.theta()/2.));

         if ( EvtHandles[i].provenance()->moduleLabel() == "VtxSmeared" )
            {
//	     cout<<" particle "<<(*Part)->pdg_id()<<" "<<(*Part)->status()<<" "<<Eta<<" "<<Mom.px()<<endl;
	     
	     nevhep++;
	     theGenPartPair[(*Part)] = nevhep;   
	     theGenPartAll.push_back(*Part);
	     
// Analysis of the particles in event
                if( !(*Part)->production_vertex() ) std::cout<<"Myparticle "<<(*Part)->pdg_id()<<std::endl;
                if( (*Part)->status()!=1 ) continue;
                if ( (*Part)->pdg_id()!=22 ) continue;
		if(fabs(Eta) > 2.8 ) continue;
		if(Mom.perp() <5. ) continue;
		
		theGenPart.push_back(*Part);
		
		if((*Part)->production_vertex()) {
		   cout<<" Mother "<<(*Part)->mother()->pdg_id()<<endl;
		 } else { cout<<" There is not mother "<<endl;
		 }  
		if((*Part)->mother()->production_vertex()) {
		   cout<<" GrandMother "<<(*Part)->mother()->mother()->pdg_id()<<endl;
		} else
		{
		    cout<<" There is not grandmother "<<endl;  
		}   
            }
            else
            {
	     cout<<" Label "<<EvtHandles[i].provenance()->moduleLabel()<<endl;
            }
         }

      }

   }
  
  } catch (std::exception& e) { // can't find it!
    if (!allowMissingInputs_) throw e;
  }
  
//  cout<<" We filled theGenPart "<<endl;
  
  (*myout_part)<<nevhep<<endl;
  int ihep = 0;
  int motherline1 = 0;
  int motherline2 = 0;
  int dauthline1 = 0;
  int dauthline2 = 0;
  
  for(vector<HepMC::GenParticle*>::iterator Part = theGenPartAll.begin(); Part != theGenPartAll.end(); Part++)
  {
    ihep++;
// If there is a mother
    motherline1 = 0;
    motherline2 = 0;
    dauthline1 = 0;
    dauthline2 = 0;
    double vx=0.; 
    double vy=0.;
    double vz=0.; 
    double vt=0.;
    
    if((*Part)->production_vertex())
    {
        if(theGenPartPair.find((*Part)->mother()) != theGenPartPair.end() ) 
	               motherline1 = (*theGenPartPair.find((*Part)->mother())).second;
	if(theGenPartPair.find((*Part)->secondMother()) != theGenPartPair.end()) 
	{
	   motherline2 = (*theGenPartPair.find((*Part)->secondMother())).second;
	} else  motherline2 = 0;
	vx = (*Part)->production_vertex()->position().x();
	vy = (*Part)->production_vertex()->position().y();
	vy = (*Part)->production_vertex()->position().z();
	vt = (*Part)->production_vertex()->position().t();
    } 
    if((*Part)->end_vertex())
    {
         if(theGenPartPair.find((*Part)->beginDaughters()) != theGenPartPair.end()) dauthline1 = (*theGenPartPair.find((*Part)->beginDaughters())).second;
	 if(theGenPartPair.find((*Part)->endDaughters()-1) != theGenPartPair.end()) dauthline2 = (*theGenPartPair.find((*Part)->endDaughters()-1)).second;
    }
    
    (*myout_part)<<ihep<<" "<<(*Part)->status()<<" "<<(*Part)->pdg_id()
                       <<" "<<motherline1<<" "<<motherline2
                       <<" "<<dauthline1<<" "<<dauthline2
		       <<" "<<(*Part)->momentum().px()<<" "<<(*Part)->momentum().py()
		       <<" "<<(*Part)->momentum().pz()
                       <<" "<<(*Part)->momentum().e()<<" "<<(*Part)->momentum().m()
                       <<" "<<vx<<" "<<vy<<" "<<vz<<" "<<vt
     <<endl;
     if(ihep-1<4000)
     {
       Status[ihep-1] = (*Part)->status();
       Code[ihep-1] = (*Part)->pdg_id();
       Mother1[ihep-1] = motherline1;
       partpx[ihep-1] = (*Part)->momentum().px();
       partpy[ihep-1] = (*Part)->momentum().py();
       partpz[ihep-1] = (*Part)->momentum().pz();
       parte[ihep-1] = (*Part)->momentum().e();
       partm[ihep-1] = (*Part)->momentum().m();
       partvx[ihep-1] = vx;
       partvy[ihep-1] = vy;
       partvz[ihep-1] = vz;
       partvt[ihep-1] = vt;
       NumPart = ihep;
       
     } 
  }
     
  int ntype=0;
  for(vector<HepMC::GenParticle*>::const_iterator Part = theGenPart.begin(); Part != theGenPart.end(); Part++)
  {
//    cout<<"Particle "<<(*Part)->pdg_id()<<endl;
//    if( (*Part)->production_vertex() ) cout<<"Mother "<<(*Part)->mother()->pdg_id()<<endl;
//    if( (*Part)->mother()->production_vertex() ) cout<<"Mother "<<(*Part)->mother()->mother()->pdg_id()<<endl;

        if( !(*Part)->production_vertex() ) 
	{
	   ntype = 1;
	   continue;
	}   

	if( (*Part)->mother()->pdg_id() == 223 || (*Part)->mother()->mother()->pdg_id() == 223 ) 
	{
		 ntype = 5;
        }
        if( (*Part)->mother()->pdg_id() == 221 || (*Part)->mother()->mother()->pdg_id() == 221 ) 
	{
		 ntype = 4;
        }
        if( (*Part)->mother()->pdg_id() == 111 )
	{
		
	   if ((*Part)->mother()->mother()->pdg_id() != 221 && (*Part)->mother()->mother()->pdg_id() != 223 ) 
	   {
		 ntype = 3;
           }
        }
	if( (*Part)->mother()->pdg_id() != 111 && (*Part)->mother()->pdg_id() != 221 	
	                         && (*Part)->mother()->pdg_id() != 223) 
        {
				  ntype = 2;
	}
  }
  (*myout_part)<<endl;
  (*myout_part)<<theGenPart.size()<<endl; 
// Load Jets Calo and Gen

    std::vector<edm::InputTag>::const_iterator ic;
    int jettype = 0;
    int jetexist = -100;
    int reco = 1;
    double etlost = -100.1;
    
    NumRecoJets = 0;
    
    for (ic=mInputCalo.begin(); ic!=mInputCalo.end(); ic++) {
     try {
       
       edm::Handle<reco::CaloJetCollection> jets;
       iEvent.getByLabel(*ic, jets);
       reco::CaloJetCollection::const_iterator jet = jets->begin ();
//       cout<<" Size of jets "<<jets->size()<<endl;
        jettype++;
       
        
       
       if(jets->size() > 0 )
       {
         int ij = 0;
         for (; jet != jets->end (); jet++) 
         {
//            cout<<*ic<<" "<<" et "<<(*jet).et()<<" "<<(*jet).eta()<<" "<<(*jet).phi()<<endl;
	    ij++;
	    if(ij<4) (*myout_jet)<<jettype<<" "<<reco<<" "<<ij<<" "<<(*jet).et()<<" "<<(*jet).eta()<<" "<<(*jet).phi()
	    <<" "<<iEvent.id().run()<<" "<<iEvent.id().event()<<endl;
	    jetexist = ij;
	    if( NumRecoJets < 10 )
	    {
	     JetRecoEt[NumRecoJets] = (*jet).et();
	     JetRecoEta[NumRecoJets] = (*jet).eta();
	     JetRecoPhi[NumRecoJets] = (*jet).phi();
	     JetRecoType[NumRecoJets] = jettype;
	     NumRecoJets++;
	    }
         }
       }
     } catch (std::exception& e) { // can't find it!
       if (!allowMissingInputs_) throw e;
     }
   }
     if( jetexist < 0 ) (*myout_jet)<<jetexist<<" "<<reco<<" "<<etlost
                         <<" "<<etlost<<" "<<etlost<<" "<<iEvent.id().run()
			 <<" "<<iEvent.id().event()<<endl;
    jettype = 0;
    jetexist = -100;
    reco = 0;
    NumGenJets = 0;
    
    for (ic=mInputGen.begin(); ic!=mInputGen.end(); ic++) {
     try {
  
       edm::Handle<reco::GenJetCollection> jets;
       iEvent.getByLabel(*ic, jets);
       reco::GenJetCollection::const_iterator jet = jets->begin ();
//       cout<<" Size of jets "<<jets->size()<<endl;
       jettype++;
       
       
       if(jets->size() > 0 )
       {
          int ij = 0;
         for (; jet != jets->end (); jet++) 
         {
	     ij++;
//            cout<<*ic<<" "<<" et "<<(*jet).et()<<" "<<(*jet).eta()<<" "<<(*jet).phi()<<endl;
	    if(ij<4) (*myout_jet)<<jettype<<" "<<reco<<" "<<ij<<" "<<(*jet).et()<<" "<<(*jet).eta()<<" "<<(*jet).phi()
	    <<" "<<iEvent.id().event()<<endl;
	    jetexist = jettype;	
	    if( NumGenJets < 10 )
	    {    
	      JetGenEt[NumGenJets] = (*jet).et();
	      JetGenEta[NumGenJets] = (*jet).eta();
	      JetGenPhi[NumGenJets] = (*jet).phi();
	      JetGenType[NumGenJets] = jettype;
	      NumGenJets++;
	    }

         }
       }
     } catch (std::exception& e) { // can't find it!
       if (!allowMissingInputs_) throw e;
     }
   }
   
     if( jetexist < 0 ) (*myout_jet)<<jetexist<<" "<<reco<<" "<<etlost
                         <<" "<<etlost<<" "<<etlost<<" "<<iEvent.id().run()
			 <<" "<<iEvent.id().event()<<endl;

// Load EcalRecHits
  
    std::vector<edm::InputTag>::const_iterator i;
    vector<CaloRecHit> theRecHits;
      
    for (i=ecalLabels_.begin(); i!=ecalLabels_.end(); i++) {
    try {
      
      edm::Handle<EcalRecHitCollection> ec;
      iEvent.getByLabel(*i,ec);
      
      (*myout_ecal)<<(*ec).size()<<endl;

       for(EcalRecHitCollection::const_iterator recHit = (*ec).begin();
                                                recHit != (*ec).end(); ++recHit)
       {
// EcalBarrel = 1, EcalEndcap = 2

	 GlobalPoint pos = geo->getPosition(recHit->detid());
         theRecHits.push_back(*recHit);

	 (*myout_ecal)<<recHit->detid().subdetId()<<" "<<(*recHit).energy()<<" "<<pos.phi()<<" "<<pos.eta()<<endl;    
       } 
      
    } catch (std::exception& e) { // can't find it!
    if (!allowMissingInputs_) throw e;
    }
    }

//  cout<<" Start to get hbhe "<<endl;
// Hcal Barrel and endcap for isolation  
    try {
      edm::Handle<HBHERecHitCollection> hbhe;
      iEvent.getByLabel(hbheLabel_,hbhe);
      
      (*myout_hcal)<<(*hbhe).size()<<endl;
  for(HBHERecHitCollection::const_iterator hbheItr = (*hbhe).begin();
  
      hbheItr != (*hbhe).end(); ++hbheItr)
      {
        DetId id = (hbheItr)->detid();
	GlobalPoint pos = geo->getPosition(hbheItr->detid());
	(*myout_hcal)<<id.subdetId()<<" "<<(*hbheItr).energy()<<" "<<pos.phi()<<" "<<pos.eta()<<endl;    
        theRecHits.push_back(*hbheItr);
	
      }
    } catch (std::exception& iEvent) { // can't find it!
      cout<<" Exception in hbhe "<<endl;
      if (!allowMissingInputs_) throw iEvent;
    }
//  }

 
 for(int i = 0; i<9; i++)
 {
    for(int j= 0; j<10; j++) GammaIsoEcal[i][j] = 0.;
 }

// Load Ecal clusters
 jetexist = -100; 
 int barrel = 1;
 NumRecoGamma = 0;
 try {
 int ij = 0;
  // Get island super clusters after energy correction
  Handle<reco::SuperClusterCollection> pCorrectedIslandBarrelSuperClusters;
  iEvent.getByLabel(correctedIslandBarrelSuperClusterProducer_, correctedIslandBarrelSuperClusterCollection_, pCorrectedIslandBarrelSuperClusters);
  const reco::SuperClusterCollection* correctedIslandBarrelSuperClusters = pCorrectedIslandBarrelSuperClusters.product();
  // loop over the super clusters and fill the histogram
  for(reco::SuperClusterCollection::const_iterator aClus = correctedIslandBarrelSuperClusters->begin();
                                                           aClus != correctedIslandBarrelSuperClusters->end(); aClus++) {
    double vet = aClus->energy()/cosh(aClus->eta());
//    std::cout<<" SuperCluster Island clusters barrel after correction "<<vet<<std::endl;
    if(vet>20.) {
    
      ij++;
//      std::cout<<" SuperCluster energy barrel  "<<ij<<" "<<vet<<std::endl;
      NumRecoGamma++;
      float gammaiso_ecal[9] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
      
      for(vector<CaloRecHit>::iterator it = theRecHits.begin(); it != theRecHits.end(); it++)
      {
           GlobalPoint pos = geo->getPosition(it->detid());
           double eta = pos.eta();
	   double phi = pos.phi();
	   double deta = fabs(eta-aClus->eta());
	   double dphi = fabs(phi-aClus->phi());
	   if(dphi>4.*atan(1.)) dphi = 8.*atan(1.)-dphi;
	   double dr = sqrt(deta*deta+dphi*dphi);
	   
	   double rmin = 0.07;
	   if( fabs(aClus->eta()) > 1.47 ) rmin = 0.07*(fabs(aClus->eta())-.47)*1.2;
	   if( fabs(aClus->eta()) > 2.2 ) rmin = 0.07*(fabs(aClus->eta())-.47)*1.4;
	   
	   int itype_ecal = 0;
	   double ecutn = 0.;
	   for (int i = 0; i<3; i++)
	   {
	   
//             std::cout<<" Cluster-rechits "<<aClus->eta()<<" "<<eta<<" "<<aClus->phi()<<" "<<phi<<" dr "<<dr
//	     <<" energy "<<(*it).energy()<<" et "<<(*it).energy()/cosh(eta)<<std::endl;	

	     for (int j = 0; j<3; j++)
	     {
	     
	        if(it->detid().det() == DetId::Ecal ) 
		{
		  if(it->detid().subdetId() == 1) ecutn = ecut[0][j];
		  if(it->detid().subdetId() == 2) ecutn = ecut[1][j];
		  if( dr>rmin && dr<risol[i])
		  {
		   if((*it).energy() > ecutn) gammaiso_ecal[itype_ecal] = gammaiso_ecal[itype_ecal]+(*it).energy()/cosh(eta);
		  } 
		}
		
		if(it->detid().det() == DetId::Hcal ) 
		{
		   ecutn = ecut[2][j];
//		       cout<<" Hcal Tower "<<eta<<" "<<phi<<" "<<dr<<" "<<(*it).energy()<<" "<<ecutn<<endl;
		   if( dr>rmin && dr<risol[i])
		   {
		     if((*it).energy() > ecutn) 
		     {
		        gammaiso_ecal[itype_ecal] = gammaiso_ecal[itype_ecal]+(*it).energy()/cosh(eta);
//			cout<<" Hcal energy "<<eta<<" "<<phi<<" "<<(*it).energy()/cosh(eta)<<" "<<gammaiso_ecal[itype_ecal]<<endl;
		     }
		   }
		} 
		jetexist = ij;
		if( ij < 10 ) GammaIsoEcal[itype_ecal][NumRecoGamma-1] = gammaiso_ecal[itype_ecal]; 
	        itype_ecal++;
		
	     } // Ecal
	       
	   } // cycle on iso radii
	   
      }
// Fill Tree      
      EcalClusDet[NumRecoGamma-1] = 1;
      GammaRecoEt[NumRecoGamma-1] = vet;
      GammaRecoEta[NumRecoGamma-1] = aClus->eta();
      GammaRecoPhi[NumRecoGamma-1] = aClus->phi();
      
      
    (*myout_photon)<<ij<<" "<<barrel<<" "<<vet<<" "<<aClus->eta()<<" "<<aClus->phi()<<" "<<iEvent.id().event()<<endl;
    (*myout_photon)<<ij<<" "<<gammaiso_ecal[0]<<" "<<gammaiso_ecal[1] <<" "<<gammaiso_ecal[2]<<" "<<gammaiso_ecal[3]
      <<" "<<gammaiso_ecal[4]<<" "<<gammaiso_ecal[5]<<" "<<gammaiso_ecal[6]<<" "<<gammaiso_ecal[7]<<" "<<gammaiso_ecal[8]<<endl;
//    (*myout_photon)<<ij<<" "<<gammaiso_hcal[0]<<" "<<gammaiso_hcal[1] <<" "<<gammaiso_hcal[2]<<" "<<gammaiso_hcal[3]
//      <<" "<<gammaiso_hcal[4]<<" "<<gammaiso_hcal[5]<<" "<<gammaiso_hcal[6]<<" "<<gammaiso_hcal[7]<<" "<<gammaiso_hcal[8]<<endl;
      
       jetexist = ij;
    }  
  }
  } catch (std::exception& e) { // can't find it!
    if (!allowMissingInputs_) throw e;
  }

  barrel = 2;
  
  try {
  int ij = 0;
  // Get island super clusters after energy correction
  Handle<reco::SuperClusterCollection> pCorrectedIslandEndcapSuperClusters;
  iEvent.getByLabel(correctedIslandEndcapSuperClusterProducer_, correctedIslandEndcapSuperClusterCollection_, pCorrectedIslandEndcapSuperClusters);
  const reco::SuperClusterCollection* correctedIslandEndcapSuperClusters = pCorrectedIslandEndcapSuperClusters.product();
  
  // loop over the super clusters and fill the histogram
  for(reco::SuperClusterCollection::const_iterator aClus = correctedIslandEndcapSuperClusters->begin();
                                                           aClus != correctedIslandEndcapSuperClusters->end(); aClus++) {
    double vet = aClus->energy()/cosh(aClus->eta());
    
//    std::cout<<" Cluster energy in endcap "<<vet<<std::endl;
    
    if(vet>20.) {
    ij++;
//    std::cout<<" Cluster energy in endcap number "<<ij<<" "<<vet<<std::endl;
    NumRecoGamma++;
    float gammaiso_ecal[9] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
    
      for(vector<CaloRecHit>::iterator it = theRecHits.begin(); it != theRecHits.end(); it++)
      {
           GlobalPoint pos = geo->getPosition(it->detid());
           double eta = pos.eta();
	   double phi = pos.phi();
	   double deta = fabs(eta-aClus->eta());
	   double dphi = fabs(phi-aClus->phi());
	   if(dphi>4.*atan(1.)) dphi = 8.*atan(1.)-dphi;
	   double dr = sqrt(deta*deta+dphi*dphi);
	   
	   double rmin = 0.07;
	   if( fabs(aClus->eta()) > 1.47 ) rmin = 0.07*(fabs(aClus->eta())-.47)*1.2;
	   if( fabs(aClus->eta()) > 2.2 ) rmin = 0.07*(fabs(aClus->eta())-.47)*1.4;
	   
	   int itype_ecal = 0;

	   double ecutn = 0.;
	   
//	     std::cout<<" Cluster-rechits "<<aClus->eta()<<" "<<eta<<" "<<aClus->phi()<<" "<<phi<<" dr "<<dr
//	     <<" energy "<<(*it).energy()<<" et "<<(*it).energy()/cosh(eta)<<std::endl;
	     
	   for (int i = 0; i<3; i++)
	   {
	     for (int j = 0; j<3; j++)
	     {		
	     
	        if(it->detid().det() == DetId::Ecal ) 
		{
		  if(it->detid().subdetId() == 1) ecutn = ecut[0][j];
		  if(it->detid().subdetId() == 2) ecutn = ecut[1][j];
		  if( dr>rmin && dr<risol[i])
		  {
		   if((*it).energy() > ecutn) gammaiso_ecal[itype_ecal] = gammaiso_ecal[itype_ecal]+(*it).energy()/cosh(eta);
		  } 
		}
		
		if(it->detid().det() == DetId::Hcal ) 
		{		   
		   if( dr>rmin && dr<risol[i])
		   {
		     ecutn = ecut[2][j];
		     if((*it).energy() > ecutn) gammaiso_ecal[itype_ecal] = gammaiso_ecal[itype_ecal]+(*it).energy()/cosh(eta);
		   }
		} 
		
		if( ij < 10 ) GammaIsoEcal[itype_ecal][NumRecoGamma-1] = gammaiso_ecal[itype_ecal]; 
	        itype_ecal++;

	     }
	       
	   }
	   
      }
//     std::cout<<" Fill endcap ecal cluster tree "<<std::endl; 
// Fill Tree    
    EcalClusDet[NumRecoGamma-1] = 2;
    GammaRecoEt[NumRecoGamma-1] = vet;
    GammaRecoEta[NumRecoGamma-1] = aClus->eta();
    GammaRecoPhi[NumRecoGamma-1] = aClus->phi();
////    
    (*myout_photon)<<ij<<" "<<barrel<<" "<<vet<<" "<<aClus->eta()<<" "<<aClus->phi()<<" "<<iEvent.id().event()<<endl;
    (*myout_photon)<<ij<<" "<<gammaiso_ecal[0]<<" "<<gammaiso_ecal[1] <<" "<<gammaiso_ecal[2]<<" "<<gammaiso_ecal[3]
      <<" "<<gammaiso_ecal[4]<<" "<<gammaiso_ecal[5]<<" "<<gammaiso_ecal[6]<<" "<<gammaiso_ecal[7]<<" "<<gammaiso_ecal[8]<<endl;
//    (*myout_photon)<<ij<<" "<<gammaiso_hcal[0]<<" "<<gammaiso_hcal[1] <<" "<<gammaiso_hcal[2]<<" "<<gammaiso_hcal[3]
//      <<" "<<gammaiso_hcal[4]<<" "<<gammaiso_hcal[5]<<" "<<gammaiso_hcal[6]<<" "<<gammaiso_hcal[7]<<" "<<gammaiso_hcal[8]<<endl;

      jetexist = ij;
    }							   
//   std::cout<<" Endcap SuperClusters after corrections "<<aClus->energy()<<std::endl;
    
  }
  // ----- hybrid 
  } catch (std::exception& e) { // can't find it!
    if (!allowMissingInputs_) throw e;
  }
    double ecluslost = -100.1;
    if(jetexist<0) (*myout_photon)<<jetexist<<" "<<barrel<<" "<<ecluslost<<" "<<ecluslost<<" "<<ecluslost<<" "<<iEvent.id().run()<<" "<<iEvent.id().event()<<endl;
  
//  }
//  if (!hoLabel_.empty()) {
// Load HORecHits
    try {
      edm::Handle<HORecHitCollection> ho;
      iEvent.getByLabel(hoLabel_,ho);
      
      (*myout_hcal)<<(*ho).size()<<endl;
      
  for(HORecHitCollection::const_iterator hoItr = (*ho).begin();
      hoItr != (*ho).end(); ++hoItr)
      {
        DetId id = (hoItr)->detid();
	GlobalPoint pos = geo->getPosition(hoItr->detid());
	(*myout_hcal)<<id.subdetId()<<" "<<(*hoItr).energy()<<" "<<pos.phi()<<" "<<pos.eta()<<endl;    

      }
    } catch (std::exception& iEvent) { // can't find it!
      if (!allowMissingInputs_) throw iEvent;
    }
    
// Load HFRecHits
    try {
      edm::Handle<HFRecHitCollection> hf;
      iEvent.getByLabel(hfLabel_,hf);
      
      (*myout_hcal)<<(*hf).size()<<endl;
  for(HFRecHitCollection::const_iterator hfItr = (*hf).begin();
      hfItr != (*hf).end(); ++hfItr)
      {  
         DetId id = (hfItr)->detid();
	GlobalPoint pos = geo->getPosition(hfItr->detid());
	(*myout_hcal)<<id.subdetId()<<" "<<(*hfItr).energy()<<" "<<pos.phi()<<" "<<pos.eta()<<endl;    

      }
    } catch (std::exception& iEvent) { // can't find it!
      if (!allowMissingInputs_) throw iEvent;
    }
//  }

// Load CaloTowers
//  Handle<CaloTowerCollection> calotowers;
//  iEvent.getByType(calotowers);
//  int idx = 0;
//  for (; idx < calotowers->size (); idx++) {
//    const CaloTower* cal = &((*calotowers) [idx]);

//           cout<<" Et = "<<cal->et()
//               <<" Number of constituents "<< cal->constituentsSize()
//               <<" eta = " << cal->eta()
//               <<" phi = " << cal->phi()
//               <<" emE = " << cal->emEnergy()
//               <<" emEt = " << cal->emEt()
//               <<" detid = " << cal->id().rawId()
//               <<" hadEt = " << cal->hadEt() << endl;
//   }
// Load tracks
//   cout<<" Fill NumPart "<<NumPart<<endl;

// Load Tracks
   
   
   myTree->Fill();
   
} // analyze method
} // namespace cms
