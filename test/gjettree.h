//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Dec  7 11:38:21 2007 by ROOT version 5.14/00f
// from TTree pippo/Analysis tree
// found on file: output.root
//////////////////////////////////////////////////////////

#ifndef gjettree_h
#define gjettree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TTree.h>

#include <TMVAPhotons_MLP.class.C>
#include <TMVAPhotons_1_MLP.class.C>
#include <TMVAPhotons_2_MLP.class.C>
#include <TMVAPhotons_3_MLP.class.C>
#include <TMVAPhotons_4_MLP.class.C>

#include "Config.hpp"

#include <map>
#include <string>

struct photoncandidate {
  double e;
  double pt;
  double eta;
  double phi;
};

struct jetcandidate {
  double e;
  double pt;
  double eta;
  double phi;
};

struct photonidcuts {
  float hcaliso;
  float ecaliso;
  float trackiso;
  int tracknb;
  float sminmin;
  float smajmaj;
};

class gjettree_histos {
  public :
    
  // CODING: L=leading/S=second; R=reco/G=gen/P=parton; P=photon/J=jet
  TH1D *lrp_pt;
  TH1D *lrp_eta;
  TH1D *lrp_phi;

  TH1D *lrj_pt;
  TH1D *lrj_eta;
  TH1D *lrj_phi;

  TH1D *lgj_pt;
  TH1D *lgj_eta;
  TH1D *lgj_phi;

  TH1D *srp_pt;
  TH1D *srp_eta;
  TH1D *srp_phi;

  TH1D *srj_pt;
  TH1D *srj_eta;
  TH1D *srj_phi;

  // Studies on additional radiation
  TH1D *srj_ptsum;
  TH1D *srj_ptvecsum;
  TH1D *srj_ptsumop;
  TH1D *srj_ptvecsumop;

  // ptop=PtOverPtPhoton
  TH1D *lrj_ptop;
  TH1D *lgj_ptop;
  TH1D *lgq_ptop;
  TH1D *srp_ptop;
  TH1D *srj_ptop;
  // ptogg=PtOverPt(Gen)Gamma
  TH1D *lrj_ptogg;
  TH1D *lgj_ptogg;
  // dphivp=DeltaPhiVsPhoton
  TH1D *lrj_dphivp;
  TH1D *lgj_dphivp;
  TH1D *srp_dphivp;
  TH1D *srj_dphivp;
  // dphivp=DeltaRVsPhoton
  TH1D *srj_drvp;
  // dphivgg=DeltaPhiVs(Gen)Gamma
  TH1D *lrj_dphivgg;
  TH1D *lgj_dphivgg;
  // detap=DeltaEtaVsPhoton
  TH1D *lrj_detavp;
  TH1D *lgj_detavp;
  TH1D *srp_detavp;
  TH1D *srj_detavp;
  // detap=DeltaEtaVsGenGamma
  TH1D *lrj_detavgg;
  TH1D *lgj_detavgg;

  // dr=DeltaR vg=VsGenObject (GenJet, parton photon) vgg=GenGamma
  TH1D *lrp_devg;
  TH1D *lrp_dptvg;
  TH1D *lrp_drvg;
  TH1D *lrp_detavg;
  TH1D *lrp_dphivg;

  TH1D *lrp_devgg;
  TH1D *lrp_dptvgg;
  TH1D *lrp_drvgg;
  TH1D *lrp_detavgg;
  TH1D *lrp_dphivgg;

  TH1D *lgg_devg;
  TH1D *lgg_dptvg;
  TH1D *lgg_drvg;
  TH1D *lgg_detavg;
  TH1D *lgg_dphivg;


  TH1D *lrj_devg;
  TH1D *lrj_dptvg;
  TH1D *lrj_drvg;
  TH1D *lrj_detavg;
  TH1D *lrj_dphivg;

  TH1D *lrj_devq;
  TH1D *lrj_dptvq;
  TH1D *lrj_drvq;
  TH1D *lrj_detavq;
  TH1D *lrj_dphivq;

  TH1D *lgj_devq;
  TH1D *lgj_dptvq;
  TH1D *lgj_drvq;
  TH1D *lgj_detavq;
  TH1D *lgj_dphivq;

  // PhotonID histograms
  /*
  TH1D *pid_ntrkiso035;
  TH1D *pid_ptiso035;
  TH1D *pid_hcalovecal04;
  TH1D *pid_sMajMaj;
  TH1D *pid_sMinMin;
  TH1D *pid_ecaliso04;
  */
  
  TH1D *lrp_ntrkiso;
  TH1D *lrp_ptiso;
  TH1D *lrp_emf;
  TH1D *lrp_sMajMaj;
  TH1D *lrp_sMinMin;
  TH1D *lrp_ecaliso_old;
  TH1D *lrp_ecaliso_new;

  TH1D *srp_ntrkiso;
  TH1D *srp_ptiso;
  TH1D *srp_emf;
  TH1D *srp_sMajMaj;
  TH1D *srp_sMinMin;
  TH1D *srp_ecaliso_old;
  TH1D *srp_ecaliso_new;

  gjettree_histos();
  ~gjettree_histos();
  void Write();
};

class gjettree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leave types
   Int_t           event;
   Int_t           nMC;
   Int_t           pdgIdMC[1500];   //[nMC] was 1446, 100 enough
   Int_t           statusMC[1500];   //[nMC] 
   Int_t           motherIDMC[1500];   //[nMC]
   Float_t         pxMC [1500];   //[nMC] 
   Float_t         pyMC [1500];   //[nMC]
   Float_t         pzMC [1500];   //[nMC]
   Float_t         eMC  [1500];   //[nMC]
   Float_t         etaMC[1500];   //[nMC]
   Float_t         phiMC[1500];   //[nMC]
   Int_t           nPhot;
   Float_t         pxPhot [25];   //[nPhot] was 25, 17 enough
   Float_t         pyPhot [25];   //[nPhot]
   Float_t         pzPhot [25];   //[nPhot]
   Float_t         ePhot  [25];   //[nPhot]
   Float_t         etaPhot[25];   //[nPhot]
   Float_t         phiPhot[25];   //[nPhot]
   Float_t         ptiso02Phot[25];   //[nPhot]
   Int_t           ntrkiso02Phot[25];   //[nPhot]
   Float_t         ptiso025Phot[25];   //[nPhot]
   Int_t           ntrkiso025Phot[25];   //[nPhot]
   Float_t         ptiso03Phot[25];   //[nPhot]
   Int_t           ntrkiso03Phot[25];   //[nPhot]
   Float_t         ptiso035Phot[25];   //[nPhot]
   Int_t           ntrkiso035Phot[25];   //[nPhot]
   Float_t         ptiso04Phot[25];   //[nPhot]
   Float_t         ptisoatecal02Phot[25];   //[nPhot]
   Int_t           ntrkisoatecal02Phot[25];   //[nPhot]
   Float_t         ptisoatecal025Phot[25];   //[nPhot]
   Int_t           ntrkisoatecal025Phot[25];   //[nPhot]
   Float_t         ptisoatecal03Phot[25];   //[nPhot]
   Int_t           ntrkisoatecal03Phot[25];   //[nPhot]
   Float_t         ptisoatecal035Phot[25];   //[nPhot]
   Int_t           ntrkisoatecal035Phot[25];   //[nPhot]
   Float_t         ptisoatecal04Phot[25];   //[nPhot]
   Int_t           ntrkiso04Phot[25];   //[nPhot]
   Float_t         hcalovecal01Phot[25];   //[nPhot]
   Float_t         hcalovecal02Phot[25];   //[nPhot]
   Float_t         hcalovecal025Phot[25];   //[nPhot]
   Float_t         hcalovecal03Phot[25];   //[nPhot]
   Float_t         hcalovecal04Phot[25];   //[nPhot]
   Float_t         hcalovecal05Phot[25];   //[nPhot]
   Float_t         ecaliso02Phot[25];   //[nPhot]
   Float_t         ecaliso03Phot[25];   //[nPhot]
   Float_t         ecaliso035Phot[25];   //[nPhot]
   Float_t         ecaliso04Phot[25];   //[nPhot]
   Float_t         ecaliso05Phot[25];   //[nPhot]
   Float_t         LATPhot[25];   //[nPhot] was 15, 17 enough
   Float_t         sMajMajPhot[25];   //[nPhot]
   Float_t         sMinMinPhot[25];   //[nPhot]
   Float_t         FisherPhot[25];   //[nPhot]
   Int_t           nJet;
   Float_t         pxJet [50];   //[nJet] was 44, 43 enough
   Float_t         pyJet [50];   //[nJet]
   Float_t         pzJet [50];   //[nJet]
   Float_t         eJet  [50];   //[nJet]
   Float_t         etaJet[50];   //[nJet]
   Float_t         phiJet[50];   //[nJet]
   Int_t           nJetGen;
   Float_t         pxJetGen [100];   //[nJetGen] was 69, 72 enough
   Float_t         pyJetGen [100];   //[nJetGen]
   Float_t         pzJetGen [100];   //[nJetGen]
   Float_t         eJetGen  [100];   //[nJetGen]
   Float_t         etaJetGen[100];   //[nJetGen]
   Float_t         phiJetGen[100];   //[nJetGen]
   Float_t         pxMet ;
   Float_t         pyMet ;
   Float_t         pzMet ;
   Float_t         eMet  ;
   Float_t         etaMet;
   Float_t         phiMet;
   Float_t         pxMetGen ;
   Float_t         pyMetGen ;
   Float_t         pzMetGen ;
   Float_t         eMetGen  ;
   Float_t         etaMetGen;
   Float_t         phiMetGen;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_nMC;   //!
   TBranch        *b_pdgIdMC;   //!
   TBranch        *b_statusMC;   //!
   TBranch        *b_motherIDMC;   //!
   TBranch        *b_pxMC ;   //!
   TBranch        *b_pyMC ;   //!
   TBranch        *b_pzMC ;   //!
   TBranch        *b_eMC  ;   //!
   TBranch        *b_etaMC;   //!
   TBranch        *b_phiMC;   //!
   TBranch        *b_nPhot;   //!
   TBranch        *b_pxPhot ;   //!
   TBranch        *b_pyPhot ;   //!
   TBranch        *b_pzPhot ;   //!
   TBranch        *b_ePhot  ;   //!
   TBranch        *b_etaPhot;   //!
   TBranch        *b_phiPhot;   //!
   TBranch        *b_ptiso02Phot;   //!
   TBranch        *b_ntrkiso02Phot;   //!
   TBranch        *b_ptiso025Phot;   //!
   TBranch        *b_ntrkiso025Phot;   //!
   TBranch        *b_ptiso03Phot;   //!
   TBranch        *b_ntrkiso03Phot;   //!
   TBranch        *b_ptiso035Phot;   //!
   TBranch        *b_ntrkiso035Phot;   //!
   TBranch        *b_ptiso04Phot;   //!
   TBranch        *b_ntrkiso04Phot;   //!
   TBranch        *b_ptisoatecal02Phot;   //!
   TBranch        *b_ntrkisoatecal02Phot;   //!
   TBranch        *b_ptisoatecal025Phot;   //!
   TBranch        *b_ntrkisoatecal025Phot;   //!
   TBranch        *b_ptisoatecal03Phot;   //!
   TBranch        *b_ntrkisoatecal03Phot;   //!
   TBranch        *b_ptisoatecal035Phot;   //!
   TBranch        *b_ntrkisoatecal035Phot;   //!
   TBranch        *b_ptisoatecal04Phot;   //!
   TBranch        *b_hcalovecal01Phot;   //!
   TBranch        *b_hcalovecal02Phot;   //!
   TBranch        *b_hcalovecal025Phot;   //!
   TBranch        *b_hcalovecal03Phot;   //!
   TBranch        *b_hcalovecal04Phot;   //!
   TBranch        *b_hcalovecal05Phot;   //!
   TBranch        *b_ecaliso02Phot;   //!
   TBranch        *b_ecaliso03Phot;   //!
   TBranch        *b_ecaliso035Phot;   //!
   TBranch        *b_ecaliso04Phot;   //!
   TBranch        *b_ecaliso05Phot;   //!
    TBranch        *b_LATPhot;   //!
   TBranch        *b_sMajMajPhot;   //!
   TBranch        *b_sMinMinPhot;   //!
   TBranch        *b_FisherPhot;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_pxJet ;   //!
   TBranch        *b_pyJet ;   //!
   TBranch        *b_pzJet ;   //!
   TBranch        *b_eJet  ;   //!
   TBranch        *b_etaJet;   //!
   TBranch        *b_phiJet;   //!
   TBranch        *b_nJetGen;   //!
   TBranch        *b_pxJetGen ;   //!
   TBranch        *b_pyJetGen ;   //!
   TBranch        *b_pzJetGen ;   //!
   TBranch        *b_eJetGen  ;   //!
   TBranch        *b_etaJetGen;   //!
   TBranch        *b_phiJetGen;   //!
   TBranch        *b_pxMet;   //!
   TBranch        *b_pyMet;   //!
   TBranch        *b_pzMet;   //!
   TBranch        *b_eMet;   //!
   TBranch        *b_etaMet;   //!
   TBranch        *b_phiMet;   //!
   TBranch        *b_pxMetGen;   //!
   TBranch        *b_pyMetGen;   //!
   TBranch        *b_pzMetGen;   //!
   TBranch        *b_eMetGen;   //!
   TBranch        *b_etaMetGen;   //!
   TBranch        *b_phiMetGen;   //!

   TH1D*    allptphot;
   TH1D*    allptjet;
   TH1D*    alldeltaphi;
   TH1D*    alldeltar;
   TH1D*    alldeltae;
   TH2D*    alldeltardeltae;

   TH1D*    mcptphot;
   TH1D*    mcetaphot;
   TH1D*    mcidphot;
   TH1D*    mcdeltaphi;
   TH1D*    mcdeltapt;
   TH1D*    mcptjet;
   TH1D*    mcetajet;
   TH1D*    mcidjet;

   TH1D*    gendeltaphi;
   TH1D*    gendeltapt;
   TH1D*    genmcdeltaphi;
   TH1D*    genmcdeltapt;
   TH1D*    genptjet;
   TH1D*    genetajet;

   TH1D*    leadptphot;
   TH1D*    leadetaphot;
   TH1D*    leaddeltaphi;
   TH1D*    leadptjet;
   TH1D*    leadetajet;

   TH1D*    secptphot;
   TH1D*    secetaphot;
   TH1D*    secdeltaphiphot;
   TH1D*    secptjet;
   TH1D*    secetajet;
   TH1D*    secdeltaphijet;

   TH1D*    matchdeltar;
   TH1D*    matchdeltae;

   photoncandidate phot;
   jetcandidate jet;

   // Configuration settings
   float _rcone;
   float _backtoback;
   //
   float _photetacut;
   float _photptcut;
   float _jetetacut;
   float _jet2_minpt;
   float _jet2_maxfrac;
   float _deltaphi;
   float _deltaeta;

   photonidcuts _looseid;
   photonidcuts _mediumid;
   photonidcuts _tightid;

   int irecphot;
   int isubjet; // substitute photon (reco jet)
   int isubphot; // substitute photon (parton jet)
   int i2ndphot;
   int igenphot;
   int imcphot;

   int irecjet;
   int i2ndjet;
   int igenjet;
   int imcjet;

   Bool_t matched;
   Float_t pt2ndjet;
   Float_t pt2ndphot;

   Int_t ntrkisos;
   Float_t pts;
   Float_t ptisos;
   Float_t ptisoatecals;
   Float_t hcalovecals;
   Float_t sMajMajs;
   Float_t sMinMins;
   Float_t ecalisos;

   Int_t ntrkisob;
   Float_t ptb;
   Float_t ptisob;
   Float_t ptisoatecalb;
   Float_t hcalovecalb;
   Float_t sMajMajb;
   Float_t sMinMinb;
   Float_t ecalisob;

   Int_t photonid;
   Int_t isgamma;
   Int_t ismatched;
   Float_t weight;
   Float_t nniso;
   Float_t nniso_int;
   Float_t ptph;
   Float_t ptj;
   Float_t etaph;
   Float_t etaj;
   Float_t phiph;
   Float_t phij;
   Float_t pt2jet;
   Float_t pt2sum;
   Float_t pt2vecsum;
   Float_t pt2phot;
   Float_t ptphottrue;
   Float_t ptjettrue;
   Float_t ptquarktrue;
   Float_t phiphottrue;
   Float_t phijettrue;
   Float_t phiquarktrue;


   TFile* hOutputFile ;
   TTree * S_tree ;
   TTree * B_tree ;
   TTree * ana_tree ;

   int npb;

   gjettree(TTree *tree=0, char * outputname = "test.root");
   virtual ~gjettree();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(double cross = 1., bool isgamma = true,
			 int NEVT = 10000000);
   virtual void     BookHistos();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   
   // Load configuration from the file determined by 'CAFE_CONFIG'
   // Set the environment variable 'CAFE_CONFIG' like this:
   //   gSystem->Setenv("CAFE_CONFIG", "gjettree.config");
   // Input parameter "tree" gives the name of the configuration.
   // Set the parameters like this:
   //   tree.DeltaPhi: 3.0
   //   tree.DeltaEta: 1.0
   // Remember to load the Config class before gjettree.C:
   //   gROOT->ProcessLine(".L Config.cpp++");
   void Configure(const char* cfg = "tree");
   void DumpConfig();

   //  calculate phi1-phi2 keeping value between 0 and pi
   inline float delta_phi(float phi1, float phi2);
   inline float delta_eta(float eta1, float eta2);
   //bool looseID(int i);
   //bool mediumID(int i);
   //bool tightID(int i);
   bool cutID(int i, photonidcuts const& pid);
   bool NNID(int i, std::vector<IClassifierReader*> const& classreaders);

   std::map<std::string, gjettree_histos*> histos;
   void FillHistos(std::string dirname);
   void WriteHistos();
};

#endif

#ifdef gjettree_cxx
gjettree::gjettree(TTree *tree, char * outputname)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("output.root");
      if (!f) {
         f = new TFile("output.root");
      }
      tree = (TTree*)gDirectory->Get("pippo");

   }
   Init(tree);
   hOutputFile   = new TFile(outputname , "RECREATE" ) ;

   Configure();
   DumpConfig();
}

gjettree::~gjettree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   S_tree->Write();
   B_tree->Write();
   ana_tree->Write();
   hOutputFile->Write() ;
   hOutputFile->Close() ;
   hOutputFile->Delete();
}

Int_t gjettree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t gjettree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void gjettree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normaly not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("nMC", &nMC, &b_nMC);
   fChain->SetBranchAddress("pdgIdMC", pdgIdMC, &b_pdgIdMC);
   fChain->SetBranchAddress("statusMC", statusMC, &b_statusMC);
   fChain->SetBranchAddress("motherIDMC", motherIDMC, &b_motherIDMC);
   fChain->SetBranchAddress("pxMC ", pxMC , &b_pxMC );
   fChain->SetBranchAddress("pyMC ", pyMC , &b_pyMC );
   fChain->SetBranchAddress("pzMC ", pzMC , &b_pzMC );
   fChain->SetBranchAddress("eMC  ", eMC  , &b_eMC  );
   fChain->SetBranchAddress("etaMC", etaMC, &b_etaMC);
   fChain->SetBranchAddress("phiMC", phiMC, &b_phiMC);
   fChain->SetBranchAddress("nPhot", &nPhot, &b_nPhot);
   fChain->SetBranchAddress("pxPhot ", pxPhot , &b_pxPhot );
   fChain->SetBranchAddress("pyPhot ", pyPhot , &b_pyPhot );
   fChain->SetBranchAddress("pzPhot ", pzPhot , &b_pzPhot );
   fChain->SetBranchAddress("ePhot  ", ePhot  , &b_ePhot  );
   fChain->SetBranchAddress("etaPhot", etaPhot, &b_etaPhot);
   fChain->SetBranchAddress("phiPhot", phiPhot, &b_phiPhot);
   fChain->SetBranchAddress("ptiso02Phot", ptiso02Phot, &b_ptiso02Phot);
   fChain->SetBranchAddress("ntrkiso02Phot", ntrkiso02Phot, &b_ntrkiso02Phot);
   fChain->SetBranchAddress("ptiso025Phot", ptiso025Phot, &b_ptiso025Phot);
   fChain->SetBranchAddress("ntrkiso025Phot", ntrkiso025Phot, &b_ntrkiso025Phot);
   fChain->SetBranchAddress("ptiso03Phot", ptiso03Phot, &b_ptiso03Phot);
   fChain->SetBranchAddress("ntrkiso03Phot", ntrkiso03Phot, &b_ntrkiso03Phot);
   fChain->SetBranchAddress("ptiso035Phot", ptiso035Phot, &b_ptiso035Phot);
   fChain->SetBranchAddress("ntrkiso035Phot", ntrkiso035Phot, &b_ntrkiso035Phot);
   fChain->SetBranchAddress("ptiso04Phot", ptiso04Phot, &b_ptiso04Phot);
   fChain->SetBranchAddress("ntrkiso04Phot", ntrkiso04Phot, &b_ntrkiso04Phot);
   fChain->SetBranchAddress("ptisoatecal02Phot", ptisoatecal02Phot, &b_ptisoatecal02Phot);
   fChain->SetBranchAddress("ntrkisoatecal02Phot", ntrkisoatecal02Phot, &b_ntrkisoatecal02Phot);
   fChain->SetBranchAddress("ptisoatecal025Phot", ptisoatecal025Phot, &b_ptisoatecal025Phot);
   fChain->SetBranchAddress("ntrkisoatecal025Phot", ntrkisoatecal025Phot, &b_ntrkisoatecal025Phot);
   fChain->SetBranchAddress("ptisoatecal03Phot", ptisoatecal03Phot, &b_ptisoatecal03Phot);
   fChain->SetBranchAddress("ntrkisoatecal03Phot", ntrkisoatecal03Phot, &b_ntrkisoatecal03Phot);
   fChain->SetBranchAddress("ptisoatecal035Phot", ptisoatecal035Phot, &b_ptisoatecal035Phot);
   fChain->SetBranchAddress("ntrkisoatecal035Phot", ntrkisoatecal035Phot, &b_ntrkisoatecal035Phot);
   fChain->SetBranchAddress("ptisoatecal04Phot", ptisoatecal04Phot, &b_ptisoatecal04Phot);
   fChain->SetBranchAddress("hcalovecal01Phot", hcalovecal01Phot, &b_hcalovecal01Phot);
   fChain->SetBranchAddress("hcalovecal02Phot", hcalovecal02Phot, &b_hcalovecal02Phot);
   fChain->SetBranchAddress("hcalovecal025Phot", hcalovecal025Phot, &b_hcalovecal025Phot);
   fChain->SetBranchAddress("hcalovecal03Phot", hcalovecal03Phot, &b_hcalovecal03Phot);
   fChain->SetBranchAddress("hcalovecal04Phot", hcalovecal04Phot, &b_hcalovecal04Phot);
   fChain->SetBranchAddress("hcalovecal05Phot", hcalovecal05Phot, &b_hcalovecal05Phot);
   fChain->SetBranchAddress("ecaliso02Phot", ecaliso02Phot, &b_ecaliso02Phot);
   fChain->SetBranchAddress("ecaliso03Phot", ecaliso03Phot, &b_ecaliso03Phot);
   fChain->SetBranchAddress("ecaliso035Phot", ecaliso035Phot, &b_ecaliso035Phot);
   fChain->SetBranchAddress("ecaliso04Phot", ecaliso04Phot, &b_ecaliso04Phot);
   fChain->SetBranchAddress("ecaliso05Phot", ecaliso05Phot, &b_ecaliso05Phot);
   fChain->SetBranchAddress("LATPhot", LATPhot, &b_LATPhot);
   fChain->SetBranchAddress("sMajMajPhot", sMajMajPhot, &b_sMajMajPhot);
   fChain->SetBranchAddress("sMinMinPhot", sMinMinPhot, &b_sMinMinPhot);
   fChain->SetBranchAddress("FisherPhot", FisherPhot, &b_FisherPhot);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("pxJet ", pxJet , &b_pxJet );
   fChain->SetBranchAddress("pyJet ", pyJet , &b_pyJet );
   fChain->SetBranchAddress("pzJet ", pzJet , &b_pzJet );
   fChain->SetBranchAddress("eJet  ", eJet  , &b_eJet  );
   fChain->SetBranchAddress("etaJet", etaJet, &b_etaJet);
   fChain->SetBranchAddress("phiJet", phiJet, &b_phiJet);
   fChain->SetBranchAddress("nJetGen", &nJetGen, &b_nJetGen);
   fChain->SetBranchAddress("pxJetGen ", pxJetGen , &b_pxJetGen );
   fChain->SetBranchAddress("pyJetGen ", pyJetGen , &b_pyJetGen );
   fChain->SetBranchAddress("pzJetGen ", pzJetGen , &b_pzJetGen );
   fChain->SetBranchAddress("eJetGen  ", eJetGen  , &b_eJetGen  );
   fChain->SetBranchAddress("etaJetGen", etaJetGen, &b_etaJetGen);
   fChain->SetBranchAddress("phiJetGen", phiJetGen, &b_phiJetGen);
   fChain->SetBranchAddress("pxMet ", &pxMet , &b_pxMet);
   fChain->SetBranchAddress("pyMet ", &pyMet , &b_pyMet);
   fChain->SetBranchAddress("pzMet ", &pzMet , &b_pzMet);
   fChain->SetBranchAddress("eMet  ", &eMet  , &b_eMet);
   fChain->SetBranchAddress("etaMet", &etaMet, &b_etaMet);
   fChain->SetBranchAddress("phiMet", &phiMet, &b_phiMet);
   fChain->SetBranchAddress("pxMetGen ", &pxMetGen , &b_pxMetGen);
   fChain->SetBranchAddress("pyMetGen ", &pyMetGen , &b_pyMetGen);
   fChain->SetBranchAddress("pzMetGen ", &pzMetGen , &b_pzMetGen);
   fChain->SetBranchAddress("eMetGen  ", &eMetGen  , &b_eMetGen);
   fChain->SetBranchAddress("etaMetGen", &etaMetGen, &b_etaMetGen);
   fChain->SetBranchAddress("phiMetGen", &phiMetGen, &b_phiMetGen);
   Notify();

   npb = 1000;

}

Bool_t gjettree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gjettree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

#endif // #ifdef gjettree_cxx
