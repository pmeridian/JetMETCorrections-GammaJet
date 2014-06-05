# $Id: common_dump_config.py,v 1.12 2013/06/08 09:16:57 meridian Exp $
#
#  common configuration to dump ntuples in MC and data
#    all changes affecting the path and additional modules msut be done here
#

import FWCore.ParameterSet.Config as cms
import RecoJets.JetProducers.JetIDParams_cfi
theJetIDParams = RecoJets.JetProducers.JetIDParams_cfi.JetIDParams.clone()
from RecoJets.Configuration.RecoGenJets_cff import *
from RecoMET.Configuration.RecoGenMET_cff import *
from RecoMET.Configuration.GenMETParticles_cff import *

process = cms.Process("myprocess")

process.load("FWCore.MessageLogger.MessageLogger_cfi")

is41X=False
doCleanMet=False
metFilters=True

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("PhysicsTools.HepMCCandAlgos.genParticleCandidates_cfi")
##process.load("Configuration.Geometry.GeometryIdeal_cff")
## The geometry sequence now generates a deprecation warning
## and so has been replaced by the one above
##process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.GeometryDB_cff") 
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration/StandardSequences/Services_cff')
#process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
#process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag = cms.string('START311_V2::All')
#process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
#process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
#process.load("Geometry.TrackerGeometryBuilder.trackerGeometry_cfi")
#process.load("Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cfi")
#process.load("Geometry.CaloEventSetup.CaloTopology_cfi")

process.load('Configuration/StandardSequences/Reconstruction_cff')

#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck")
process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(
    'file:events_GluGluToHToZZTo2L2Q_M-550_7TeV-powheg-pythia6_Summer11_PROVA.root'
)

)

from HiggsAnalysis.HiggsTo2photons.hggPhotonIDCuts_cfi import *

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
)

process.MessageLogger.cerr.FwkReport.reportEvery = 10

#monster track event cleaning
process.monster = cms.EDFilter(
   "FilterOutScraping",
   applyfilter = cms.untracked.bool(True),
   debugOn = cms.untracked.bool(False),
   numtrack = cms.untracked.uint32(10),
   thresh = cms.untracked.double(0.2)
)


###########  EB SPIKE CLEANING BEGIN #####################


###########  EB SPIKE CLEANING END   #####################

## produce JPT jets
#process.load('RecoJets.Configuration.RecoJPTJets_cff')

#############   Include the corrections ##########
process.load("RecoMET.Configuration.RecoGenMET_cff")
process.load("RecoMET.Configuration.GenMETParticles_cff")
process.load("JetMETCorrections.Configuration.DefaultJEC_cff")
#process.load("JetMETCorrections.Type1MET.MetType1Corrections_cff")    # new for 52X
process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")        # chiara: controlla la storia dell'ordine

#process.ak5CaloL1Offset.useCondDB = False
#process.ak5PFL1Offset.useCondDB = False
#process.ak5JPTL1Offset.useCondDB = False

#process.ak5CaloL1Fastjet.useCondDB = False
#process.ak5PFL1Fastjet.useCondDB = False
#process.ak5JPTL1Fastjet.useCondDB = False

process.ak5PFJets.doAreaFastjet = True
#process.ak7PFJets.doAreaFastjet = True

process.load("JetMETCorrections.GammaJet.NoPileUp_cff")

# new for 52X
#process.metMuonJESCorAK5 = process.metJESCorAK5CaloJet.clone()
#process.metMuonJESCorAK5.inputUncorJetsLabel = "ak5CaloJets"
#process.metMuonJESCorAK5.corrector = "ak5CaloL2L3"
#process.metMuonJESCorAK5.inputUncorMetLabel = "corMetGlobalMuons"

# new for 52X 
#* process.producePFMETCorrections)
#from JetMETCorrections.Type1MET.MetType1Corrections_cff import metJESCorAK5PFJet
#process.metJESCorPFAK5 = metJESCorAK5PFJet.clone()
##process.metJESCorPFAK5.inputUncorJetsLabel = "pfJets"  
#process.metJESCorPFAK5.corrector = cms.string('ak5PFL1FastL2L3')
#process.metCorSequence = cms.Sequence( process.metJESCorPFAK5 )

# configure B-tagging to be run on ak5PFJets
process.myBtag = cms.Sequence(process.ak5JetTracksAssociatorAtVertex*process.btagging)
process.ak5JetTracksAssociatorAtVertex.jets = cms.InputTag("ak5PFJets")
process.softMuonTagInfos.jets =  cms.InputTag("ak5PFJets")
process.softElectronTagInfos.jets =  cms.InputTag("ak5PFJets")

# mva jetId
from CMGTools.External.pujetidsequence_cff import puJetMva

# chiara - init!
# this is the PF candidate isolation with pfnopu input and custom vetoes
from CommonTools.ParticleFlow.pfNoPileUp_cff import *
pfPileUp.PFCandidates = "particleFlow"
pfNoPileUp.bottomCollection = "particleFlow"
process.myPfPUSequence = cms.Sequence( pfPileUp * pfNoPileUp )

# here the isolations
from MyAnalysis.IsolationTools.electronPFIsolations_cff import *
from MyAnalysis.IsolationTools.muonPFIsolations_cff import *
process.myElectronPFIsoChHad03  = electronPFIsoChHad03.clone()
process.myElectronPFIsoNHad03   = electronPFIsoNHad03.clone();
process.myElectronPFIsoPhoton03 = electronPFIsoPhoton03.clone();
process.myElectronPFIsoChHad04  = electronPFIsoChHad04.clone()
process.myElectronPFIsoNHad04   = electronPFIsoNHad04.clone();
process.myElectronPFIsoPhoton04 = electronPFIsoPhoton04.clone();
process.myElectronPFIsoChHad05  = electronPFIsoChHad05.clone()
process.myElectronPFIsoNHad05   = electronPFIsoNHad05.clone();
process.myElectronPFIsoPhoton05 = electronPFIsoPhoton05.clone();
#
process.myPfIsolationSingleType = cms.Sequence(process.myElectronPFIsoChHad03 * process.myElectronPFIsoNHad03 * process.myElectronPFIsoPhoton03 * process.myElectronPFIsoChHad04 * process.myElectronPFIsoNHad04 * process.myElectronPFIsoPhoton04 * process.myElectronPFIsoChHad05 * process.myElectronPFIsoNHad05 * process.myElectronPFIsoPhoton05 );
# chiara - end!





#from HiggsAnalysis.HiggsToGammaGamma.PhotonFixParams4_2_cfi import *
## dumper module
process.myanalysis = cms.EDAnalyzer("GammaJetAnalyzer",
    debug = cms.bool(False),
#    PFParameters = PhotonFixParameters,
    outFileName = cms.untracked.string("output.root"),                                    
    dumpBeamHaloInformations = cms.untracked.bool(True),
    dumpAKT5Jets=cms.untracked.bool(True),
    dumpAKT7Jets=cms.untracked.bool(True),
    dumpJPTAKT5Jets=cms.untracked.bool(False),
    dumpJPTAKT7Jets=cms.untracked.bool(False),
    dumpPFAKT5Jets=cms.untracked.bool(True),
    dumpPFAKT7Jets=cms.untracked.bool(True),                                                                        
    dumpKT4Jets=cms.untracked.bool(False),
    dumpKT6Jets=cms.untracked.bool(False),                                    
    recoProducer = cms.string('ecalRecHit'),
    PUSummaryInfoCollection = cms.InputTag("addPileupInfo"),
    MCTruthCollection = cms.untracked.InputTag("source"),
    genMet = cms.untracked.InputTag("genMetTrue"),
    met = cms.untracked.InputTag("met"),
    tracks = cms.untracked.InputTag("generalTracks"),
    Electronsrc = cms.untracked.InputTag("gsfElectrons"),
    Photonsrc = cms.untracked.InputTag("photons"),
    Conversionsrc = cms.untracked.InputTag("allConversions"),
    recoCollection = cms.string('reducedEcalRecHitsEB'),
    JetCorrectionService_akt5 = cms.string('ak5CaloL1L2L3'),
    JetCorrectionService_akt7 = cms.string('ak7CaloL1L2L3'),
    JetCorrectionService_jptak5 = cms.string('ak5JPTL1L2L3'),
    JetCorrectionService_jptak7 = cms.string('ak7JPTL1L2L3'),
    JetCorrectionService_pfakt5 = cms.string('ak5PFL1FastL2L3'),
    JetCorrectionService_pfakt5_nopu = cms.string('ak5PFL1FastL2L3NoPU'),
    JetCorrectionService_pfakt7 = cms.string('ak7PFL1FastL2L3'),
    jetskt4 = cms.untracked.InputTag("kt4CaloJets"),
    jetskt6 = cms.untracked.InputTag("kt6CaloJets"),
    jetsakt5 = cms.untracked.InputTag("ak5CaloJets"),
    jetsakt7 = cms.untracked.InputTag("ak7CaloJets"),
    jetsjptak5 = cms.untracked.InputTag("JetPlusTrackZSPCorJetAntiKt5"),
    jetspfkt4 = cms.untracked.InputTag("kt4PFJets"),
    jetspfkt6 = cms.untracked.InputTag("kt6PFJets"),
    jetspfakt5 = cms.untracked.InputTag("ak5PFJets"),
    jetspfakt5_nopu = cms.untracked.InputTag("ak5PFJetsNoPU"),
    jetspfakt7 = cms.untracked.InputTag("ak7PFJets"),
    hbhits = cms.untracked.InputTag("hbhereco"),
    jetsgenkt4 = cms.untracked.InputTag("kt4GenJets"),
    jetsgenkt6 = cms.untracked.InputTag("kt6GenJets"),
    jetsgenakt5 = cms.untracked.InputTag("ak5GenJets"),
    jetsgenakt7 = cms.untracked.InputTag("ak7GenJets"),
    TriggerTag = cms.untracked.InputTag("TriggerResults::HLT"),
    vertices = cms.untracked.InputTag("offlinePrimaryVerticesWithBS"),
                                    
    hggPhotonIDConfiguration = cms.PSet(hggPhotonIDCuts),
    puJetIDAlgos = puJetMva.algos,
                                    
    genjetptthr = cms.double(5.),
    calojetptthr = cms.double(5.),
    pfjetptthr = cms.double(6.),
    jptjetptthr = cms.double(6.),
    genjetnmin = cms.int32(10),
    pfjetnmin = cms.int32(5),
    jptjetnmin = cms.int32(5),
    JetIDParams = theJetIDParams,
    Xsec = cms.double(1.)
                                    
)

# compute rho with PF candidates
process.load('RecoJets.Configuration.RecoPFJets_cff')
#process.kt6PFJets = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True,  doAreaFastjet = cms.bool(True), voronoiRfact = cms.double(0.9) )
#process.kt6PFJets = process.kt6PFJets.clone( rParam = 0.6, doRhoFastjet = True )

#process.load('RecoJets.JetProducers.kt6PFJets_cfi')
#process.kt6PFJetsForIso = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True,  doAreaFastjet = cms.bool(True), voronoiRfact = cms.double(0.9) )
process.kt6PFJetsForIso = process.kt6PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIso.Rho_EtaMax = cms.double(2.5)

# compute rho with calo towers
#process.load('RecoJets.JetProducers.kt6CaloJets_cfi')
#process.kt6CaloJetsForIso = process.kt4CaloJets.clone( rParam = 0.6, doRhoFastjet = True, doAreaFastjet = cms.bool(True), voronoiRfact = cms.double(0.9) )
process.kt6CaloJetsForIso = process.kt6CaloJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6CaloJetsForIso.Rho_EtaMax = cms.double(2.5)

# re-reconstructing the primary vertices with the Deterministic Annealing (DA) vertex finder
# from B. Mangano studies
if (is41X):
    from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesDA_cfi import *
    import RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesDA_cfi


# associated met producers (G. Cerminara P. Silva)
from CommonTools.ClusteredPFMetProducer.ClusteredPFMetProducer_cfi import ClusteredPFMetProducer
process.ClusteredPFMetProducerStd = ClusteredPFMetProducer.clone()

process.ak5PFJetsCorrected = process.ak5PFJetsL1FastL2L3.clone()
              
# cleaned MET (M. Marionneau)
process.load("CommonTools.ParticleFlow.PF2PAT_cff")
process.ak5PFJetsL1FastL2L3   = cms.EDProducer('PFJetCorrectionProducer',
                                               src         = cms.InputTag('pfJets'),
                                               correctors  = cms.vstring('ak5PFL1FastL2L3')
                                               )
if (doCleanMet):
    process.load("MMarionneau.CleanMETProducer.cleanMETProducer_cfi")
    process.cleanMETProducer.pfMETInput = cms.InputTag("metJESCorPFAK5")
    process.cleanMETProducer.pfJETInput = cms.InputTag("ak5PFJetsL1FastL2L3")
    process.cleanMETProducer.vertexInput = cms.InputTag("offlinePrimaryVerticesWithBS")
    process.cleanMETProducer.vtxUserDef = cms.untracked.bool(False)


if (is41X):
    print "ADDING DA VERTICES NOT DEFAULT IN 41X RELEASE"
    process.offlinePrimaryVerticesDA = RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesDA_cfi.offlinePrimaryVerticesDA.clone()
    process.offlinePrimaryVerticesDA.useBeamConstraint = cms.bool(True)
    process.offlinePrimaryVerticesDA.TkClusParameters.TkDAClusParameters.Tmin = cms.double(4.)
    process.offlinePrimaryVerticesDA.TkClusParameters.TkDAClusParameters.vertexSize = cms.double(0.01) 

# high purity tracks
process.highPurityTracks = cms.EDFilter("TrackSelector",
    src = cms.InputTag("generalTracks"),
    cut = cms.string('quality("highPurity")')
)

#process.analysisSequence = cms.Sequence( process.monster )
process.analysisSequence = cms.Sequence(  )
if (is41X):
    process.analysisSequence *= process.offlinePrimaryVerticesDA
    process.myanalysis.vertices = cms.untracked.InputTag("offlinePrimaryVerticesDA")

## The ECAL laser calibration filter
if (metFilters):
    process.load('RecoMET.METFilters.ecalLaserCorrFilter_cfi')
    process.ecalLaserCorrFilter.taggingMode = cms.bool(True)

## The good primary vertex filter ____________________________________________||
    process.primaryVertexFilter = cms.EDFilter(
        "VertexSelector",
        src = cms.InputTag("offlinePrimaryVertices"),
        cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
        filter = cms.bool(True)
        )
## The beam scraping filter __________________________________________________||
    process.noscraping = cms.EDFilter(
        "FilterOutScraping",
        applyfilter = cms.untracked.bool(True),
        debugOn = cms.untracked.bool(False),
        numtrack = cms.untracked.uint32(10),
        thresh = cms.untracked.double(0.25)
        )
    
## The iso-based HBHE noise filter ___________________________________________||
    process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
    
## The CSC beam halo tight filter ____________________________________________||
##    process.load('RecoMET.METAnalyzers.CSCHaloFilter_cfi')
    
## The HCAL laser filter _____________________________________________________||
    process.load("RecoMET.METFilters.hcalLaserEventFilter_cfi")
    process.hcalLaserEventFilter.vetoByRunEventNumber=cms.untracked.bool(False)
    process.hcalLaserEventFilter.vetoByHBHEOccupancy=cms.untracked.bool(True)
    process.hcalLaserEventFilter.taggingMode = cms.bool(True)
## The ECAL dead cell trigger primitive filter _______________________________||
    process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
## For AOD and RECO recommendation to use recovered rechits
    process.EcalDeadCellTriggerPrimitiveFilter.tpDigiCollection = cms.InputTag("ecalTPSkimNA")
    process.EcalDeadCellTriggerPrimitiveFilter.taggingMode = cms.bool(True)
## The EE bad SuperCrystal filter ____________________________________________||
    process.load('RecoMET.METFilters.eeBadScFilter_cfi')
    process.eeBadScFilter.taggingMode = cms.bool(True)

## The Good vertices collection needed by the tracking failure filter ________||
    process.goodVertices = cms.EDFilter(
        "VertexSelector",
        filter = cms.bool(False),
        src = cms.InputTag("offlinePrimaryVertices"),
        cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
        )
## The tracking failure filter _______________________________________________||
    process.load('RecoMET.METFilters.trackingFailureFilter_cfi')

    process.analysisSequence *= process.HBHENoiseFilterResultProducer
 #   process.analysisSequence *= process.CSCTightHaloFilter  
    process.analysisSequence *= process.hcalLaserEventFilter 
    process.analysisSequence *= process.EcalDeadCellTriggerPrimitiveFilter  
#    process.analysisSequence *=  process.goodVertices * process.trackingFailureFilter 
    process.analysisSequence *=  process.eeBadScFilter 
    process.analysisSequence *=  process.ecalLaserCorrFilter 




if (doCleanMet):    
    process.analysisSequence *=  (process.highPurityTracks * process.ak5PFJets *process.kt6PFJetsForIso * process.kt6CaloJetsForIso * process.myBtag * process.PF2PAT * process.ak5PFJetsL1FastL2L3 * process.metCorSequence * process.ClusteredPFMetProducerStd  * process.cleanMETProducer *  process.producePFNoPileUp * process.ak5PFJetsCorrected * process.myanalysis)
else:
    process.analysisSequence *=  (process.highPurityTracks * process.ak5PFJets * process.kt6PFJetsForIso * process.kt6CaloJetsForIso * process.myBtag * process.ClusteredPFMetProducerStd * process.pfJetMETcorr * process.pfchsMETcorr * process.pfType1CorrectedMet * process.myPfPUSequence * process.myPfIsolationSingleType * process.ak5PFJetsCorrected * process.myanalysis)

#process.out = cms.OutputModule("PoolOutputModule",
#                               fileName = cms.untracked.string("azzo.root")
#                               )
process.mypath = cms.Path(process.analysisSequence)
#process.outpath = cms.EndPath(process.out)
