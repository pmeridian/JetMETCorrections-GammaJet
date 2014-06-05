import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

process = cms.Process("PAT")
options = VarParsing ('standard')

globalTags = {"PromptReco": "GR_P_V40_AN1::All", "Aug24ReReco": "FT_53_V10_AN2::All" , "Aug06ReReco": "FT_53_V6C_AN2::All", "July13ReReco": "FT_R_53_V6::All",
              "pre526FastSim": "START53_V7F::All", "526andLaterFastSim": "START53_V7F::All", "MC": "START53_V7F::All"}
## ************************************************************
## The following line must be configured properly
## datasetType must be one of the options above in "globalTags"
## ************************************************************
datasetType = "July13ReReco"

## options for testing
#options.files='file:/cmsdata/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_PU_S10_START53_V7A_AODSIM/ECDEFDB7-AAE1-E111-B576-003048C68A88.root'
options.files='file:highMetEvents.root'
maxEvents=-1

## Print out the cfA configuration information
print "Using global tag " + globalTags[datasetType] + " selected from datasetType=" + datasetType

#-- Message Logger ------------------------------------------------------------
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
        limit = cms.untracked.int32(-1),
            reportEvery = cms.untracked.int32(100)
        )
#process.MessageLogger.suppressInfo = cms.untracked.vstring('ecallaserfilter')

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
#-- Source information ------------------------------------------------------
process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(options.files)
                            )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(maxEvents) )
#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange(
#  '190645:10-190645:110',
#)
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.Geometry.GeometryIdeal_cff")
## The geometry sequence now generates a deprecation warning
## and so has been replaced by the one above
#process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")


process.GlobalTag.globaltag = globalTags[datasetType]

## The ECAL laser calibration filter
process.load('RecoMET.METFilters.ecalLaserCorrFilter_cfi')
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
process.load('CommonTools.RecoAlgos.HBHENoiseFilter_cfi')

## The CSC beam halo tight filter ____________________________________________||
process.load('RecoMET.METAnalyzers.CSCHaloFilter_cfi')

## The HCAL laser filter _____________________________________________________||
process.load("RecoMET.METFilters.hcalLaserEventFilter_cfi")
process.hcalLaserEventFilter.vetoByRunEventNumber=cms.untracked.bool(False)
process.hcalLaserEventFilter.vetoByHBHEOccupancy=cms.untracked.bool(True)

## The ECAL dead cell trigger primitive filter _______________________________||
process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
## For AOD and RECO recommendation to use recovered rechits
process.EcalDeadCellTriggerPrimitiveFilter.tpDigiCollection = cms.InputTag("ecalTPSkimNA")

## The EE bad SuperCrystal filter ____________________________________________||
process.load('RecoMET.METFilters.eeBadScFilter_cfi')

## The Good vertices collection needed by the tracking failure filter ________||
process.goodVertices = cms.EDFilter(
      "VertexSelector",
        filter = cms.bool(False),
        src = cms.InputTag("offlinePrimaryVertices"),
        cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
      )

## The tracking failure filter _______________________________________________||
process.load('RecoMET.METFilters.trackingFailureFilter_cfi')

process.HBHENoise = cms.Path(process.HBHENoiseFilter)
process.CSCTightHalo = cms.Path(process.CSCTightHaloFilter)
process.hcalLaser = cms.Path(process.hcalLaserEventFilter) 
process.ecalDeadCell  = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter )
process.trackingFailure = cms.Path( process.goodVertices * process.trackingFailureFilter )
process.eeBadSc = cms.Path ( process.eeBadScFilter )
process.ecallaserfilter = cms.Path(process.ecalLaserCorrFilter)

