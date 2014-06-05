#
# $Id: dump_DATA_53X.py,v 1.4 2013/06/18 17:19:05 meridian Exp $
#
#  configuration to dump ntuples in data
#   the only diff should be for jetmet corrections
#
import FWCore.ParameterSet.Config as cms
import sys
import os
import imp

# read in process from file
filename  = 'common_dump_config.py'
handle = open(filename, 'r')
cfo = imp.load_source("pycfg", filename, handle)
process = cfo.process
handle.close()

import HLTrigger.HLTfilters.hltHighLevel_cfi
process.DiPhotonHltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.DiPhotonHltFilter.throw = cms.bool(False)
process.DiPhotonHltFilter.HLTPaths = ["HLT_Photon*_Photon*"]
process.p = cms.Path( process.DiPhotonHltFilter* process.analysisSequence )

## DO NOT CHANGE THE PATH HERE! New modules should be added ONLY in the common configuration 
#  only paramaters should be changes for data and MC
process.source.fileNames = cms.untracked.vstring(
    #'file:/tmp/capalmer/electronevents//tmp/testPhotonAOD.root'
    #    '/store/data/Run2012A/Photon/AOD/PromptReco-v1/000/191/226/9E7EF5CF-DA87-E111-8BC6-5404A63886C6.root'
    #    '/store/data/Run2012A/Photon/RECO/PromptReco-v1/000/191/226/78D726A3-F887-E111-9A04-001D09F27003.root'
    #    '/store/data/Run2012C/DoublePhoton/AOD/PromptReco-v2/000/202/016/16610E56-FBF4-E111-A6A6-001D09F28EA3.root'
    #    'file:/afs/cern.ch/user/m/meridian/work/CMSSW532p4/src/PickEvents/highMetEvents_2012C.root'
    #'/store/data/Run2012C/DoubleElectron/AOD/22Jan2013-v1/20002/F6ABB3A5-E968-E211-ABEF-003048FFD7BE.root'#
    '/store/data/Run2012C/DoubleElectron/AOD/22Jan2013-v1/20001/F8980718-4968-E211-907F-00261894393F.root'
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Global tag
#process.GlobalTag.globaltag = cms.string("FT_R_53_V6::All")
process.GlobalTag.globaltag = cms.string("FT_53_V21_AN4::All")
#process.load("CondCore.DBCommon.CondDBCommon_cfi")
#from CondCore.DBCommon.CondDBSetup_cfi import *
#process.jec = cms.ESSource("PoolDBESSource",
#                           DBParameters = cms.PSet(
#    messageLevel = cms.untracked.int32(0)
#    ),
#                           timetype = cms.string('runnumber'),
#                           toGet = cms.VPSet(
#    cms.PSet(
#    record = cms.string('JetCorrectionsRecord'),
#    tag    = cms.string('JetCorrectorParametersCollection_Summer12_V7_DATA_AK5PF'),
#    label  = cms.untracked.string('AK5PF')
#    ),
#    ## here you add as many jet types as you need
#    ## note that the tag name is specific for the particular sqlite file
#    ),
#                           connect = cms.string('sqlite:Summer12_V7_DATA.db')
#                           )
### add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
#process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

process.ak5PFJetsCorrected = process.ak5PFJetsL1FastL2L3Residual.clone()

process.myanalysis.JetCorrectionService_akt5 = cms.string('ak5CaloL1FastL2L3Residual')
process.myanalysis.JetCorrectionService_akt7 = cms.string('ak7CaloL1FastL2L3Residual')
process.myanalysis.JetCorrectionService_jptak5 = cms.string('ak5JPTL1L2L3Residual')
process.myanalysis.JetCorrectionService_jptak7 = cms.string('ak7JPTL1L2L3Residual')
process.myanalysis.JetCorrectionService_pfakt5 = cms.string('ak5PFL1FastL2L3Residual')
process.myanalysis.JetCorrectionService_pfakt7 = cms.string('ak7PFL1FastL2L3Residual')




