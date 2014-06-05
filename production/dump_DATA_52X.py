#
# $Id: dump_DATA_52X.py,v 1.3 2012/09/11 12:31:33 meridian Exp $
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
    'file:/afs/cern.ch/work/p/pandolf/public/DoubleMu_DC0254E0-8097-E111-9EC0-003048F11112.root'
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Global tag
process.GlobalTag.globaltag = cms.string("GR_R_52_V8::All")
process.load("CondCore.DBCommon.CondDBCommon_cfi")
from CondCore.DBCommon.CondDBSetup_cfi import *
process.jec = cms.ESSource("PoolDBESSource",
                           DBParameters = cms.PSet(
    messageLevel = cms.untracked.int32(0)
    ),
                           timetype = cms.string('runnumber'),
                           toGet = cms.VPSet(
    cms.PSet(
    record = cms.string('JetCorrectionsRecord'),
    tag    = cms.string('JetCorrectorParametersCollection_Summer12_V7_DATA_AK5PF'),
    label  = cms.untracked.string('AK5PF')
    ),
    ## here you add as many jet types as you need
    ## note that the tag name is specific for the particular sqlite file
    ),
                           connect = cms.string('sqlite:Summer12_V7_DATA.db')
                           )
## add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

process.myanalysis.JetCorrectionService_akt5 = cms.string('ak5CaloL1L2L3Residual')
process.myanalysis.JetCorrectionService_akt7 = cms.string('ak7CaloL1L2L3Residual')
process.myanalysis.JetCorrectionService_jptak5 = cms.string('ak5JPTL1L2L3Residual')
process.myanalysis.JetCorrectionService_jptak7 = cms.string('ak7JPTL1L2L3Residual')
process.myanalysis.JetCorrectionService_pfakt5 = cms.string('ak5PFL1FastL2L3Residual')
process.myanalysis.JetCorrectionService_pfakt7 = cms.string('ak7PFL1FastL2L3Residual')
