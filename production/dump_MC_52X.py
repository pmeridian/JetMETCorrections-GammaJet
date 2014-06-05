#
# $Id: dump_MC_52X.py,v 1.4 2013/06/08 09:16:57 meridian Exp $
#
#  configuration to dump ntuples in MC
#   the only diff should be for jetmet corrections
#
import FWCore.ParameterSet.Config as cms
import sys
import os
import imp

# contains all details of what to run and defines the process and the path
#from  JetMETCorrections.GammaJet.dumper_process_cff import process

# read in process from file
filename  = 'common_dump_config.py'
handle = open(filename, 'r')
cfo = imp.load_source("pycfg", filename, handle)
process = cfo.process
handle.close()

process.p = cms.Path(process.analysisSequence)

## DO NOT CHANGE THE PATH HERE! New modules should be added ONLY in the common configuration 
#  only paramaters should be changes for data and MC
<<<<<<< dump_MC_52X.py
process.source.fileNames = cms.untracked.vstring('file:mcPool.root')
=======
process.source.fileNames = cms.untracked.vstring('root://pccmsrm27.cern.ch///cms/local/meridian/data/Summer12/vbf1208tev.root')
>>>>>>> 1.3

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Global tag#
# Remember to change it and adapt to your needs #
process.GlobalTag.globaltag = cms.string('START52_V9::All')

##  apply only L2 and L3 jet corrections in MC
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
    tag    = cms.string('JetCorrectorParametersCollection_Summer12_V7_MC_AK5PF'),
    label  = cms.untracked.string('AK5PF')
    ),
    ## here you add as many jet types as you need
    ## note that the tag name is specific for the particular sqlite file
    ),
                           connect = cms.string('sqlite:Summer12_V7_MC.db')
                           )
## add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

process.myanalysis.JetCorrectionService_akt5 = cms.string('ak5CaloL1L2L3')
process.myanalysis.JetCorrectionService_akt7 = cms.string('ak7CaloL1L2L3')
process.myanalysis.JetCorrectionService_jptak5 = cms.string('ak5JPTL1L2L3')
process.myanalysis.JetCorrectionService_jptak7 = cms.string('ak7JPTL1L2L3')
process.myanalysis.JetCorrectionService_pfakt5 = cms.string('ak5PFL1FastL2L3')
process.myanalysis.JetCorrectionService_pfakt7 = cms.string('ak7PFL1FastL2L3')
#process.myanalysis.TriggerTag = cms.untracked.InputTag("TriggerResults::REDIGI311X")

#print process.dumpPython()
