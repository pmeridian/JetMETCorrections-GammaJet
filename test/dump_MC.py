#
# $Id: $
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

## DO NOT CHANGE THE PATH HERE! New modules should be added ONLY in the common configuration 
#  only paramaters should be changes for data and MC
process.source.fileNames = cms.untracked.vstring(
'/store/mc/Spring11/GluGluToHToGG_M-100_7TeV-powheg-pythia6/AODSIM/PU_S1_START311_V1G1-v1/0000/FE78A3B2-B04F-E011-95D7-0025B3E0638E.root'
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

# Global tag
process.GlobalTag.globaltag = cms.string('START311_V2::All')

##  apply only L2 and L3 jet corrections in MC
process.metMuonJESCorAK5.corrector = "ak5CaloL2L3"

process.myanalysis.JetCorrectionService_akt5 = cms.string('ak5CaloL1L2L3')
process.myanalysis.JetCorrectionService_akt7 = cms.string('ak7CaloL2L3')
process.myanalysis.JetCorrectionService_jptak5 = cms.string('ak5JPTL1L2L3')
process.myanalysis.JetCorrectionService_jptak7 = cms.string('ak7JPTL2L3')
process.myanalysis.JetCorrectionService_pfakt5 = cms.string('ak5PFL1L2L3')
process.myanalysis.JetCorrectionService_pfakt7 = cms.string('ak7PFL2L3')

print process.dumpPython()
