#
# $Id: dump_DATA_DiPhotonSkim_42X.py,v 1.2 2011/05/22 12:57:08 meridian Exp $
#
#  configuration to dump ntuples in data
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

#Switch for diPhoton Skim

diPhotonSkim=True
skim2010=True
if (diPhotonSkim!=True):
    process.p = cms.Path(process.analysisSequence)
else:
    process.load('Configuration.Skimming.PDWG_DiPhoton_SD_cff')
    process.CaloIdIsoPath = cms.Path( process.CaloIdIsoPhotonPairsFilter * process.analysisSequence)
    process.R9IdPath = cms.Path( process.R9IdPhotonPairsFilter * process.analysisSequence)
    if (skim2010==False):
        process.schedule = cms.Schedule(process.CaloIdIsoPath, process.R9IdPath)
    else:
        process.DiPhotonHltFilter.HLTPaths = ["HLT_DoublePhoton*"]
        process.hltDiPhotonCaloIdIsoObjectProducer.triggerName = cms.string("HLT_DoublePhoton.*")
        process.schedule = cms.Schedule(process.CaloIdIsoPath)
                                   

## DO NOT CHANGE THE PATH HERE! New modules should be added ONLY in the common configuration 
#  only paramaters should be changes for data and MC
process.source.fileNames = cms.untracked.vstring(
    'file:/cmshome/meridian/data/Photon_42_2010ReReco_AOD.root'
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
)

# Global tag
# Remember to change it and adapt to your needs #
process.GlobalTag.globaltag = cms.string('GR_R_42_V12::All')

##  apply only L2 and L3 jet corrections in MC
process.metMuonJESCorAK5.corrector = "ak5CaloL2L3"

process.myanalysis.JetCorrectionService_akt5 = cms.string('ak5CaloL1L2L3')
process.myanalysis.JetCorrectionService_akt7 = cms.string('ak7CaloL1L2L3')
process.myanalysis.JetCorrectionService_jptak5 = cms.string('ak5JPTL1L2L3')
process.myanalysis.JetCorrectionService_jptak7 = cms.string('ak7JPTL1L2L3')
process.myanalysis.JetCorrectionService_pfakt5 = cms.string('ak5PFL1FastL2L3')
process.myanalysis.JetCorrectionService_pfakt7 = cms.string('ak7PFL1FastL2L3')


#print process.dumpPython()
