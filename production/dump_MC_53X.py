#
# $Id: dump_MC_53X.py,v 1.2 2013/06/14 06:31:10 meridian Exp $
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

process.source.fileNames = cms.untracked.vstring(
#'/store/group/phys_higgs/meridian/HGGProd/VBF_HToGG_M-125_8TeV-powheg-pythia6-TuneZ2Star-PU_S10_START53_V7A/meridian/VBF_HToGG_M-125_8TeV-powheg-LHE_v1/VBF_HToGG_M-125_8TeV-powheg-pythia6-TuneZ2Star-PU_S10_START53_V7A/989286ed06004706c6519652f2332a71/POWHEG_PYTHIA6_Tauola_H_2gamma_8TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_PU_65_1_Ty7.root',
#'/store/group/phys_higgs/meridian/HGGProd/VBF_HToGG_M-125_8TeV-powheg-pythia6-TuneZ2Star-PU_S10_START53_V7A/meridian/VBF_HToGG_M-125_8TeV-powheg-LHE_v1/VBF_HToGG_M-125_8TeV-powheg-pythia6-TuneZ2Star-PU_S10_START53_V7A/989286ed06004706c6519652f2332a71/POWHEG_PYTHIA6_Tauola_H_2gamma_8TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_PU_66_1_q5x.root',
#'/store/group/phys_higgs/meridian/HGGProd/VBF_HToGG_M-125_8TeV-powheg-pythia6-TuneZ2Star-PU_S10_START53_V7A/meridian/VBF_HToGG_M-125_8TeV-powheg-LHE_v1/VBF_HToGG_M-125_8TeV-powheg-pythia6-TuneZ2Star-PU_S10_START53_V7A/989286ed06004706c6519652f2332a71/POWHEG_PYTHIA6_Tauola_H_2gamma_8TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_PU_67_1_51D.root',
#'/store/group/phys_higgs/meridian/HGGProd/VBF_HToGG_M-125_8TeV-powheg-pythia6-TuneZ2Star-PU_S10_START53_V7A/meridian/VBF_HToGG_M-125_8TeV-powheg-LHE_v1/VBF_HToGG_M-125_8TeV-powheg-pythia6-TuneZ2Star-PU_S10_START53_V7A/989286ed06004706c6519652f2332a71/POWHEG_PYTHIA6_Tauola_H_2gamma_8TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_PU_6_1_wYc.root ',
#'/store/group/phys_higgs/meridian/HGGProd/VBF_HToGG_M-125_8TeV-powheg-pythia6-TuneZ2Star-PU_S10_START53_V7A/meridian/VBF_HToGG_M-125_8TeV-powheg-LHE_v1/VBF_HToGG_M-125_8TeV-powheg-pythia6-TuneZ2Star-PU_S10_START53_V7A/989286ed06004706c6519652f2332a71/POWHEG_PYTHIA6_Tauola_H_2gamma_8TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_PU_70_1_cD3.root',
#'/store/group/phys_higgs/meridian/HGGProd/VBF_HToGG_M-125_8TeV-powheg-pythia6-TuneZ2Star-PU_S10_START53_V7A/meridian/VBF_HToGG_M-125_8TeV-powheg-LHE_v1/VBF_HToGG_M-125_8TeV-powheg-pythia6-TuneZ2Star-PU_S10_START53_V7A/989286ed06004706c6519652f2332a71/POWHEG_PYTHIA6_Tauola_H_2gamma_8TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_PU_71_1_5D1.root'
    "file:/afs/cern.ch/work/p/pandolf/public/events_GluGluHToGG_M-125_8TeV-pythia6_Summer12_DR53X.root"
#    "root://cmseos.fnal.gov///eos/uscms/store/user/lpcpjm/PrivateMC/FastSim/539p3/ewkinoHiggs/Matched/hW/aaw/gen/WinoNLSP_chargino350_bino1_5_hw_aaw_EV20051.root"
#    '/store/group/phys_higgs/fmargaro/CRAB_thq_leptonic_test_LHE/tHqLeptonic_mH125_8TeV_testtest_SIM/b6a3cff37c57cbbd32fc33adf9268b96/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_tauola_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_PU_100_1_nHJ.root'
    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.myanalysis.outFileName = cms.untracked.string("output.root")
# Global tag#
# Remember to change it and adapt to your needs #
process.GlobalTag.globaltag = cms.string('START53_V23::All')

##  apply only L2 and L3 jet corrections in MC
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
#    tag    = cms.string('JetCorrectorParametersCollection_Summer12_V7_MC_AK5PF'),
#    label  = cms.untracked.string('AK5PF')
#    ),
#    ## here you add as many jet types as you need
#    ## note that the tag name is specific for the particular sqlite file
#    ),
#                           connect = cms.string('sqlite:Summer12_V7_MC.db')
#                           )
### add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
#process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

process.myanalysis.JetCorrectionService_akt5 = cms.string('ak5CaloL1FastL2L3')
process.myanalysis.JetCorrectionService_akt7 = cms.string('ak7CaloL1FastL2L3')
process.myanalysis.JetCorrectionService_jptak5 = cms.string('ak5JPTL1L2L3')
process.myanalysis.JetCorrectionService_jptak7 = cms.string('ak7JPTL1L2L3')
process.myanalysis.JetCorrectionService_pfakt5 = cms.string('ak5PFL1FastL2L3')
process.myanalysis.JetCorrectionService_pfakt7 = cms.string('ak7PFL1FastL2L3')
#process.myanalysis.TriggerTag = cms.untracked.InputTag("TriggerResults::REDIGI311X")

#print process.dumpPython()
