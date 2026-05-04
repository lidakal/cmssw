### HiForest Configuration
# Input: miniAOD
# Type: data

import FWCore.ParameterSet.Config as cms
#Switch the era to nbe Era_Run3_2023_cff - Run3_2023
from Configuration.Eras.Era_Run3_2023_cff import Run3_2023
process = cms.Process('HiForest',Run3_2023)

###############################################################################

# HiForest info
process.load("HeavyIonsAnalysis.EventAnalysis.HiForestInfo_cfi")
process.HiForestInfo.info = cms.vstring("HiForest, miniAOD, 132X, data")

# import subprocess, os
# version = subprocess.check_output(
#     ['git', '-C', os.path.expandvars('$CMSSW_BASE/src'), 'describe', '--tags'])
# if version == '':
#     version = 'no git info'
# process.HiForestInfo.HiForestVersion = cms.string(version)

###############################################################################

# input files
process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring(
        # 2023 data with prompt reco
        '/store/hidata/HIRun2023A/HIForward0/MINIAOD/PromptReco-v2/000/374/810/00000/563ffd45-c6d0-45d1-adc2-46424283e330.root'

        # 2023 data with 2024 upc rereco
        # 'root://cms-xrd-global.cern.ch//store/hidata/HIRun2023A/HIForward0/MINIAOD/16Jan2024-v1/2810000/0640e99d-84f5-49fd-a012-48be69252e99.root'
    ), 
)

# number of events to process, set to -1 to process all events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

import FWCore.PythonUtilities.LumiList as LumiList
process.source.lumisToProcess = LumiList.LumiList(filename = 'Cert_Collisions2023HI_374288_375823_Good_ZDC_Golden.json').getVLuminosityBlockRange()

# Select specific event
# process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange('1:172')
# process.source.eventsToProcess = cms.untracked.VEventRange('1:181724455')


###############################################################################

# load Global Tag, geometry, etc.
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100


from Configuration.AlCa.GlobalTag import GlobalTag
#SWITCHING THE GT TO PROMPT RECO of PbPb
process.GlobalTag = GlobalTag(process.GlobalTag, '132X_dataRun3_Prompt_v4', '')
process.HiForestInfo.GlobalTagLabel = process.GlobalTag.globaltag

###############################################################################

# root output
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("HiForestMiniAOD.root"))

# # edm output for debugging purposes
# process.output = cms.OutputModule(
#     "PoolOutputModule",
#     fileName = cms.untracked.string('HiForestEDM.root'),
#     outputCommands = cms.untracked.vstring(
#         'keep *',
#         )
#     )

# process.output_path = cms.EndPath(process.output)

###############################################################################

# event analysis
#CM TEMP EDIT
process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_data_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.skimanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hltobject_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.l1object_cfi')

process.hiEvtAnalyzer.doCentrality = cms.bool(False)
process.hiEvtAnalyzer.doHFfilters = cms.bool(False)

# from HeavyIonsAnalysis.EventAnalysis.hltobject_cfi import trigger_list_data
# process.hltobject.triggerNames = trigger_list_data
# process.hltobject.triggerNames = cms.vstring(['HLT_HIUPC_SingleJet8_ZDC1nXOR_MaxPixelCluster50000_v2'])

process.load('HeavyIonsAnalysis.EventAnalysis.particleFlowAnalyser_cfi')
################################
# electrons, photons, muons
process.load('HeavyIonsAnalysis.EGMAnalysis.ggHiNtuplizer_cfi')
process.ggHiNtuplizer.doMuons = cms.bool(False)
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
################################
# jet reco sequence
#CM Edit - switch to ak4PFJetSeequence
#process.load('HeavyIonsAnalysis.JetAnalysis.akCs4PFJetSequence_pponPbPb_data_cff')
process.load("HeavyIonsAnalysis.JetAnalysis.ak2PFJetSequence_ppref_data_cff")
process.load("HeavyIonsAnalysis.JetAnalysis.ak3PFJetSequence_ppref_data_cff")
process.load("HeavyIonsAnalysis.JetAnalysis.ak4PFJetSequence_ppref_data_cff")
process.load('HeavyIonsAnalysis.JetAnalysis.ak4CaloJetSequence_pp_data_cff')


#The following series of analyzers is to hack in a calorimeter jet correction
process.hltAK4CaloRelativeCorrector = cms.EDProducer("LXXXCorrectorProducer",
                                                     algorithm = cms.string('AK4Calo'),
                                                     level = cms.string('L2Relative')
)
process.hltAK4CaloAbsoluteCorrector = cms.EDProducer("LXXXCorrectorProducer",
                                                     algorithm = cms.string('AK4Calo'),
                                                     level = cms.string('L3Absolute')
                                                 )
process.hltAK4CaloCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
                                             correctors = cms.VInputTag("hltAK4CaloRelativeCorrector", "hltAK4CaloAbsoluteCorrector")
)
process.hltAK4CaloJetsCorrected = cms.EDProducer("CorrectedCaloJetProducer",
                                                 correctors = cms.VInputTag("hltAK4CaloCorrector"),
                                                 src = cms.InputTag("slimmedCaloJets")
)

process.ak4CaloJetAnalyzer.jetTag = cms.InputTag("hltAK4CaloJetsCorrected")
#End calorimeter jet correction hack


################################
# tracks
process.load("HeavyIonsAnalysis.TrackAnalysis.TrackAnalyzers_cff")
#process.load("HeavyIonsAnalysis.MuonAnalysis.unpackedMuons_cfi")
#process.load("HeavyIonsAnalysis.MuonAnalysis.muonAnalyzer_cfi")
###############################################################################

# ZDC RecHit Producer
process.load('HeavyIonsAnalysis.ZDCAnalysis.QWZDC2018Producer_cfi')
process.load('HeavyIonsAnalysis.ZDCAnalysis.QWZDC2018RecHit_cfi')
process.load('HeavyIonsAnalysis.ZDCAnalysis.zdcanalyzer_cfi')

process.zdcanalyzer.doZDCRecHit = False
process.zdcanalyzer.doZDCDigi = True
process.zdcanalyzer.zdcRecHitSrc = cms.InputTag("QWzdcreco")
process.zdcanalyzer.zdcDigiSrc = cms.InputTag("hcalDigis", "ZDC")
process.zdcanalyzer.calZDCDigi = False
process.zdcanalyzer.verbose = False

from CondCore.CondDB.CondDB_cfi import *
process.es_pool = cms.ESSource("PoolDBESSource",
   timetype = cms.string('runnumber'),
   toGet = cms.VPSet(
       cms.PSet(
           record = cms.string("HcalElectronicsMapRcd"),
           tag = cms.string("HcalElectronicsMap_2021_v2.0_data")
       )
   ),
   connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
       authenticationMethod = cms.untracked.uint32(1)
   )

process.es_prefer = cms.ESPrefer('HcalTextCalibrations', 'es_ascii')
process.es_ascii = cms.ESSource(
   'HcalTextCalibrations',
   input = cms.VPSet(
       cms.PSet(

           object = cms.string('ElectronicsMap'),
           file = cms.FileInPath("HeavyIonsAnalysis/Configuration/data/emap_2023_newZDC_v3.txt")

            )
       )
   )

process.zdc = cms.Sequence(process.zdcdigi * process.QWzdcreco * process.zdcanalyzer)

###############################################################################
# main forest sequence
process.forest = cms.Path(
    process.HiForestInfo 
    + process.hiEvtAnalyzer
    + process.hltanalysis
    # + process.hltobject
    # + process.l1object
    # + process.trackSequencePP
    # + process.hltAK4CaloRelativeCorrector 
    # + process.hltAK4CaloAbsoluteCorrector
    # + process.hltAK4CaloCorrector
    # + process.hltAK4CaloJetsCorrected
    # + process.ak4CaloJetAnalyzer
    # + process.particleFlowAnalyser
    # + process.ggHiNtuplizer
    + process.zdc
    # + process.unpackedMuons
    # + process.muonAnalyzer
    )

#customisation
##########################
## Event Selection -> add the needed filters here
##########################
#
process.load('HeavyIonsAnalysis.EventAnalysis.collisionEventSelection_cff')
# process.pclusterCompatibilityFilter = cms.Path(process.clusterCompatibilityFilter)
process.pprimaryVertexFilter = cms.Path(process.primaryVertexFilter)
# process.load('HeavyIonsAnalysis.EventAnalysis.hffilter_cfi')
# process.pphfCoincFilter4Th2 = cms.Path(process.phfCoincFilter4Th2)
# process.pphfCoincFilter1Th3 = cms.Path(process.phfCoincFilter1Th3)
# process.pphfCoincFilter2Th3 = cms.Path(process.phfCoincFilter2Th3)
# process.pphfCoincFilter3Th3 = cms.Path(process.phfCoincFilter3Th3)
# process.pphfCoincFilter4Th3 = cms.Path(process.phfCoincFilter4Th3)
# process.pphfCoincFilter5Th3 = cms.Path(process.phfCoincFilter5Th3)
# process.pphfCoincFilter1Th4 = cms.Path(process.phfCoincFilter1Th4)
# process.pphfCoincFilter2Th4 = cms.Path(process.phfCoincFilter2Th4)
# process.pphfCoincFilter3Th4 = cms.Path(process.phfCoincFilter3Th4)
# process.pphfCoincFilter4Th4 = cms.Path(process.phfCoincFilter4Th4)
# process.pphfCoincFilter5Th4 = cms.Path(process.phfCoincFilter5Th4)
# process.pphfCoincFilter1Th5 = cms.Path(process.phfCoincFilter1Th5)
# process.pphfCoincFilter2Th5 = cms.Path(process.phfCoincFilter2Th5)
# process.pphfCoincFilter3Th5 = cms.Path(process.phfCoincFilter3Th5)
# process.pphfCoincFilter4Th5 = cms.Path(process.phfCoincFilter4Th5)
# process.pphfCoincFilter5Th5 = cms.Path(process.phfCoincFilter5Th5)
# process.pphfCoincFilter1Th6 = cms.Path(process.phfCoincFilter1Th6)
# process.pphfCoincFilter2Th6 = cms.Path(process.phfCoincFilter2Th6)
# process.pphfCoincFilter3Th6 = cms.Path(process.phfCoincFilter3Th6)
# process.pphfCoincFilter4Th6 = cms.Path(process.phfCoincFilter4Th6)
# process.pphfCoincFilter5Th6 = cms.Path(process.phfCoincFilter5Th6)

# process.NoScraping = cms.EDFilter("FilterOutScraping",
# applyfilter = cms.untracked.bool(True),
# debugOn = cms.untracked.bool(False),
# numtrack = cms.untracked.uint32(10),
# thresh = cms.untracked.double(0.25)
# )
# process.pBeamScrapingFilter=cms.Path(process.NoScraping)

# process.pAna = cms.EndPath(process.skimanalysis)

# from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
# process.hltfilter = hltHighLevel.clone(
#    HLTPaths = [
#        "HLT_HIUPC_SingleJet8_ZDC1nXOR_MaxPixelCluster50000_v2",                                                     
#        "HLT_HIMinimumBiasHF1AND_v3",
#    ]
# )
# process.filterSequence = cms.Sequence(
#    process.hltfilter
# )

# process.superFilterPath = cms.Path(process.filterSequence)
# process.skimanalysis.superFilters = cms.vstring("superFilterPath")

# for path in process.paths:
#    getattr(process, path)._seq = process.filterSequence * getattr(process,path)._seq


# Select the types of jets filled
addR2Jets = False
addR3Jets = False
addR3FlowJets = False
addR4Jets = True
addR4FlowJets = False
addUnsubtractedR4Jets = False

# Choose which additional information is added to jet trees
doHIJetID = True             # Fill jet ID and composition information branches
doWTARecluster = False        # Add jet phi and eta for WTA axis

# this is only for non-reclustered jets
addCandidateTagging = False

# Add addional jet collections
process.load("HeavyIonsAnalysis.JetAnalysis.extraJets_cff")
from HeavyIonsAnalysis.JetAnalysis.clusterJetsFromMiniAOD_cff import setupPprefJets

if addR2Jets :
    process.jetsR2 = cms.Sequence()
    setupPprefJets('ak2PF', process.jetsR2, process, isMC = 0, radius = 0.20, JECTag = '\
AK2PF')
    process.ak2PFpatJetCorrFactors.levels = ['L2Relative', 'L3Absolute']
    process.ak2PFpatJetCorrFactors.primaryVertices = "offlineSlimmedPrimaryVertices"
    process.load("HeavyIonsAnalysis.JetAnalysis.candidateBtaggingMiniAOD_cff")
    process.ak2PFJetAnalyzer.jetTag = 'ak2PFpatJets'
    process.ak2PFJetAnalyzer.jetName = 'ak2PF'
    process.ak2PFJetAnalyzer.doSubEvent = False # Need to disable this, since there is some issue with the gen jet constituents. More debugging needed is want to use constituents. 
    process.forest += process.extraJets * process.jetsR2 * process.ak2PFJetAnalyzer

if addR3Jets :
    process.jetsR3 = cms.Sequence()
    setupPprefJets('ak3PF', process.jetsR3, process, isMC = 0, radius = 0.30, JECTag = '\
AK3PF')
    process.ak3PFpatJetCorrFactors.levels = ['L2Relative', 'L3Absolute']
    process.ak3PFpatJetCorrFactors.primaryVertices = "offlineSlimmedPrimaryVertices"
    process.load("HeavyIonsAnalysis.JetAnalysis.candidateBtaggingMiniAOD_cff")
    process.ak3PFJetAnalyzer.jetTag = 'ak3PFpatJets'
    process.ak3PFJetAnalyzer.jetName = 'ak3PF'
    process.ak3PFJetAnalyzer.doSubEvent = False # Need to disable this, since there is some issue with the gen jet constituents. More debugging needed is want to use constituents. 
    process.forest += process.extraJets * process.jetsR3 * process.ak3PFJetAnalyzer

if addR4Jets :
    # Recluster using an alias "0" in order not to get mixed up with the default AK4 collections                                                                                    
    print("ADD R4 JETS")
    process.jetsR4 = cms.Sequence()
    setupPprefJets('ak04PF', process.jetsR4, process, isMC = 0, radius = 0.40, JECTag = 'AK4PF')
    process.ak04PFpatJetCorrFactors.levels = ['L2Relative', 'L3Absolute']
    process.ak04PFpatJetCorrFactors.primaryVertices = "offlineSlimmedPrimaryVertices"
    process.load("HeavyIonsAnalysis.JetAnalysis.candidateBtaggingMiniAOD_cff")
    process.ak4PFJetAnalyzer.jetTag = 'ak04PFpatJets'
    process.ak4PFJetAnalyzer.jetName = 'ak04PF'
    process.ak4PFJetAnalyzer.jetPtMin = cms.double(5.0)
    process.ak4PFJetAnalyzer.doSubEvent = False # Need to disable this, since there is some issue with the gen jet constituents. More debugging needed is want to use constituents. 
    process.forest += process.extraJets * process.jetsR4 * process.ak4PFJetAnalyzer
