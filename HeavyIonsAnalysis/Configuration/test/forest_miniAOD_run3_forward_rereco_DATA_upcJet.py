### HiForest Configuration
# Input: miniAOD
# Type: data

import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run3_2023_UPC_cff import Run3_2023_UPC
process = cms.Process('HiForest',Run3_2023_UPC)

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
        'root://xrootd-cms.infn.it//store/hidata/HIRun2023A/HIForward13/MINIAOD/16Jan2024-v1/2810000/b54df9bf-425b-4546-a07c-49eb475d5089.root'
        # 'file:/u/user/lidakal/SE_UserHome/input_files/2023/forward/b54df9bf-425b-4546-a07c-49eb475d5089.root'
        # '/store/hidata/HIRun2023A/HIForward2/MINIAOD/PromptReco-v2/000/375/531/00000/077422bb-4945-4e2c-9254-d07d5f7f8816.root'
        # HIEmptyBX
        # '/store/hidata/HIRun2023A/HIEmptyBX/MINIAOD/16Jan2024-v1/40000/c2e833a3-e8f0-451e-aa4c-e93557902996.root',
        # '/store/hidata/HIRun2023A/HIEmptyBX/MINIAOD/16Jan2024-v1/40000/782cdbae-6d30-457f-b1a1-d277e4ea89f9.root',
        # '/store/hidata/HIRun2023A/HIEmptyBX/MINIAOD/16Jan2024-v1/40000/6c4021f6-009d-4589-9ec1-a8aec84cad4b.root',
        # '/store/hidata/HIRun2023A/HIEmptyBX/MINIAOD/16Jan2024-v1/40000/26855567-b19a-4fc7-a9d8-d49150c44c77.root',
        # '/store/hidata/HIRun2023A/HIEmptyBX/MINIAOD/16Jan2024-v1/40000/5e1b3d75-a62d-4202-92b5-19ede9070881.root',
        # '/store/hidata/HIRun2023A/HIEmptyBX/MINIAOD/16Jan2024-v1/40000/93c45921-dbf5-446e-b7d0-8c347c6d9be1.root',
        # '/store/hidata/HIRun2023A/HIEmptyBX/MINIAOD/16Jan2024-v1/40000/9961be13-b6ef-4833-91af-7b4f7322e52b.root',
        # '/store/hidata/HIRun2023A/HIEmptyBX/MINIAOD/16Jan2024-v1/40000/a0411f51-3044-4ba8-ae5e-546d08626a22.root',
        # '/store/hidata/HIRun2023A/HIEmptyBX/MINIAOD/16Jan2024-v1/40000/7bf5251d-20cf-4a19-87a6-560c736e320d.root',
        # '/store/hidata/HIRun2023A/HIEmptyBX/MINIAOD/16Jan2024-v1/40000/bc34120f-d95b-48a1-928b-a14c5f30e069.root',
        # '/store/hidata/HIRun2023A/HIEmptyBX/MINIAOD/16Jan2024-v1/2810000/3f1dde54-ce9f-4201-8093-bb76a63c689a.root',
        # '/store/hidata/HIRun2023A/HIEmptyBX/MINIAOD/16Jan2024-v1/2810000/42e2049a-ea8b-48a8-ac6b-ac1943c1c09c.root',
        # '/store/hidata/HIRun2023A/HIEmptyBX/MINIAOD/16Jan2024-v1/2810000/e3eb78ef-2cae-4c8e-ae8e-076b8c0414fb.root'
    ),
)

# lumi list
import FWCore.PythonUtilities.LumiList as LumiList
# process.source.lumisToProcess = LumiList.LumiList(filename = '/u/user/lidakal/UPC2023/CMSSW_13_2_15/src/HeavyIonsAnalysis/Configuration/test/Cert_Collisions2023HI_374288_375823_Good_ZDC_Golden.json').getVLuminosityBlockRange()

# number of events to process, set to -1 to process all events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
)

###############################################################################

# load Global Tag, geometry, etc.
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '132X_dataRun3_Prompt_v7', '')
process.HiForestInfo.GlobalTagLabel = process.GlobalTag.globaltag

###############################################################################

# root output
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("HiForestMiniAOD.root"))

# edm output for debugging purposes
# process.output = cms.OutputModule(
#     "PoolOutputModule",
#     fileName = cms.untracked.string('HiForestEDM.root'),
#     outputCommands = cms.untracked.vstring(
#         'keep *',
#         "drop *_*ParticleTransformerAK4*_*_*",
#         )
#     )
# process.output_path = cms.EndPath(process.output)

###############################################################################

# event analysis
process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_data_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.skimanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hltobject_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.l1object_cfi')
process.hiEvtAnalyzer.doCentrality = cms.bool(False)
process.hiEvtAnalyzer.doHFfilters = cms.bool(False)

from HeavyIonsAnalysis.EventAnalysis.hltobject_cfi import trigger_list_data_2023
process.hltobject.triggerNames = trigger_list_data_2023

process.load('HeavyIonsAnalysis.EventAnalysis.particleFlowAnalyser_cfi')
process.particleFlowAnalyser.ptMin = cms.double(0.0)
process.particleFlowAnalyser.absEtaMax = cms.double(5.2)

################################
# electrons, photons, muons
process.load('HeavyIonsAnalysis.EGMAnalysis.ggHiNtuplizer_cfi')
process.ggHiNtuplizer.doMuons = cms.bool(False) # unpackedMuons collection not found from file
process.ggHiNtuplizer.useValMapIso = cms.bool(False) # True here causes seg fault
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
################################
# jet reco sequence - in forward/UPC we don't want subtracted jets
process.load("HeavyIonsAnalysis.JetAnalysis.ak4PFJetSequence_ppref_data_cff")

################################
# tracks
process.load("HeavyIonsAnalysis.TrackAnalysis.TrackAnalyzers_cff")
# muons
process.load("HeavyIonsAnalysis.MuonAnalysis.muonAnalyzer_cfi")
###############################################################################

# ZDC RecHit producer + analyzer
process.load('HeavyIonsAnalysis.ZDCAnalysis.ZDCAnalyzersHC2023_cff')

###############################################################################
# main forest sequence
process.forest = cms.Path(
    process.HiForestInfo +
    process.hiEvtAnalyzer +
    process.hltanalysis +
    # process.hltobject +
    # process.l1object +
    # process.trackSequencePP +
    process.particleFlowAnalyser +
    # process.ggHiNtuplizer +
    # process.muonSequencePP +
    process.zdcSequence 
    )

#customisation

#####################################################################################
# Select the types of jets filled
matchJets = True             # Enables q/g and heavy flavor jet identification in MC
jetPtMin = 5
jetAbsEtaMax = 2.5

# Choose which additional information is added to jet trees
doHIJetID = True             # Fill jet ID and composition information branches
doWTARecluster = False        # Add jet phi and eta for WTA axis
doBtagging  =  False         # Note that setting to True increases computing time a lot

# 0 means use original mini-AOD jets, otherwise use R value, e.g., 3,4,8
jetLabel = "0"

# add candidate tagging, copy/paste to add other jet radii
from HeavyIonsAnalysis.JetAnalysis.setupJets_ppRef_cff import candidateBtaggingMiniAOD
candidateBtaggingMiniAOD(process, isMC = False, jetPtMin = jetPtMin, jetCorrLevels = ['L2Relative', 'L3Absolute'], doBtagging = doBtagging, labelR = jetLabel)
# setup jet analyzer

setattr(process,"ak"+jetLabel+"PFJetAnalyzer",process.ak4PFJetAnalyzer.clone())
getattr(process,"ak"+jetLabel+"PFJetAnalyzer").jetTag = 'slimmedJets'
getattr(process,"ak"+jetLabel+"PFJetAnalyzer").jetName = 'ak'+jetLabel+'PF'
getattr(process,"ak"+jetLabel+"PFJetAnalyzer").matchJets = matchJets
getattr(process,"ak"+jetLabel+"PFJetAnalyzer").matchTag = 'patJetsAK'+jetLabel+'PFUnsubJets'
getattr(process,"ak"+jetLabel+"PFJetAnalyzer").doHiJetID = doHIJetID
getattr(process,"ak"+jetLabel+"PFJetAnalyzer").doWTARecluster = doWTARecluster
getattr(process,"ak"+jetLabel+"PFJetAnalyzer").jetPtMin = jetPtMin
getattr(process,"ak"+jetLabel+"PFJetAnalyzer").jetAbsEtaMax = cms.untracked.double(jetAbsEtaMax)
getattr(process,"ak"+jetLabel+"PFJetAnalyzer").rParam = int(jetLabel)*0.1
if doBtagging:
    getattr(process,"ak"+jetLabel+"PFJetAnalyzer").pfJetProbabilityBJetTag = cms.untracked.string("pfJetProbabilityBJetTagsDeepFlavour")
    getattr(process,"ak"+jetLabel+"PFJetAnalyzer").pfUnifiedParticleTransformerAK4JetTags = cms.untracked.string("pfUnifiedParticleTransformerAK4JetTagsDeepFlavour")
process.forest += getattr(process,"ak"+jetLabel+"PFJetAnalyzer")

#########################
# Event Selection -> add the needed filters here
#########################
process.load('HeavyIonsAnalysis.EventAnalysis.collisionEventSelection_cff')
process.pclusterCompatibilityFilter = cms.Path(process.clusterCompatibilityFilter)
process.pprimaryVertexFilter = cms.Path(process.primaryVertexFilter)
process.load('HeavyIonsAnalysis.ZDCAnalysis.HiZDCfilter_cfi')

process.pAna = cms.EndPath(process.skimanalysis)

# process.NoScraping = cms.EDFilter("FilterOutScraping",
# applyfilter = cms.untracked.bool(True),
# debugOn = cms.untracked.bool(False),
# numtrack = cms.untracked.uint32(10),
# thresh = cms.untracked.double(0.25)
# )
# process.pBeamScrapingFilter=cms.Path(process.NoScraping)

from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltfilter = hltHighLevel.clone(
   HLTPaths = [
    #    "HLT_HIUPC_*",
       "HLT_HIUPC_SingleJet*_NotMBHF2AND_MaxPixelCluster50000_v*",
   ]
)
# process.hltfilter.andOr = cms.bool(True)  # True = OR, False = AND between the HLT paths
# process.hltfilter.throw = cms.bool(False) # throw exception on unknown path names

process.filterSequence = cms.Sequence(
    process.hltfilter
    * process.primaryVertexFilter 
    * process.clusterCompatibilityFilter
    # * (process.zdcreco2023HardCode + process.zdcEnergyFilter0nAnd)
)
process.prefilter = cms.Path(process.filterSequence)
process.skimanalysis.superFilters = cms.vstring("prefilter")
for path in process.paths:
      getattr(process, path)._seq = process.filterSequence * getattr(process,path)._seq

## Sometimes needed
# process.options = cms.untracked.PSet(
#     SkipEvent = cms.untracked.vstring('ProductNotFound')
# )

### Message Loggers

process.MessageLogger.cerr.FwkReport.reportEvery = 10000

# # (debug) See all paths in the process
# for name, path in process.paths_().items():
#     print(name, ":", path)

# # Find where the offending module came from
# if hasattr(process, 'patJetGenJetMatchAK0PFUnsubJets'):
#     print(process.patJetGenJetMatchAK0PFUnsubJets.dumpPython())

# if hasattr(process, 'unsubUpdatedPatJetsAK4PFCHS'):
#     print(process.unsubUpdatedPatJetsAK4PFCHS.dumpPython())

print(process.forest.dumpPython())