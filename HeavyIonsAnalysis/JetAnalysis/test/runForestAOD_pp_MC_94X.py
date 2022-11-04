### HiForest Configuration
# Collisions: pp
# Type: MC
# Input: AOD

import FWCore.ParameterSet.Config as cms
process = cms.Process('HiForest')
process.options = cms.untracked.PSet()

#####################################################################################
# HiForest labelling info
#####################################################################################

process.load("HeavyIonsAnalysis.JetAnalysis.HiForest_cff")
process.HiForest.inputLines = cms.vstring("HiForest 94X",)
import subprocess, os
version = subprocess.check_output(['git',
    '-C', os.path.expandvars('$CMSSW_BASE/src'), 'describe', '--tags'])
if version == '':
    version = 'no git info'
process.HiForest.HiForestVersion = cms.string(version)

#####################################################################################
# Input source
#####################################################################################

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        "/store/himc/RunIIpp5Spring18DR/QCD_pThat-15_bJet_TuneCP5_5p02TeV_pythia8/AODSIM/94X_mc2017_realistic_forppRef5TeV_v1-v1/230000/92D50574-667C-E911-B0CA-1866DAEECF18.root"
        #'file:/data_CMS/cms/mnguyen/00D396CF-878A-E911-88E8-14187741278B.root'
    ),
                            #skipEvents = cms.untracked.uint32(1242)
)
# Select a specific event (propagated to CRAB => REMOVE AFTERWARDS!)
# process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange('1:12227')
# process.source.eventsToProcess = cms.untracked.VEventRange('1:148128582')

# Number of events we want to process, -1 = all events
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

#####################################################################################
# Load Global Tag, Geometry, etc.
#####################################################################################

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_forppRef5TeV', '')
process.HiForest.GlobalTagLabel = process.GlobalTag.globaltag

process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
             tag = cms.string("JPcalib_MC94X_2017pp_v2"),
             connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS")

         )
      ])

#####################################################################################
# Define tree output
#####################################################################################

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("HiForestAOD.root"))

#####################################################################################
# Additional Reconstruction and Analysis: Main Body
#####################################################################################

#############################
# Jets
#############################
process.load("HeavyIonsAnalysis.JetAnalysis.fullJetSequence_pp_mc_cff")

# Use this version for JEC
# process.load("HeavyIonsAnalysis.JetAnalysis.fullJetSequence_pp_jec_cff")

from HeavyIonsAnalysis.Configuration.CommonFunctions_cff import overrideJEC_MC_pp5020_2017  
process = overrideJEC_MC_pp5020_2017(process)
#####################################################################################

############################
# Event Analysis
############################
# use data version to avoid PbPb MC
process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_data_cfi')
process.hiEvtAnalyzer.Vertex = cms.InputTag("offlinePrimaryVertices")
process.hiEvtAnalyzer.doCentrality = cms.bool(False)
process.hiEvtAnalyzer.doEvtPlane = cms.bool(False)
process.hiEvtAnalyzer.doMC = cms.bool(True) # general MC info
process.hiEvtAnalyzer.doHiMC = cms.bool(False) # HI specific MC info

process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cff')
process.load('HeavyIonsAnalysis.EventAnalysis.hltobject_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.l1object_cfi')

process.load('HeavyIonsAnalysis.JetAnalysis.HiGenAnalyzer_cfi')
#process.HiGenParticleAna.genParticleSrc = cms.untracked.InputTag("genPartonsForJets")
process.HiGenParticleAna.genParticleRVSrc = cms.untracked.InputTag("genPartonsForJets")
process.HiGenParticleAna.doHI = False
process.HiGenParticleAna.etaMax = cms.untracked.double(2.5) # default is 2 
process.HiGenParticleAna.stableOnly = False
process.HiGenParticleAna.useRefVector = cms.untracked.bool(True)

process.bHadronAna = process.HiGenParticleAna.clone(
    genParticleRVSrc = cms.untracked.InputTag("selectedHadronsAndPartons","bHadrons"),
    useRefVector = True
)
process.cHadronAna = process.bHadronAna.clone(
    genParticleRVSrc = cms.untracked.InputTag("selectedHadronsAndPartons","cHadrons")
)
process.leptonAna = process.HiGenParticleAna.clone(
    genParticleRVSrc = cms.untracked.InputTag("selectedHadronsAndPartons","leptons")
)
process.outgoingPartonAna = process.bHadronAna.clone(
    genParticleRVSrc = cms.untracked.InputTag("selectedHadronsAndPartons","algorithmicPartons")
)
process.incomingPartonAna = process.bHadronAna.clone(
    genParticleSrc = cms.untracked.InputTag("genParticles"),
    useRefVector = cms.untracked.bool(False),
    partonMEOnly = cms.untracked.bool(True)
)

process.bDecayAna = process.incomingPartonAna.clone(
    #genParticleSrc = cms.untracked.InputTag("HFdecayProductTagger"),
    genParticleSrc = cms.untracked.InputTag("HFdecayProductTagger"),
    partonMEOnly = False,
    chargedOnly = True
)

process.load('HeavyIonsAnalysis.EventAnalysis.runanalyzer_cff')

process.load("HeavyIonsAnalysis.JetAnalysis.pfcandAnalyzer_pp_cfi")
process.pfcandAnalyzer.skipCharged      = False
process.pfcandAnalyzer.pfPtMin          = 0
process.pfcandAnalyzer.pfCandidateLabel = cms.InputTag("particleFlow")
process.pfcandAnalyzer.doVS             = cms.untracked.bool(False)
process.pfcandAnalyzer.doUEraw_         = cms.untracked.bool(False)
process.pfcandAnalyzer.genLabel         = cms.InputTag("genParticles")

#####################################################################################

#########################
# Track Analyzer
#########################
process.load('HeavyIonsAnalysis.JetAnalysis.ExtraTrackReco_cff')
process.load('HeavyIonsAnalysis.JetAnalysis.TrkAnalyzers_cff')

# Use this instead for track corrections
# process.load('HeavyIonsAnalysis.JetAnalysis.TrkAnalyzers_Corr_cff')

#####################################################################################

#####################
# photons
######################
process.load('HeavyIonsAnalysis.PhotonAnalysis.ggHiNtuplizer_cfi')

####################################################################################

#####################
# Electron ID
#####################

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format to be processed
# DataFormat.AOD or DataFormat.MiniAOD
dataFormat = DataFormat.AOD
switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce. https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2#Recipe_for_regular_users_for_7_4
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff']

# add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

#####################################################################################



#########################
# Main analysis list
#########################

process.ana_step = cms.Path(
    process.hltanalysis *
    process.hiEvtAnalyzer *
    #process.hltobject +
    # process.l1object +
    process.genJetSequence +
    process.jetSequence +
    process.bHadronAna +
    process.cHadronAna +
    #process.leptonAna +
    #process.outgoingPartonAna + 
    process.incomingPartonAna + 
    process.bDecayAna + 
    # Should be added in the path for VID module
    # process.egmGsfElectronIDSequence +
    #process.ggHiNtuplizer +
    #process.ggHiNtuplizerGED +
    #process.pfcandAnalyzer +
    process.HiForest +
    #process.trackSequencesPP +
    #process.HiGenParticleAna*
    process.runAnalyzer
)

#####################################################################################
'''
# edm output for debugging purposes
process.output = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('HiForestEDM.root'),
    outputCommands = cms.untracked.vstring(
        'keep *',
        # drop aliased products
        'drop *_akULPu3PFJets_*_*',
        'drop *_akULPu4PFJets_*_*',
        )
    )

process.output_path = cms.EndPath(process.output)
'''

#########################
# Event Selection
#########################

process.load('HeavyIonsAnalysis.JetAnalysis.EventSelection_cff')
process.pHBHENoiseFilterResultProducer = cms.Path(process.HBHENoiseFilterResultProducer)
process.HBHENoiseFilterResult = cms.Path(process.fHBHENoiseFilterResult)
process.HBHENoiseFilterResultRun1 = cms.Path(process.fHBHENoiseFilterResultRun1)
process.HBHENoiseFilterResultRun2Loose = cms.Path(process.fHBHENoiseFilterResultRun2Loose)
process.HBHENoiseFilterResultRun2Tight = cms.Path(process.fHBHENoiseFilterResultRun2Tight)
process.HBHEIsoNoiseFilterResult = cms.Path(process.fHBHEIsoNoiseFilterResult)

process.PAprimaryVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && abs(z) <= 25 && position.Rho <= 2 && tracksSize >= 2"),
    filter = cms.bool(True), # otherwise it won't filter the events
)

process.NoScraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)

process.pPAprimaryVertexFilter = cms.Path(process.PAprimaryVertexFilter)
process.pBeamScrapingFilter=cms.Path(process.NoScraping)

process.load("HeavyIonsAnalysis.VertexAnalysis.PAPileUpVertexFilter_cff")
process.pVertexFilterCutG = cms.Path(process.pileupVertexFilterCutG)
process.pVertexFilterCutGloose = cms.Path(process.pileupVertexFilterCutGloose)
process.pVertexFilterCutGtight = cms.Path(process.pileupVertexFilterCutGtight)
process.pVertexFilterCutGplus = cms.Path(process.pileupVertexFilterCutGplus)
process.pVertexFilterCutE = cms.Path(process.pileupVertexFilterCutE)
process.pVertexFilterCutEandG = cms.Path(process.pileupVertexFilterCutEandG)

process.pAna = cms.EndPath(process.skimanalysis)

# Customization
# write heavy stuff for debugging b-tagging
process.ak4PFJetAnalyzer.doLifeTimeTaggingExtra = True
process.ak4PFJetAnalyzer.doTrackMatching = True

#option to do late soft drop
process.dynGroomedGenJets.doLateSD = False
process.dynGroomedPFJets.doLateSD = False

## don't change these
process.genParticlesForJets.undecayHF = cms.bool(False)
process.genPartonsForJets.undecayHF = cms.bool(False)
process.genParticlesForJets.chargedOnly = cms.bool(False)  
process.genPartonsForJets.chargedOnly = cms.bool(False) 
## change these to run charged only or full jet declustering
# model_path = cms.string("$CMSSW_BASE/src/RecoHI/HiJetAlgos/data/trained_bst.model")
# model_path = cms.string("../../../RecoHI/HiJetAlgos/data/trained_bst.model")
# model_path = cms.string("/grid_mnt/vol_home/llr/cms/kalipoliti/CMSSW_9_4_10/src/RecoHI/HiJetAlgos/data/trained_bst.model")
model_path = cms.FileInPath("RecoHI/HiJetAlgos/data/trained_bst.model")
process.dynGroomedGenJets.chargedOnly = cms.bool(True)
process.dynGroomedGenJets.aggregateHF = cms.bool(True)
process.dynGroomedGenJets.model_path = model_path

process.dynGroomedPFJets.chargedOnly = cms.bool(True)
process.dynGroomedPFJets.aggregateHF = cms.bool(True)
process.dynGroomedPFJets.model_path = model_path


## replace b-hadron decays by the parent
process.load("RecoHI.HiJetAlgos.HFdecayProductTagger_cfi")
# option to tag B or C daughter in HFdecayProductTagger: True for B, False for C
# if aggregateHF == False, the particles will be tagged but not aggregated 
process.HFdecayProductTagger.tagBorC = cms.bool(True)

process.genJetSequence.insert(0,process.HFdecayProductTagger)
#process.genParticlesForJets.src = 'HFdecayProductTagger'

#process.load("RecoHI.HiJetAlgos.RecoHFHadronReplacer_cfi")
#process.jetSequence+=process.RecoHFHadronReplacer
