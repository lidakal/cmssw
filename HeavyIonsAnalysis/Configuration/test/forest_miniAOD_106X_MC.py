### HiForest Configuration
# Input: miniAOD
# Type: mc

import FWCore.ParameterSet.Config as cms
process = cms.Process('HiForest')

###############################################################################

# HiForest info
process.load("HeavyIonsAnalysis.EventAnalysis.HiForestInfo_cfi")
process.HiForestInfo.info = cms.vstring("HiForest, miniAOD, 106X, mc")

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
        '/store/mc/RunIISummer19UL17MiniAOD/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/100000/BFAAC85A-F5C5-8843-8D2A-76A9E873E24B.root'
        ),
    )

# number of events to process, set to -1 to process all events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
    )

# To skip events with 'Found zero products matching all criteria' error
# process.options = cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound'))

###############################################################################

# load Global Tag, geometry, etc.
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')


from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mc2017_realistic_v8', '')
process.HiForestInfo.GlobalTagLabel = process.GlobalTag.globaltag
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")


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

#############################
# Gen Analyzer
#############################
process.load('HeavyIonsAnalysis.EventAnalysis.HiGenAnalyzer_cfi')
# making cuts looser so that we can actually check dNdEta
process.HiGenParticleAna.ptMin = cms.untracked.double(0.4) # default is 5
process.HiGenParticleAna.etaMax = cms.untracked.double(5.) # default is 2.5

# event analysis
process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.particleFlowAnalyser_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_mc_cfi')
#process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.skimanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hltobject_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.l1object_cfi')

from HeavyIonsAnalysis.EventAnalysis.hltobject_cfi import trigger_list_mc
process.hltobject.triggerNames = trigger_list_mc

################################
# electrons, photons, muons
SS2018PbPbMC = "HeavyIonsAnalysis/EGMAnalysis/data/SS2018PbPbMC.dat"
process.load('HeavyIonsAnalysis.EGMAnalysis.correctedElectronProducer_cfi')
process.correctedElectrons.correctionFile = SS2018PbPbMC

process.load('HeavyIonsAnalysis.EGMAnalysis.ggHiNtuplizer_cfi')
process.ggHiNtuplizer.doGenParticles = cms.bool(True)
process.ggHiNtuplizer.doMuons = cms.bool(False)
process.ggHiNtuplizer.electronSrc = "correctedElectrons"
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
################################
# jet reco sequence
process.load('HeavyIonsAnalysis.JetAnalysis.ak4PFJetSequence_pponPbPb_mc_cff')
################################
# tracks
process.load("HeavyIonsAnalysis.TrackAnalysis.TrackAnalyzers_cff")
#muons
process.load("HeavyIonsAnalysis.MuonAnalysis.unpackedMuons_cfi")
process.load("HeavyIonsAnalysis.MuonAnalysis.muonAnalyzer_cfi")
process.muonAnalyzer.doGen = cms.bool(True)

###############################################################################



###############################################################################
# main forest sequence
process.forest = cms.Path(
    process.HiForestInfo +
    # process.hltanalysis +
    #process.hltobject +
    #process.l1object +
    process.trackSequencePbPb +
    process.particleFlowAnalyser +
    process.hiEvtAnalyzer +
    process.HiGenParticleAna +
    #process.unpackedMuons +
    #process.correctedElectrons #+
    #process.ggHiNtuplizer +
    #process.muonAnalyzer + 
    process.ak4PFJetAnalyzer # cms.EDAnalyzer("HiInclusiveJetAnalyzer")
    )

#customisation

addCandidateTagging = True

if addCandidateTagging:
    process.load("HeavyIonsAnalysis.JetAnalysis.candidateBtaggingMiniAOD_cff")
    bTagPrefix_ = 'TEST'

    from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
    updateJetCollection(
        process,
        jetSource = cms.InputTag('slimmedJets'),
        jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
        btagDiscriminators = ['pfCombinedSecondaryVertexV2BJetTags', 
                              'pfDeepCSVJetTags:probb',
                              'pfDeepCSVJetTags:probc',
                              'pfDeepCSVJetTags:probudsg',
                              'pfDeepCSVJetTags:probbb',
                              'pfDeepFlavourJetTags:probb',
                              'pfDeepFlavourJetTags:probbb',
                              'pfDeepFlavourJetTags:problepb',
                              'pfDeepFlavourJetTags:probc',
                              'pfDeepFlavourJetTags:probuds',
                              'pfDeepFlavourJetTags:probg'
                               ], ## to add discriminators,
        btagPrefix = bTagPrefix_,
    )

    process.updatedPatJets.addJetCorrFactors = False

    ## NOT NEEDED if the btagDiscriminators are defined correctly in updateJetCollection
    # process.updatedPatJets.addBTagInfo = True
    # process.updatedPatJets.addDiscriminators = True
    # process.updatedPatJets.discriminatorSources = cms.VInputTag(
        # cms.InputTag('pfDeepCSVJetTags:probb'),
        # cms.InputTag('pfDeepCSVJetTags:probc'),
        # cms.InputTag('pfDeepCSVJetTags:probudsg'),
        # cms.InputTag('pfDeepCSVJetTags:probbb'),
        # cms.InputTag('pfDeepFlavourJetTags:probb')
    # )

    ## NEEDED for track aggregation
    process.updatedPatJets.addTagInfos = True
    process.updatedPatJets.tagInfoSources = cms.VInputTag(
        cms.InputTag('pfImpactParameterTagInfos'),
        cms.InputTag('pfSecondaryVertexTagInfos'),
    )

    # process.pfSecondaryVertexTagInfos.trackSelection.jetDeltaRMax = 0.4

    process.forest.insert(1,process.candidateBtagging*process.updatedPatJets)
    # Note: candidateBtagging includes pfImpactParameterTagInfos

    process.ak4PFJetAnalyzer.jetTag = "updatedPatJets"
    process.ak4PFJetAnalyzer.bTagJetName = cms.untracked.string("")
    process.ak4PFJetAnalyzer.doTracks = cms.untracked.bool(True)
    process.ak4PFJetAnalyzer.doLegacyBtagging = cms.untracked.bool(False)
    process.ak4PFJetAnalyzer.doSvtx = cms.untracked.bool(True)

#########################
# Event Selection -> add the needed filters here
#########################

#process.load('HeavyIonsAnalysis.EventAnalysis.collisionEventSelection_cff')
#process.pclusterCompatibilityFilter = cms.Path(process.clusterCompatibilityFilter)
#process.pprimaryVertexFilter = cms.Path(process.primaryVertexFilter)
#process.pAna = cms.EndPath(process.skimanalysis)
