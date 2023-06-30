### HiForest Configuration
# Input: miniAOD
# Type: mc

import FWCore.ParameterSet.Config as cms
process = cms.Process('HiForest')
process.options = cms.untracked.PSet()

###############################################################################

# HiForest info
process.load("HeavyIonsAnalysis.EventAnalysis.HiForestInfo_cfi")
process.HiForestInfo.info = cms.vstring("HiForest, miniAOD, 106X, data")

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
        '/store/data/Run2017G/LowEGJet/MINIAOD/UL2017_MiniAODv2-v2/2810000/01869167-7867-434D-A952-5BEC77B73ABA.root'
        ),
    )


# Select specific event
# process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange('1:7897')
# process.source.eventsToProcess = cms.untracked.VEventRange('1:110394257')

# number of events to process, set to -1 to process all events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

# Multi-thread 
# process.options.numberOfThreads = cms.untracked.int32(8)

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
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v35', '')
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

# event analysis
process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.particleFlowAnalyser_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_data_cfi')
process.hiEvtAnalyzer.Vertex = cms.InputTag("offlineSlimmedPrimaryVertices")
process.hiEvtAnalyzer.doCentrality = cms.bool(False)
process.hiEvtAnalyzer.doEvtPlane = cms.bool(False)
process.load('HeavyIonsAnalysis.EventAnalysis.skimanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hltobject_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.l1object_cfi')

# from HeavyIonsAnalysis.EventAnalysis.hltobject_cfi import trigger_list_mc
# process.hltobject.triggerNames = trigger_list_mc
# process.hltobject.triggerNames = cms.vstring('HLT_HIAK4CaloJet30_v', 'HLT_HIAK4CaloJet40_v', 'HLT_HIAK4CaloJet60_v')

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
process.load('HeavyIonsAnalysis.JetAnalysis.akCs4PFJetSequence_pponPbPb_data_cff')
process.tagInfoSequence = cms.Sequence()
process.recoJetSequence = cms.Sequence()
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
    process.hltanalysis +
    # process.hltobject +
    #process.l1object +
    # process.trackSequencePbPb +
    # process.particleFlowAnalyser +
    process.hiEvtAnalyzer +
    # process.updatePATJetSequence + 
    process.tagInfoSequence +
    process.recoJetSequence + 
    #process.unpackedMuons +
    #process.correctedElectrons #+
    #process.ggHiNtuplizer +
    #process.muonAnalyzer + 
    process.akCs4PFJetAnalyzer # cms.EDAnalyzer("HiInclusiveJetAnalyzer")
    )

## Customization

ipTagInfoLabel_ = "pfImpactParameter"
svTagInfoLabel_ = "pfInclusiveSecondaryVertexFinder"
addTagInfos = True
if addTagInfos:
    ## Impact parameter tag infos
    process.load("RecoBTag.ImpactParameter.pfImpactParameterTagInfos_cfi")
    process.pfImpactParameterTagInfos.candidates  = "packedPFCandidates"
    process.pfImpactParameterTagInfos.primaryVertex = "offlineSlimmedPrimaryVertices"
    process.pfImpactParameterTagInfos.jets = "slimmedJets"

    ## Secondary vertex tag infos
    process.load("RecoBTag.SecondaryVertex.pfInclusiveSecondaryVertexFinderTagInfos_cfi")
    process.pfInclusiveSecondaryVertexFinderTagInfos.extSVCollection = "slimmedSecondaryVertices"

    _btagDiscriminators = cms.PSet(
        names = cms.vstring(
            'pfDeepFlavourJetTags:probb',
            'pfDeepFlavourJetTags:probbb',
            'pfDeepFlavourJetTags:problepb'
        )
    )
    from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
    updateJetCollection(
        process,
        jetSource = cms.InputTag('slimmedJets'),
        jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
        )
    process.updatedPatJets.addJetCorrFactors = False
    process.updatedPatJets.addTagInfos = True
    process.updatedPatJets.addBTagInfo = cms.bool(True) ## Needed to add tag infos
    process.updatedPatJets.tagInfoSources = cms.VInputTag(
        cms.InputTag(ipTagInfoLabel_ + "TagInfos"),
        cms.InputTag(svTagInfoLabel_ + "TagInfos"),
    )
    process.tagInfoSequence.insert(0, process.pfImpactParameterTagInfos *
                                   process.pfInclusiveSecondaryVertexFinderTagInfos *
                                   process.updatedPatJets)
    process.akCs4PFJetAnalyzer.jetTag = "updatedPatJets"

process.akCs4PFJetAnalyzer.doSubJetsNew = cms.untracked.bool(True)

doTracks = True
if doTracks:
    process.akCs4PFJetAnalyzer.doTracks = cms.untracked.bool(True)
    process.akCs4PFJetAnalyzer.ipTagInfoLabel = cms.untracked.string(ipTagInfoLabel_)
    # process.akCs4PFJetAnalyzer.trkPtCut = cms.untracked.double(1.)

doSvtx = True
if doSvtx:
    process.akCs4PFJetAnalyzer.doSvtx = cms.untracked.bool(True)
    process.akCs4PFJetAnalyzer.svTagInfoLabel = cms.untracked.string(svTagInfoLabel_)

doDeclustering = True
doAggregation = False
doChargedOnly = True
doLatekt_ = False

tmva_variables = ["trkIp3dSig", "trkIp2dSig", "trkDistToAxis",
                  "svtxdls", "svtxdls2d", "svtxm", "svtxmcorr",
                  "svtxnormchi2", "svtxNtrk", "svtxTrkPtOverSv",
                  "jtpt"]
if doDeclustering:
    process.load("RecoHI.HiJetAlgos.dynGroomedPATJets_cfi")
    process.dynGroomedPFJets = process.dynGroomedPATJets.clone(
        isMC = cms.bool(False),
        chargedOnly = cms.bool(doChargedOnly),
        aggregateHF = cms.bool(doAggregation),
        # aggregateHF = cms.bool(False),
        jetSrc = cms.InputTag("updatedPatJets"),
        constitSrc = cms.InputTag("packedPFCandidates"),
        doGenJets = cms.bool(False),
        # candToGenParticleMap = cms.InputTag("TrackToGenParticleMapProducer", "trackToGenParticleMap"),
        aggregateWithTruthInfo = cms.bool(False),
        aggregateWithXGB = cms.bool(False),
        aggregateWithTMVA = cms.bool(True),
        aggregateWithCuts = cms.bool(False),
        xgb_path = cms.FileInPath("RecoHI/HiJetAlgos/data/sig_vs_bkg.model"),
        tmva_path = cms.FileInPath("RecoHI/HiJetAlgos/data/TMVAClassification_BDTG.weights.xml"),
        tmva_variables = cms.vstring(tmva_variables),
        doLateKt = doLatekt_
    )
    process.recoJetSequence += process.dynGroomedPFJets
    process.akCs4PFJetAnalyzer.groomedJets = cms.untracked.InputTag("dynGroomedPFJets")
    ## creates the reco subjets
    


#########################
# Event Selection -> add the needed filters here
#########################

#process.load('HeavyIonsAnalysis.EventAnalysis.collisionEventSelection_cff')
#process.pclusterCompatibilityFilter = cms.Path(process.clusterCompatibilityFilter)
#process.pprimaryVertexFilter = cms.Path(process.primaryVertexFilter)
#process.pAna = cms.EndPath(process.skimanalysis)
