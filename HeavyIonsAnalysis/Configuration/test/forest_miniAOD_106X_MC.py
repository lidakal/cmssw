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
        # '/store/mc/RunIISummer19UL17MiniAOD/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/100000/BFAAC85A-F5C5-8843-8D2A-76A9E873E24B.root'
        # '/store/mc/RunIISummer20UL16MiniAODAPVv2/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_preVFP_v11-v1/120000/0230C4F8-6445-F74C-8409-27F77DFDE107.root'
        # '/store/mc/RunIILowPUSummer20UL17MiniAODv2/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/pilot_106X_mc2017_realistic_v9For2017H_v1-v2/2560000/D9746C19-3FD6-B246-95F4-E7DBD33ED768.root'
        '/store/mc/RunIISummer20UL17MiniAODv2/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_mc2017_realistic_v9-v1/240000/DC372BF6-08B2-9C4A-AF9E-69E80E066E4F.root'
        ),
    )

# number of events to process, set to -1 to process all events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
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
process.tagInfoSequence = cms.Sequence()
process.genJetSequence = cms.Sequence()
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
    # process.hltanalysis +
    #process.hltobject +
    #process.l1object +
    # process.trackSequencePbPb +
    process.particleFlowAnalyser +
    process.hiEvtAnalyzer +
    process.HiGenParticleAna +
    # process.updatePATJetSequence + 
    process.tagInfoSequence +
    process.genJetSequence + 
    process.recoJetSequence + 
    #process.unpackedMuons +
    #process.correctedElectrons #+
    #process.ggHiNtuplizer +
    #process.muonAnalyzer + 
    process.ak4PFJetAnalyzer # cms.EDAnalyzer("HiInclusiveJetAnalyzer")
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

    from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
    updateJetCollection(
        process,
        jetSource = cms.InputTag('slimmedJets'),
        jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None')
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
    process.ak4PFJetAnalyzer.jetTag = "updatedPatJets"

process.ak4PFJetAnalyzer.doSubJetsNew = cms.untracked.bool(True)

doTracks = True
if doTracks:
    process.ak4PFJetAnalyzer.doTracks = cms.untracked.bool(True)
    process.ak4PFJetAnalyzer.ipTagInfoLabel = cms.untracked.string(ipTagInfoLabel_)
    # process.ak4PFJetAnalyzer.trkPtCut = cms.untracked.double(1.)

doSvtx = True
if doSvtx:
    process.ak4PFJetAnalyzer.doSvtx = cms.untracked.bool(True)
    process.ak4PFJetAnalyzer.svTagInfoLabel = cms.untracked.string(svTagInfoLabel_)

doDeclustering = True
doAggregation = True
doChargedOnly = True
if doDeclustering:
    process.load("GeneratorInterface.RivetInterface.mergedGenParticles_cfi")
    process.genJetSequence += process.mergedGenParticles
    ## Produces a reco::GenParticleCollection named mergedGenParticles

    process.load("RecoHI.HiJetAlgos.HFdecayProductTagger_cfi")
    process.HFdecayProductTagger.genParticles = cms.InputTag("mergedGenParticles")
    process.HFdecayProductTagger.tagBorC = cms.bool(True) # tag B
    process.genJetSequence.insert(1, process.HFdecayProductTagger)
    taggedGenParticlesName_ = "HFdecayProductTagger"
    ## Produces a std::vector<pat::PackedGenParticle> named HFdecayProductTagger

    process.ak4PFJetAnalyzer.genParticles = cms.untracked.InputTag(taggedGenParticlesName_)

    process.bDecayAna = process.HiGenParticleAna.clone(
        genParticleSrc = cms.InputTag(taggedGenParticlesName_),
        useRefVector = cms.untracked.bool(False),
        partonMEOnly = cms.untracked.bool(False),
        chargedOnly = doChargedOnly,
        doHI = False,
        etaMax = cms.untracked.double(10),
        ptMin = cms.untracked.double(0),
        stableOnly = False
    )
    process.genJetSequence += process.bDecayAna
    ## Creates the gen particle ntuple bDecayAna/hi

    process.load("RecoHI.HiJetAlgos.TrackToGenParticleMapProducer_cfi")
    process.TrackToGenParticleMapProducer.jetSrc = cms.InputTag("updatedPatJets")
    process.TrackToGenParticleMapProducer.genParticleSrc = cms.InputTag(taggedGenParticlesName_)
    process.genJetSequence += process.TrackToGenParticleMapProducer
    ## Creates the genConstitToGenParticleMap and trackToGenParticleMap

    process.load("RecoHI.HiJetAlgos.dynGroomedPATJets_cfi")
    process.dynGroomedGenJets = process.dynGroomedPATJets.clone(
        chargedOnly = cms.bool(doChargedOnly),
        aggregateHF = cms.bool(doAggregation),
        jetSrc = cms.InputTag("updatedPatJets"),
        constitSrc = cms.InputTag("packedGenParticles"),
        doGenJets = cms.bool(True),
        candToGenParticleMap = cms.InputTag("TrackToGenParticleMapProducer", "genConstitToGenParticleMap")
    )
    process.genJetSequence += process.dynGroomedGenJets
    process.ak4PFJetAnalyzer.groomedGenJets = cms.untracked.InputTag("dynGroomedGenJets")
    ## Creates the gen jet subjets

    process.dynGroomedPFJets = process.dynGroomedPATJets.clone(
        chargedOnly = cms.bool(doChargedOnly),
        aggregateHF = cms.bool(doAggregation),
        jetSrc = cms.InputTag("updatedPatJets"),
        constitSrc = cms.InputTag("packedPFCandidates"),
        doGenJets = cms.bool(False),
        candToGenParticleMap = cms.InputTag("TrackToGenParticleMapProducer", "trackToGenParticleMap"),
        aggregateWithTruthInfo = cms.bool(True),
        aggregateWithBDT = cms.bool(False),
        model_path = cms.FileInPath("RecoHI/HiJetAlgos/data/trained_bst_30_pt_700.model")
    )
    process.recoJetSequence += process.dynGroomedPFJets
    process.ak4PFJetAnalyzer.groomedJets = cms.untracked.InputTag("dynGroomedPFJets")
    ## creates the reco subjets
    


#########################
# Event Selection -> add the needed filters here
#########################

#process.load('HeavyIonsAnalysis.EventAnalysis.collisionEventSelection_cff')
#process.pclusterCompatibilityFilter = cms.Path(process.clusterCompatibilityFilter)
#process.pprimaryVertexFilter = cms.Path(process.primaryVertexFilter)
#process.pAna = cms.EndPath(process.skimanalysis)
