### HiForest Configuration
# Input: miniAOD
# Type: mc
import sys

# fname = sys.argv[2]
# print("file in:", fname)

import FWCore.ParameterSet.Config as cms
process = cms.Process('HiForest')
process.options = cms.untracked.PSet()

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
        # '/store/himc/RunIISummer20UL17pp5TeVMiniAODv2/QCD_pThat-15_bJet_TuneCP5_5p02TeV-pythia8/MINIAODSIM/106X_mc2017_realistic_forppRef5TeV_v3-v3/2530000/071C3B54-D788-DB4A-8FDD-FBCE0A911D55.root'
        # '/store/himc/RunIISummer20UL17pp5TeVMiniAODv2/QCD_pThat-15_Dijet_TuneCP5_5p02TeV-pythia8/MINIAODSIM/106X_mc2017_realistic_forppRef5TeV_v3-v3/50000/ED89F8C3-56E9-AF48-9D64-7D9CF448BE0C.root'
        # '/store/mc/RunIISummer20UL16MiniAODAPVv2/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_preVFP_v11-v1/120000/0230C4F8-6445-F74C-8409-27F77DFDE107.root'
        # '/store/mc/RunIILowPUSummer20UL17MiniAODv2/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/pilot_106X_mc2017_realistic_v9For2017H_v1-v2/2560000/D9746C19-3FD6-B246-95F4-E7DBD33ED768.root'
        # '/store/mc/RunIISummer20UL17MiniAODv2/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_mc2017_realistic_v9-v1/240000/DC372BF6-08B2-9C4A-AF9E-69E80E066E4F.root'
        # '/store/mc/RunIISummer20UL17MiniAODv2/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_mc2017_realistic_v9-v1/240000/003621AC-81B6-444A-A43A-E5D0921C7356.root'
        # "/store/himc/RunIISummer20UL17pp5TeVMiniAODv2/QCD_pThat-15_Dijet_TuneCP5_5p02TeV-pythia8/MINIAODSIM/106X_mc2017_realistic_forppRef5TeV_v3-v3/40000/7DC3DCCF-F5E0-C04B-ABD9-21E6631E185D.root"
        # "file:/data_CMS/cms/mnguyen//forLida/miniAOD_PAT_94X.root"
        # 'file:/home/llr/cms/kalipoliti/rootFiles/071C3B54-D788-DB4A-8FDD-FBCE0A911D55.root'
        # '/store/user/mnguyen/Herwig_CH3_bjet_5TeV/Herwig_CH3_bjet_5TeV_MINI_v7/240123_124848/0000/mini_PAT_1.root' # herwig bjet
        # '/store/user/mnguyen/Herwig_CH3_qcd/Herwig_CH3_qcd_MINIAOD/231017_075304/0000/recopat_RAW2DIGI_L1Reco_RECO_RECOSIM_EI_PAT_1.root' # herwig qcd
        # '/store/data/Run2017G/LowEGJet/MINIAOD/17Nov2017-v1/100000/0048A471-EE2D-E811-BB4D-0CC47AD98F70.root'
        # '/store/user/mnguyen/Herwig_CH3_qcd_5TeV/Herwig_CH3_qcd_5TeV_MINI_v7/240110_112918/0000/mini_PAT_1.root'
        # '/store/himc/RunIISummer20UL17pp5TeVMiniAODv2/QCD_pThat-15_Mujet_TuneCP5_5p02TeV-pythia8/MINIAODSIM/106X_mc2017_realistic_forppRef5TeV_v3-v2/2530000/098E0FF6-2717-EF4E-9E93-CC946DB97402.root'
        # 'file:/data_CMS/cms/mnguyen///mini_PAT.root'
        # fname
        'file:/data_CMS/cms/kalipoliti/8B655CAA-8B9D-E04F-AE94-54227BD26D96.root'
        ),
    )

# Select specific event
# process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange('1:21907')
# process.source.eventsToProcess = cms.untracked.VEventRange('1:21906437')

# number of events to process, set to -1 to process all events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
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
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mc2017_realistic_v10', '')
process.HiForestInfo.GlobalTagLabel = process.GlobalTag.globaltag
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")

process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
             tag = cms.string("JPcalib_MC94X_2017pp_v2"), # mc tag
            # tag = cms.string("JPcalib_Data4MCSwap_94X_2017pp_v3"), # data tag for swap 
            connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS")
         )
      ])

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
# )
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
process.load('HeavyIonsAnalysis.EventAnalysis.skimanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hltobject_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.l1object_cfi')

from HeavyIonsAnalysis.EventAnalysis.hltobject_cfi import trigger_list_mc
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

#################################
# rho
process.load("RecoJets.JetProducers.fixedGridRhoProducerFastjet_cfi")
process.fixedGridRhoFastjetAll.pfCandidatesTag = cms.InputTag("packedPFCandidates")
process.rhoSequence = cms.Sequence(
    process.fixedGridRhoFastjetAll
    # process.hiFJGridEmptyAreaCalculator +
    # process.hiFJRhoProducer +
    # process.hiFJRhoAnalyzer
)

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
    # process.HiGenParticleAna + // for dN/deta checks
    # process.updatePATJetSequence + 
    process.tagInfoSequence +
    process.genJetSequence + 
    process.recoJetSequence + 
    #process.unpackedMuons +
    #process.correctedElectrons #+
    # process.ggHiNtuplizer +
    #process.muonAnalyzer + 
    process.rhoSequence +
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

    # Rerun JetProbability
    process.load("RecoBTag.ImpactParameter.pfJetProbabilityBJetTags_cfi")
    process.pfJetProbabilityBJetTags = process.pfJetProbabilityBJetTags.clone()

    from RecoBTag.ONNXRuntime.pfParticleNetAK4_cff import _pfParticleNetAK4JetTagsAll as pfParticleNetAK4JetTagsAll

    from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
    updateJetCollection(
        process,
        jetSource = cms.InputTag('slimmedJets'),
        jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
        btagDiscriminators = ['jetProbabilityBJetTags'],
        btagPrefix = "rerun"
        )
    process.updatedPatJets.addJetCorrFactors = False   
    process.updatedPatJets.addTagInfos = True
    process.updatedPatJets.addBTagInfo = cms.bool(True) ## Needed to add tag infos
    process.updatedPatJets.tagInfoSources = cms.VInputTag(
        cms.InputTag(ipTagInfoLabel_ + "TagInfos"),
        cms.InputTag(svTagInfoLabel_ + "TagInfos"),
    )
    process.updatedPatJets.discriminatorSources = cms.VInputTag("pfJetProbabilityBJetTags")

    process.tagInfoSequence.insert(0, process.pfImpactParameterTagInfos *
                                   process.pfInclusiveSecondaryVertexFinderTagInfos *
                                   process.pfJetProbabilityBJetTags *
                                   process.updatedPatJets)
    
    process.ak4PFJetAnalyzer.jetTag = cms.InputTag("updatedPatJets")
    process.ak4PFJetAnalyzer.rhoSrc = cms.InputTag("fixedGridRhoFastjetAll")

    ## Number of b and c hadrons
    from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartons
    process.selectedHadronsAndPartons = selectedHadronsAndPartons.clone(
        particles = cms.InputTag("mergedGenParticles")
        )   
    from PhysicsTools.JetMCAlgos.AK4PFJetsMCFlavourInfos_cfi import ak4JetFlavourInfos
    process.jetFlavourInfosAK4PFJets = ak4JetFlavourInfos.clone(
        jets = cms.InputTag("slimmedJets")
    )

    process.recoJetSequence += (process.selectedHadronsAndPartons * process.jetFlavourInfosAK4PFJets)
    process.ak4PFJetAnalyzer.jetFlavourInfos = cms.InputTag("jetFlavourInfosAK4PFJets")


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
doLatekt_ = False

tmva_variables = ["trkIp3dSig", "trkIp2dSig", "trkDistToAxis",
                  "svtxdls", "svtxdls2d", "svtxm", "svtxmcorr",
                  "svtxnormchi2", "svtxNtrk", "svtxTrkPtOverSv",
                  "jtpt"]
if doDeclustering:
    process.load("GeneratorInterface.RivetInterface.mergedGenParticles_cfi")
    process.genJetSequence += process.mergedGenParticles
    ## Produces a reco::GenParticleCollection named mergedGenParticles

    process.load("RecoHI.HiJetAlgos.HFdecayProductTagger_cfi")
    process.HFdecayProductTagger.genParticles = cms.InputTag("mergedGenParticles")
    process.HFdecayProductTagger.tagBorC = cms.bool(True) # tag B
    process.genJetSequence += process.HFdecayProductTagger
    taggedGenParticlesName_ = "HFdecayProductTagger"
    ## Produces a std::vector<pat::PackedGenParticle> named HFdecayProductTagger

    process.ak4PFJetAnalyzer.genParticles = cms.untracked.InputTag(taggedGenParticlesName_, "patPackedGenParticles")

    process.bDecayAna = process.HiGenParticleAna.clone(
        genParticleSrc = cms.InputTag(taggedGenParticlesName_, "patPackedGenParticles"),
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
    process.TrackToGenParticleMapProducer.genParticleSrc = cms.InputTag(taggedGenParticlesName_, "patPackedGenParticles")
    process.TrackToGenParticleMapProducer.chargedOnly = doChargedOnly
    process.genJetSequence += process.TrackToGenParticleMapProducer
    ## Creates the genConstitToGenParticleMap and trackToGenParticleMap

    process.load("RecoHI.HiJetAlgos.dynGroomedPATJets_cfi")
    process.dynGroomedGenJets = process.dynGroomedPATJets.clone(
        chargedOnly = cms.bool(doChargedOnly),
        aggregateHF = cms.bool(doAggregation),
        # aggregateHF = cms.bool(True),
        jetSrc = cms.InputTag("updatedPatJets"),
        constitSrc = cms.InputTag("packedGenParticles"),
        doGenJets = cms.bool(True),
        candToGenParticleMap = cms.InputTag("TrackToGenParticleMapProducer", "genConstitToGenParticleMap"),
        doLateKt = cms.bool(doLatekt_),
    )
    process.genJetSequence += process.dynGroomedGenJets
    process.ak4PFJetAnalyzer.groomedGenJets = cms.untracked.InputTag("dynGroomedGenJets")
    ## Creates the gen jet subjets

    process.dynGroomedPFJets = process.dynGroomedPATJets.clone(
        chargedOnly = cms.bool(doChargedOnly),
        aggregateHF = cms.bool(doAggregation),
        # aggregateHF = cms.bool(False),
        jetSrc = cms.InputTag("updatedPatJets"),
        constitSrc = cms.InputTag("packedPFCandidates"),
        doGenJets = cms.bool(False),
        candToGenParticleMap = cms.InputTag("TrackToGenParticleMapProducer", "trackToGenParticleMap"),
        aggregateWithTruthInfo = cms.bool(False),
        aggregateWithXGB = cms.bool(False),
        aggregateWithTMVA = cms.bool(True),
        aggregateWithCuts = cms.bool(False),
        xgb_path = cms.FileInPath("RecoHI/HiJetAlgos/data/sig_vs_bkg.model"),
        tmva_path = cms.FileInPath("RecoHI/HiJetAlgos/data/TMVAClassification_BDTG.weights.xml"),
        tmva_variables = cms.vstring(tmva_variables),
        doLateKt = cms.bool(doLatekt_),
        trkInefRate = cms.double(0.)
    )
    process.recoJetSequence += process.dynGroomedPFJets
    process.ak4PFJetAnalyzer.groomedJets = cms.untracked.InputTag("dynGroomedPFJets")
    ## creates the reco subjets
    
#########################
# Jet Selection
#########################
    
# for b tagging SF
# process.mujetSelector = cms.EDFilter("PatJetXSelector",
#                              src = cms.InputTag("slimmedJets"),
#                              offPV = cms.InputTag("offlineSlimmedPrimaryVertices"),
#                              cut = cms.string("pt > 5.0 && abs(rapidity()) < 3."),
#                              dummy = cms.bool(False)
#                          )
# process.recoJetSequence += process.mujetSelector
# process.ak4PFJetAnalyzer.mujetTag = cms.InputTag("mujetSelector")

#########################
# Event Selection -> add the needed filters here
#########################
    


#process.load('HeavyIonsAnalysis.EventAnalysis.collisionEventSelection_cff')
#process.pclusterCompatibilityFilter = cms.Path(process.clusterCompatibilityFilter)
#process.pprimaryVertexFilter = cms.Path(process.primaryVertexFilter)
#process.pAna = cms.EndPath(process.skimanalysis)
