import FWCore.ParameterSet.Config as cms

process = cms.Process('RECLUSTER')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/mc/RunIISummer19UL17MiniAOD/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/100000/BFAAC85A-F5C5-8843-8D2A-76A9E873E24B.root'),
    secondaryFileNames = cms.untracked.vstring(),
                            #skipEvents = cms.untracked.uint32(62),
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('test_dynGroomedJets_cfg.py nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# edm output for debugging purposes
process.output = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('JetClusterEDM.root'),
    outputCommands = cms.untracked.vstring(
        'keep *',
        # drop aliased products
        'drop *_akULPu3PFJets_*_*',
        'drop *_akULPu4PFJets_*_*',
    )
)


# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mc2017_realistic_v8', '')


process.load("GeneratorInterface.RivetInterface.mergedGenParticles_cfi")
process.load("PhysicsTools.JetMCAlgos.HFdecayProductTagger_cfi")
process.HFdecayProductTagger.genParticles = 'mergedGenParticles'
process.load("RecoBTag.ImpactParameter.impactParameter_EventSetup_cff")
process.load("RecoBTag.ImpactParameter.pfImpactParameterTagInfos_cfi")
process.pfImpactParameterTagInfos.candidates  = 'packedPFCandidates'
process.pfImpactParameterTagInfos.primaryVertex = 'offlineSlimmedPrimaryVertices'
process.pfImpactParameterTagInfos.jets = 'slimmedJets'
process.load("RecoBTag.SecondaryVertex.pfInclusiveSecondaryVertexFinderTagInfos_cfi")
process.pfInclusiveSecondaryVertexFinderTagInfos.extSVCollection = 'slimmedSecondaryVertices'
process.load("RecoJets.JetProducers.dynGroomedGenJets_cfi")
process.load("RecoJets.JetProducers.dynGroomedPFJets_cfi")

process.jetCluster = cms.Sequence(
    process.mergedGenParticles+
    process.HFdecayProductTagger+
    process.pfImpactParameterTagInfos+
    process.pfInclusiveSecondaryVertexFinderTagInfos+
    process.dynGroomedGenJets+
    process.dynGroomedPFJets
)

# Path and EndPath definitions
process.ana_step = cms.Path(process.jetCluster)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.output_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(process.ana_step,process.endjob_step,process.output_step)

#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(1)
process.options.numberOfStreams=cms.untracked.uint32(0)
process.options.numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1)

# customisation of the process.


# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
