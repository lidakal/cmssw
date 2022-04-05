

import FWCore.ParameterSet.Config as cms
from HeavyIonsAnalysis.JetAnalysis.patHeavyIonSequences_cff import patJetGenJetMatch, patJetPartonMatch, patJetCorrFactors, patJets
from HeavyIonsAnalysis.JetAnalysis.inclusiveJetAnalyzer_cff import *
from HeavyIonsAnalysis.JetAnalysis.bTaggers_cff import *
from RecoJets.JetProducers.JetIDParams_cfi import *

ak4PFmatch = patJetGenJetMatch.clone(
    src = "ak4PFJets",
    matched = "ak4GenJets",
    resolveByMatchQuality = False,
    maxDeltaR = 0.4
    )

ak4PFmatchGroomed = patJetGenJetMatch.clone(
    src = "ak4GenJets",
    matched = "ak4GenJets",
    resolveByMatchQuality = False,
    maxDeltaR = 0.4
    )

ak4PFparton = patJetPartonMatch.clone(src = "ak4PFJets")

ak4PFcorr = patJetCorrFactors.clone(
    useNPV = False,
    useRho = False,
    levels   = cms.vstring('L2Relative'),
    src = "ak4PFJets",
    payload = "AK4PF"
    )

ak4PFJetID= cms.EDProducer('JetIDProducer', JetIDParams, src = cms.InputTag('ak4CaloJets'))

ak4PFbTagger = bTaggers("ak4PF",0.4,True,False)

ak4PFparton = patJetPartonMatch.clone(src = "ak4PFJets", matched = "genParticles")
ak4PFPatJetFlavourAssociationLegacy = ak4PFbTagger.PatJetFlavourAssociationLegacy
#ak4PFPatJetPartons = ak4PFbTagger.PatJetPartons
ak4PFJetTracksAssociatorAtVertex = ak4PFbTagger.JetTracksAssociatorAtVertex
ak4PFJetTracksAssociatorAtVertex.tracks = "highPurityTracks"
ak4PFCombinedSecondaryVertexV2BJetTags = ak4PFbTagger.CombinedSecondaryVertexV2BJetTags
ak4PFPatJetPartonAssociationLegacy = ak4PFbTagger.PatJetPartonAssociationLegacy

ak4PFImpactParameterTagInfos = ak4PFbTagger.ImpactParameterTagInfos
ak4PFPfImpactParameterTagInfos = ak4PFbTagger.PfImpactParameterTagInfos
ak4PFJetProbabilityBJetTags = ak4PFbTagger.JetProbabilityBJetTags
#ak4PFJetProbabilityBJetTags.jetTagComputer = cms.string('candidateJetProbabilityComputer')
#ak4PFJetProbabilityBJetTags.tagInfos = cms.VInputTag("ak4PFPfImpactParameterTagInfos")


ak4PFSecondaryVertexTagInfos = ak4PFbTagger.SecondaryVertexTagInfos
ak4PFCombinedSecondaryVertexV2BJetTags = ak4PFbTagger.CombinedSecondaryVertexV2BJetTags
ak4PFPfInclusiveSecondaryVertexFinderTagInfos = ak4PFbTagger.PfInclusiveSecondaryVertexFinderTagInfos
ak4PFPfDeepCSVTagInfos = ak4PFbTagger.PfDeepCSVTagInfos
ak4PFPfDeepCSVJetTags = ak4PFbTagger.PfDeepCSVJetTags
ak4PFPrimaryVertexAssociation = ak4PFbTagger.PrimaryVertexAssociation
ak4PFPfDeepFlavourTagInfos = ak4PFbTagger.PfDeepFlavourTagInfos
ak4PFPfDeepFlavourJetTags = ak4PFbTagger.PfDeepFlavourJetTags


ak4PFPatJetFlavourIdLegacy = cms.Sequence(ak4PFPatJetPartonAssociationLegacy*ak4PFPatJetFlavourAssociationLegacy)

#ak4PFPatJetFlavourAssociation = ak4PFbTagger.PatJetFlavourAssociation
ak4PFPatJetFlavourAssociation = cms.EDProducer("JetFlavourClustering",
        bHadrons = cms.InputTag("selectedHadronsAndPartons:bHadrons"),
        cHadrons = cms.InputTag("selectedHadronsAndPartons:cHadrons"),
        ghostRescaling = cms.double(1e-18),
       # groomedJets = cms.InputTag('akSoftDrop4PFJets'),
        hadronFlavourHasPriority = cms.bool(False),
        jetAlgorithm = cms.string('AntiKt'),
        jets = cms.InputTag("ak4PFJets"),
        leptons = cms.InputTag("selectedHadronsAndPartons:leptons"),
        partons = cms.InputTag("selectedHadronsAndPartons:algorithmicPartons"),
        rParam = cms.double(0.4),
        #subjets = cms.InputTag("akSoftDrop4PFJets:SubJets")
)

#ak4PFPatJetFlavourId = cms.Sequence(ak4PFPatJetPartons*ak4PFPatJetFlavourAssociation)
ak4PFPatJetFlavourId = cms.Sequence(ak4PFPatJetFlavourAssociation) # moved partons to global gen stuff (not for each jet algo)

ak4PFJetBtaggingIP       = cms.Sequence(ak4PFImpactParameterTagInfos *
                                        ak4PFPfImpactParameterTagInfos * 
                                        ak4PFJetProbabilityBJetTags 
)

ak4PFJetBtaggingSV = cms.Sequence(ak4PFImpactParameterTagInfos *
    ak4PFSecondaryVertexTagInfos *
    ak4PFCombinedSecondaryVertexV2BJetTags *
    ak4PFPfImpactParameterTagInfos *
    ak4PFPfInclusiveSecondaryVertexFinderTagInfos *
    ak4PFPfDeepCSVTagInfos *
    ak4PFPfDeepCSVJetTags *
    ak4PFPrimaryVertexAssociation *
    ak4PFPfDeepFlavourTagInfos *
    ak4PFPfDeepFlavourJetTags 
)


ak4PFJetBtagging = cms.Sequence(ak4PFJetBtaggingIP
            *ak4PFJetBtaggingSV
            )




ak4PFpatJetsWithBtagging = patJets.clone(jetSource = "ak4PFJets",
        genJetMatch          = "ak4PFmatch",
        genPartonMatch       = "ak4PFparton",
        jetCorrFactorsSource = ["ak4PFcorr"],
        JetPartonMapSource   = "ak4PFPatJetFlavourAssociation",
	JetFlavourInfoSource   = "ak4PFPatJetFlavourAssociation",
        trackAssociationSource = "ak4PFJetTracksAssociatorAtVertex",
	useLegacyJetMCFlavour = False,
        discriminatorSources = [
            "ak4PFCombinedSecondaryVertexV2BJetTags",
            "ak4PFJetProbabilityBJetTags",
            "ak4PFPfDeepCSVJetTags:probb",
            "ak4PFPfDeepCSVJetTags:probbb",
            "ak4PFPfDeepCSVJetTags:probc",
            "ak4PFPfDeepFlavourJetTags:probb",
            "ak4PFPfDeepFlavourJetTags:probbb",
            "ak4PFPfDeepFlavourJetTags:problepb",
            "ak4PFPfDeepFlavourJetTags:probc",
        ],
        tagInfoSources = ["ak4PFPfImpactParameterTagInfos","ak4PFPfInclusiveSecondaryVertexFinderTagInfos","ak4PFPfDeepCSVTagInfos","ak4PFPfDeepFlavourTagInfos"],
        jetIDMap = "ak4PFJetID",
        addBTagInfo = True,
        addTagInfos = True,
        addDiscriminators = True,
        addAssociatedTracks = True,
        addJetCharge = False,
        addJetID = False,
        getJetMCFlavour = True,
        addGenPartonMatch = True,
        addGenJetMatch = True,
        embedGenJetMatch = True,
        embedGenPartonMatch = True,
)



'''
ak4PFSubJetPfImpactParameterTagInfos = cms.EDProducer("CandIPProducer",
    candidates = cms.InputTag("particleFlow"),
    computeGhostTrack = cms.bool(True),
    computeProbabilities = cms.bool(True),
    explicitJTA = cms.bool(True),
    ghostTrackPriorDeltaR = cms.double(0.03),
    jetDirectionUsingGhostTrack = cms.bool(False),
    jetDirectionUsingTracks = cms.bool(False),
    jets = cms.InputTag("akSoftDrop4PFJets","SubJets"),
    maxDeltaR = cms.double(0.2),
    maximumChiSquared = cms.double(5.0),
    maximumLongitudinalImpactParameter = cms.double(17.0),
    maximumTransverseImpactParameter = cms.double(0.2),
    minimumNumberOfHits = cms.int32(0),
    minimumNumberOfPixelHits = cms.int32(1),
    minimumTransverseMomentum = cms.double(1.0),
    primaryVertex = cms.InputTag("offlinePrimaryVertices"),
    useTrackQuality = cms.bool(False)
)

ak4PFSubJetPfInclusiveSecondaryVertexFinderTagInfos = ak4PFPfInclusiveSecondaryVertexFinderTagInfos.clone(
    trackIPTagInfos = cms.InputTag("ak4PFSubJetPfImpactParameterTagInfos"),
    fatJets = cms.InputTag("ak4PFJets"),
    groomedFatJets = cms.InputTag("akSoftDrop4PFJets"),
    rParam = cms.double(0.4),
    jetAlgorithm = cms.string('AntiKt'),
    useSVClustering = cms.bool(True),
)

ak4PFSubJetPfDeepCSVTagInfos = ak4PFPfDeepCSVTagInfos.clone(
    fatJets = cms.InputTag("ak4PFJets"),
    groomedFatJets = cms.InputTag("akSoftDrop4PFJets"),
    jetAlgorithm = cms.string('AntiKt'),
    rParam = cms.double(0.4),
    svTagInfos = cms.InputTag("ak4PFSubJetPfInclusiveSecondaryVertexFinderTagInfos"),
    useSVClustering = cms.bool(True)
)

ak4PFPfDeepCSVSubJetTags = ak4PFPfDeepCSVJetTags.clone(src = "ak4PFSubJetPfDeepCSVTagInfos")
'''
ak4PFpatSubJetsWithBtagging = patJets.clone(jetSource = "akSoftDrop4PFJets:SubJets",
        genJetMatch          = cms.InputTag("NULL"),
        genPartonMatch       = cms.InputTag("NULL"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("NULL")),
                                            #JetPartonMapSource   = "ak4PFPatJetFlavourAssociation:SubJets",
                                            #JetFlavourInfoSource   = "ak4PFPatJetFlavourAssociation:SubJets",
        JetPartonMapSource   = "NULL",
        JetFlavourInfoSource   = "NULL",
        trackAssociationSource = cms.InputTag("NULL"),
	useLegacyJetMCFlavour = False,
        #discriminatorSources = cms.VInputTag(
        #    cms.InputTag("ak4PFPfDeepCSVSubJetTags:probb"),
        #    cms.InputTag("ak4PFPfDeepCSVSubJetTags:probbb"),
        #    cms.InputTag("ak4PFPfDeepCSVSubJetTags:probc"),
        #    cms.InputTag("ak4PFPfDeepCSVSubJetTags:probudsg"),            
                                            #),
        discriminatorSources = cms.VInputTag(),
        #tagInfoSources = cms.VInputTag("ak4PFSubJetPfImpactParameterTagInfos","ak4PFSubJetPfInclusiveSecondaryVertexFinderTagInfos","ak4PFSubJetPfDeepCSVTagInfos"),
        tagInfoSources = cms.VInputTag(),
        jetIDMap = cms.InputTag("NULL"),
        #addBTagInfo = True,
        #addTagInfos = True,
        #addDiscriminators = True,
        addBTagInfo = False,
        addTagInfos = False,
        addDiscriminators = False,
        addAssociatedTracks = False,
        addJetCharge = False,
        addJetID = False,
        getJetMCFlavour = False,
        addGenPartonMatch = False,
        addGenJetMatch = False,
        embedGenJetMatch = False,
        embedGenPartonMatch = False,
        addJetCorrFactors = False
)



ak4PFJetAnalyzer = inclusiveJetAnalyzer.clone(jetTag = cms.InputTag("ak4PFpatJetsWithBtagging"),
                                              genjetTag = 'ak4GenJets',
                                              rParam = 0.4,
                                              matchJets = cms.untracked.bool(False),
                                              matchTag = 'patJetsWithBtagging',
                                              pfCandidateLabel = cms.untracked.InputTag('particleFlow'),
                                              fillGenJets = True,
                                              isMC = True,
                                              genParticles = cms.untracked.InputTag("genParticles"),
                                              doLifeTimeTagging = cms.untracked.bool(True),
                                              bTagJetName = cms.untracked.string("ak4PF"),
                                              jetName = cms.untracked.string("ak4PF"),
                                              genPtMin = cms.untracked.double(20.),
                                              doSubJets = cms.untracked.bool(False),
                                              doGenSubJets = cms.untracked.bool(False),     
                                              subjetGenTag = cms.untracked.InputTag("ak4GenJets"),
                                              doExtendedFlavorTagging = cms.untracked.bool(True),
                                              jetFlavourInfos = cms.InputTag("ak4PFPatJetFlavourAssociation"),
                                              subjetFlavourInfos = cms.InputTag("ak4PFPatJetFlavourAssociation","SubJets"),
                                              groomedJets = cms.InputTag("ak4PFJets"),
                                              jetPtMin = 20.
)



ak4PFJetSequence_mc = cms.Sequence(
    ak4PFmatch*
    ak4PFparton*
    ak4PFcorr*
    ak4PFPatJetFlavourId*  
    ak4PFJetTracksAssociatorAtVertex*                                                 
    ak4PFJetBtagging*
    ak4PFpatJetsWithBtagging*
    #ak4PFPatSubJetFlavourAssociation*
    #ak4PFSubJetPfImpactParameterTagInfos *
    #ak4PFSubJetPfInclusiveSecondaryVertexFinderTagInfos *
    #ak4PFSubJetPfDeepCSVTagInfos *
    #ak4PFPfDeepCSVSubJetTags *                                                                        
    ak4PFpatSubJetsWithBtagging*
    ak4PFJetAnalyzer)

ak4PFJetSequence_data = cms.Sequence(ak4PFcorr*
                                     ak4PFJetTracksAssociatorAtVertex*
                                     ak4PFJetBtagging*
                                     ak4PFpatJetsWithBtagging*
                                     ak4PFJetAnalyzer)

ak4PFJetSequence_jec = cms.Sequence(ak4PFJetSequence_mc)
ak4PFJetSequence_mb = cms.Sequence(ak4PFJetSequence_mc)

ak4PFJetSequence = cms.Sequence(ak4PFJetSequence_mc)

### extra substructure stuff
ak4PFJetAnalyzer.doSubJets = True
ak4PFJetAnalyzer.doGenSubJets = True
ak4PFJetAnalyzer.groomedJets = 'akSoftDrop4PFJets'
ak4PFJetAnalyzer.groomedGenJets = cms.InputTag('akSoftDrop4GenJets')
