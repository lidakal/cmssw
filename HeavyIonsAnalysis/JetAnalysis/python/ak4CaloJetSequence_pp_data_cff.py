import FWCore.ParameterSet.Config as cms

from HeavyIonsAnalysis.JetAnalysis.caloJetAnalyzer_cff import *

ak4CaloJetAnalyzer = caloJetAnalyzer.clone(
    jetTag = cms.InputTag("slimmedCaloJets"),
    rParam = 0.4,
    fillGenJets = False,
    isMC = False,
    jetName = cms.untracked.string("ak4Calo"),
    hltTrgResults = cms.untracked.string('TriggerResults::'+'HISIGNAL')
)
