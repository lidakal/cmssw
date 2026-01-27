import FWCore.ParameterSet.Config as cms

from HeavyIonsAnalysis.JetAnalysis.inclusiveJetAnalyzer_cff import *

ak3PFJetAnalyzer = inclusiveJetAnalyzer.clone(
    jetTag = cms.InputTag("slimmedJets"),
    rParam = 0.3,
    fillGenJets = False,
    isMC = False,
    jetName = cms.untracked.string("ak3PF"),
    hltTrgResults = cms.untracked.string('TriggerResults::'+'HISIGNAL')
    )
