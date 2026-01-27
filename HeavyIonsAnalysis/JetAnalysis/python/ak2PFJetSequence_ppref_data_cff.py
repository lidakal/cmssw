import FWCore.ParameterSet.Config as cms

from HeavyIonsAnalysis.JetAnalysis.inclusiveJetAnalyzer_cff import *

ak2PFJetAnalyzer = inclusiveJetAnalyzer.clone(
    jetTag = cms.InputTag("slimmedJets"),
    rParam = 0.2,
    fillGenJets = False,
    isMC = False,
    jetName = cms.untracked.string("ak2PF"),
    hltTrgResults = cms.untracked.string('TriggerResults::'+'HISIGNAL')
    )
