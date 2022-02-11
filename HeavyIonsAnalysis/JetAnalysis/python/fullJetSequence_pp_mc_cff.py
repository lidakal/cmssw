import FWCore.ParameterSet.Config as cms

from HeavyIonsAnalysis.JetAnalysis.rerecoGen_cff import *
from HeavyIonsAnalysis.JetAnalysis.rerecoJets_cff import *
from HeavyIonsAnalysis.JetAnalysis.rerecoTracks_cff import *

#from HeavyIonsAnalysis.JetAnalysis.jets.ak4CaloJetSequence_pp_mc_cff import *

#from HeavyIonsAnalysis.JetAnalysis.jets.ak3PFJetSequence_pp_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.ak4PFSubJetSequence_pp_mc_cff import *
#from HeavyIonsAnalysis.JetAnalysis.jets.ak5PFJetSequence_pp_mc_cff import *
#from HeavyIonsAnalysis.JetAnalysis.jets.akSoftDrop4PFJetSequence_pp_mc_cff import *

genParticlesForJets.ignoreParticleIDs += [12,14,16]
genPartonsForJets = genParticlesForJets.clone(partonicFinalState = True)
# not clear whether status=1 are turned off by default w/ partonicFinalState = true
genPartonsForJets.ignoreParticleIDs += [11,13,22,130,211,310,321,2212,2112,2114]

ak4PartonJets = ak4GenJets.clone(src = 'genPartonsForJets')
akSoftDrop4PartonJets = akSoftDrop4GenJets.clone(src = 'genPartonsForJets')

genJetSequence = cms.Sequence(
    genParticlesForJets +
    genPartonsForJets +
    #ak3GenJets +
    ak4GenJets +  # need to be reclustered to drop nu's
    ak4PartonJets + 
    akSoftDrop4GenJets +
    akSoftDrop4PartonJets #+
    #ak3GenNjettiness +
    #ak4GenNjettiness
)

jetSequence = cms.Sequence(
    # ak4CaloJets +

    #ak3PFJets +
    # ak4PFJets +
    akSoftDrop4PFJets +
    highPurityTracks +
    ak4PFJetSequence #+
    #akSoftDrop4PFJetSequence 
   #ak4CaloJetSequence +

    #ak3PFJetSequence +
)
