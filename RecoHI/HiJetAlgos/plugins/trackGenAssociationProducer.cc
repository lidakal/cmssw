// WIP 

// -*- C++ -*-
//
// Package:    trackGenAssociationProducer
// Class:      trackGenAssociationProducer
//

/**\class trackGenAssociationProducer trackGenAssociationProducer.cc
* @brief Produce a mapping between the tracks and the gen particles in the event
* given the IPTagInfo and a list of final state generated particles
*
* The description of the run-time parameters can be found at fillDescriptions()
* The description of the products can be found at trackGenAssociationProducer()
*/

// Original Author:  Lida Kalipoliti, LLR


// system include files
#include <memory>
#include <utility>
#include <vector>
#include <algorithm>
#include <unordered_map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/EDPutToken.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "CommonTools/CandUtils/interface/pdgIdUtils.h"
#include "PhysicsTools/JetMCUtils/interface/CandMCTag.h"
#include "TLorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/PatCandidates/interface/Jet.h"


class trackGenAssociationProducer : public edm::global::EDProducer<> {
public:
  explicit trackGenAssociationProducer(const edm::ParameterSet& cfg);
  ~trackGenAssociationProducer() override = default;

  void produce(edm::StreamID, edm::Event &, const edm::EventSetup &) const override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  // Definitions 
  typedef std::unordered_map<const reco::Candidate *, const pat::PackedGenParticle *> trackGenMap;
  // Configurables
  edm::EDGetTokenT<std::vector<pat::PackedGenParticle>> genParticlesToken_;
  edm::EDGetTokenT<std::vector<pat::Jet>> inputJetsToken_;

  double maxDR2_;
  double minRelPt_;
  double maxRelPt_;
};

trackGenAssociationProducer::trackGenAssociationProducer(const edm::ParameterSet& cfg) {
  // Initialize configurables
  inputJetsToken_ = consumes<std::vector<pat::Jet>>(cfg.getParameter<edm::InputTag>("jetSrc"));
  genParticlesToken_ = consumes<std::vector<pat::PackedGenParticle>>(cfg.getParameter<edm::InputTag>("genParticleSrc"));
  maxDR2_ = cfg.getUntrackedParameter<double>("maxDR2", 0.0004);
  minRelPt_ = cfg.getUntrackedParameter<double>("minRelPt", 0.8);
  maxRelPt_ = cfg.getUntrackedParameter<double>("maxRelPt", 1.2);

  produces<std::unique_ptr<trackGenMap>>("trackGenMap");
}


// ------------ method called to produce the data  ------------
void trackGenAssociationProducer::produce(edm::StreamID, edm::Event &evt, const edm::EventSetup &setup) const {
  //std::cout << "Event being processed..." << std::endl;
  std::cout << "trackGenAssociationProducer produce" << std::endl;

  // Grab the jets
  edm::Handle<std::vector<pat::Jet>> inputJets;
  evt.getByToken(inputJetsToken_, inputJets);

  // Grab the gen particles
  edm::Handle<std::vector<pat::PackedGenParticle>> genParticles;
  evt.getByToken(genParticlesToken_, genParticles);

  std::cout << "nb of input particles, trackGenAssociationProducer: " << (*genParticles).size() << std::endl;

  // Create output collection
  std::unique_ptr<trackGenMap> trackGenMapPtr = std::make_unique<trackGenMap>();

  // Go over the jets
  for (const pat::Jet& jet : *inputJets) {
    // Go over the jet constituents
    const std::vector<reco::CandidatePtr> constituents = jet.getJetConstituents();
    
    for (const reco::CandidatePtr constituent : constituents) {
      if (constituent->charge() == 0) continue;

      // [DEBUG]: 
      // std::cout << "edm::Ptr<reco::Candidate> " << consituent << std::endl;
      // std::cout << "reco::Candidate * " << consituent << std::endl;


      // Look for match in gen particles 
      double minDR2 = std::numeric_limits<double>::max();

      for (const pat::PackedGenParticle& genParticle : *genParticles) {
        if (genParticle.charge() == 0) continue;

        double DR2 = reco::deltaR2(*constituent, genParticle);
        if (DR2 > maxDR2_) continue;

        double relPt = (*constituent).pt() / genParticle.pt();
        if (relPt < minRelPt_ || relPt > maxRelPt_) continue;

        if (DR2 < minDR2) {
          minDR2 = DR2;
          (*trackGenMapPtr)[constituent.get()] = &genParticle;
        }
      } // end gen particle loop
    } // end consituent loop   
  } // end jet loop

  std::cout << "entries in map: " << (*trackGenMapPtr).size() << std::endl;
  
  evt.put(std::move(trackGenMapPtr));
}

void trackGenAssociationProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setComment("track to gen particle association map");
  desc.add<edm::InputTag>("jetSrc", edm::InputTag("slimmedJets"));
  desc.add<edm::InputTag>("genParticleSrc", edm::InputTag("packedGenParticles"));

  desc.addUntracked<double>("maxDR2", 0.0004);
  desc.addUntracked<double>("minRelPt", 0.8);
  desc.addUntracked<double>("maxRelPt", 1.2);

  descriptions.add("trackGenAssociationProducer", desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(trackGenAssociationProducer);
