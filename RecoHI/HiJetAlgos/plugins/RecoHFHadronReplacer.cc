// -*- C++ -*-
//
// Package:    RecoHFHadronReplacer
// Class:      RecoHFHadronReplacer
//

/**\class RecoHFHadronReplacer RecoHFHadronReplacer.cc
* @brief Given a list of heavy flavour hadrons, produce a genParticle collection replacing their decay products
*
* optionally use a pseudo-particle for the HF hadron, composed of only its charged decay products 
*
* The description of the run-time parameters can be found at fillDescriptions()
*
* The description of the products can be found at RecoHFHadronReplacer()
*/

// Original Author:  Matthew Nguyen, LLR


// system include files
#include <memory>
#include <utility>
#include <vector>
#include <algorithm>

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

#include "CommonTools/CandUtils/interface/pdgIdUtils.h"
#include "PhysicsTools/JetMCUtils/interface/CandMCTag.h"
#include "TLorentzVector.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/ShallowTagInfo.h"
#include "CommonTools/UtilAlgos/interface/DeltaR.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "DataFormats/TrackReco/interface/Track.h"



class RecoHFHadronReplacer : public edm::global::EDProducer<> {
public:
  explicit RecoHFHadronReplacer(const edm::ParameterSet &);
  ~RecoHFHadronReplacer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
  
private:
  void produce(edm::StreamID, edm::Event &, const edm::EventSetup &) const override;
  
  // ----------member data ---------------------------
  const edm::EDGetTokenT<edm::View<reco::BaseTagInfo>> tagInfoToken_;
};


RecoHFHadronReplacer::RecoHFHadronReplacer(const edm::ParameterSet& cfg)
  : tagInfoToken_(consumes<edm::View<reco::BaseTagInfo>>(cfg.getParameter<edm::InputTag>("ak4PFPfImpactParameterTagInfos"))){
  produces<reco::GenParticleCollection>();
}


// ------------ method called to produce the data  ------------
void RecoHFHadronReplacer::produce(edm::StreamID, edm::Event &evt, const edm::EventSetup &setup) const {


  edm::Handle<edm::View<reco::BaseTagInfo>> ipTagInfo;
  evt.getByToken(tagInfoToken_, ipTagInfo);

  auto outputCollection = std::make_unique<reco::GenParticleCollection>();

  for (auto const& ipJet : *ipTagInfo) {
    std::cout<<ipJet.jet()->eta()<< ipJet.jet()->phi() <<std::endl;

    const reco::IPTagInfo<std::vector<reco::CandidatePtr>, reco::JetTagInfo>* tagInfo =
      dynamic_cast<const reco::IPTagInfo<std::vector<reco::CandidatePtr>, reco::JetTagInfo>*>(&ipJet);
      
    //const std::vector<reco::btag::TrackIPData>& ip = tagInfo->impactParameterData();
    
    for (const auto & trk : tagInfo->selectedTracks()){
      std::cout<<"track: "<<trk->pt()<<" "<<trk->eta()<<std::endl;
      
    }
  }
  
  evt.put(std::move(outputCollection));
}



void RecoHFHadronReplacer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("ak4PFPfImpactParameterTagInfos", edm::InputTag("ak4PFPfImpactParameterTagInfos"));
  descriptions.add("RecoHFHadronReplacer", desc);
}


DEFINE_FWK_MODULE(RecoHFHadronReplacer);


