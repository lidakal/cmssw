// -*- C++ -*-
//
// Package:    GenHFHadronReplacer
// Class:      GenHFHadronReplacer
//

/**\class GenHFHadronReplacer GenHFHadronReplacer.cc
* @brief Given a list of heavy flavour hadrons, produce a genParticle collection replacing their decay products
*
* optionally use a pseudo-particle for the HF hadron, composed of only its charged decay products 
*
* The description of the run-time parameters can be found at fillDescriptions()
*
* The description of the products can be found at GenHFHadronReplacer()
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



class GenHFHadronReplacer : public edm::global::EDProducer<> {
public:
  explicit GenHFHadronReplacer(const edm::ParameterSet &);
  ~GenHFHadronReplacer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
  
private:
  void produce(edm::StreamID, edm::Event &, const edm::EventSetup &) const override;
  
  // ----------member data ---------------------------
  const edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  bool isFinalB( const reco::Candidate &particle) const;  
  bool isFromB( const reco::Candidate &particle) const;  
  void visible( reco::Candidate::PolarLorentzVector &v, const reco::Candidate &particle, bool doCharge) const;
};


GenHFHadronReplacer::GenHFHadronReplacer(const edm::ParameterSet& cfg)
  : genParticlesToken_(consumes<reco::GenParticleCollection>(cfg.getParameter<edm::InputTag>("genParticles"))){
  produces<reco::GenParticleCollection>();
		 }


// ------------ method called to produce the data  ------------
void GenHFHadronReplacer::produce(edm::StreamID, edm::Event &evt, const edm::EventSetup &setup) const {


  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByToken(genParticlesToken_, genParticles);

  auto outputCollection = std::make_unique<reco::GenParticleCollection>();

  for (const reco::GenParticle& genPart : *genParticles) {
    if(isFinalB(genPart)){
      //std::cout<<" b hadron = "<<genPart.pt()<<" "<<genPart.eta()<<" "<<genPart.phi()<<" "<<genPart.mass()<<" "<<genPart.status()<<" "<<genPart.px()<<" "<<genPart.p()<<std::endl;
      reco::Candidate::PolarLorentzVector chargedB(0.,0.,0.,0.);
      reco::Candidate::PolarLorentzVector neutralB(0.,0.,0.,0.);
      visible(chargedB, genPart,true);
      visible(neutralB, genPart,false);
      reco::GenParticle newGenCharged(1., chargedB, genPart.vertex(), genPart.pdgId(), 1, true);
      reco::GenParticle newGenNeutral(0., neutralB, genPart.vertex(), genPart.pdgId(), 1, true);
      //std::cout<<" charged  = "<<newGenCharged.pt()<<" "<<newGenCharged.eta()<<" "<<newGenCharged.phi()<<" "<<newGenCharged.mass()<<" "<<newGenCharged.status()<<" "<<newGenCharged.px()<<" "<<newGenCharged.p()<<std::endl;
      //std::cout<<" neutral  = "<<newGenNeutral.pt()<<" "<<newGenNeutral.eta()<<" "<<newGenNeutral.phi()<<" "<<newGenNeutral.mass()<<" "<<newGenNeutral.status()<<" "<<newGenNeutral.px()<<" "<<newGenNeutral.p()<<std::endl;
      outputCollection->push_back(newGenCharged);
      outputCollection->push_back(newGenNeutral);
    }
    if(genPart.status()!=1) continue;
    if(isFromB(genPart)) continue;
    outputCollection->push_back(reco::GenParticle(genPart.charge(), genPart.p4(), genPart.vertex(), genPart.pdgId(), genPart.status(), true));
  }
  
  evt.put(std::move(outputCollection));
}

bool GenHFHadronReplacer::isFinalB( const reco::Candidate &particle) const{
  if (!CandMCTagUtils::hasBottom(particle) ) return false;

  // check if any of the daughters is also a b hadron
    unsigned int npart=particle.numberOfDaughters();
    
    for (size_t i = 0; i < npart; ++i) {
      if (CandMCTagUtils::hasBottom(*particle.daughter(i))) return false;
    }
    
    return true;
}


bool GenHFHadronReplacer::isFromB( const reco::Candidate &particle) const{

  unsigned int npart=particle.numberOfMothers();
  for (size_t i = 0; i < npart; ++i) {
    const reco::Candidate &mom = *particle.mother(i);
    if (CandMCTagUtils::hasBottom(mom)) return true;
    isFromB(mom);
  }
  return false;
}

void GenHFHadronReplacer::visible
(
 reco::Candidate::PolarLorentzVector &v, const reco::Candidate &particle, bool doCharge) const
{


  unsigned int npart=particle.numberOfDaughters();
  
  for (size_t i = 0; i < npart; ++i) {
    if(particle.daughter(i)->status()==1){
      int charge = particle.daughter(i)->charge();
      if(doCharge && charge ==0) continue;
      int pdgid = abs(particle.daughter(i)->pdgId());
      if(!doCharge && (charge !=0 || pdgid == 12 || pdgid == 14 || pdgid == 16) ) continue;     

      reco::Candidate::PolarLorentzVector vTemp(0.,0.,0.,0.);
      vTemp.SetPt(particle.daughter(i)->pt());
      vTemp.SetEta(particle.daughter(i)->eta());
      vTemp.SetPhi(particle.daughter(i)->phi());
      vTemp.SetM(particle.daughter(i)->mass());
      v+=vTemp;

    }
    else{
      visible(v,*particle.daughter(i),doCharge);
    }
  }
  
}



void GenHFHadronReplacer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("genParticles", edm::InputTag("genParticles"));
  descriptions.add("GenHFHadronReplacer", desc);
}


DEFINE_FWK_MODULE(GenHFHadronReplacer);


