// -*- C++ -*-
//
// Package:    CheatHFHadronReplacer
// Class:      CheatHFHadronReplacer
//

/**\class CheatHFHadronReplacer CheatHFHadronReplacer.cc
* @brief Given a list of heavy flavour hadrons, produce a genParticle collection replacing their decay products
*
* optionally use a pseudo-particle for the HF hadron, composed of only its charged decay products 
*
* The description of the run-time parameters can be found at fillDescriptions()
*
* The description of the products can be found at CheatHFHadronReplacer()
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



class CheatHFHadronReplacer : public edm::global::EDProducer<> {
public:
  explicit CheatHFHadronReplacer(const edm::ParameterSet &);
  ~CheatHFHadronReplacer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
  
private:
  void produce(edm::StreamID, edm::Event &, const edm::EventSetup &) const override;
  
  // ----------member data ---------------------------
  const edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  bool isFinalB(const reco::Candidate &particle) const;  
  bool isFromB(const reco::Candidate &particle) const;  
  std::unique_ptr<reco::GenParticleCollection> addDaughters(const reco::Candidate &particle, int bCode) const;
  //void visible( reco::Candidate::PolarLorentzVector &v, const reco::Candidate &particle, int doCharge) const;
};

CheatHFHadronReplacer::CheatHFHadronReplacer(const edm::ParameterSet& cfg)
  : genParticlesToken_(consumes<reco::GenParticleCollection>(cfg.getParameter<edm::InputTag>("genParticles"))){
  produces<reco::GenParticleCollection>();
		 }


// ------------ method called to produce the data  ------------
void CheatHFHadronReplacer::produce(edm::StreamID, edm::Event &evt, const edm::EventSetup &setup) const {
  std::cout << "Event being processed..." << std::endl;

  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByToken(genParticlesToken_, genParticles);
    
  // Create output collection
  auto outputCollection = std::make_unique<reco::GenParticleCollection>();

  // Kep track of B's
  int bCode = 100;

  for (const reco::GenParticle& genPart : *genParticles) {
    if(isFinalB(genPart)){
      // Add the B
      outputCollection->push_back(genPart);

      // Add the daughters to the output collection
      auto daughterCollection = addDaughters(genPart, bCode);
      (*outputCollection).insert((*outputCollection).end(), (*daughterCollection).begin(), (*daughterCollection).end());
      bCode += 100;
    }

    // Skip b daughters
    if (isFromB(genPart)) continue;	

    // Add the rest of the particles			 
    if (genPart.status()!=1) continue;
    outputCollection->push_back(genPart);
  }
  evt.put(std::move(outputCollection));
}

bool CheatHFHadronReplacer::isFinalB(const reco::Candidate &particle) const {
  if (!CandMCTagUtils::hasBottom(particle) ) return false;

  // check if any of the daughters is also a b hadron
    unsigned int npart=particle.numberOfDaughters();
    
    for (size_t i = 0; i < npart; ++i) {
      if (CandMCTagUtils::hasBottom(*particle.daughter(i))) return false;
    }
    
    return true;
}


bool CheatHFHadronReplacer::isFromB( const reco::Candidate &particle) const{

  bool fromB = false;

  unsigned int npart=particle.numberOfMothers();
  for (size_t i = 0; i < npart; ++i) {
    const reco::Candidate &mom = *particle.mother(i);
    if (CandMCTagUtils::hasBottom(mom)){
      if(isFinalB(mom)){
	fromB = true;
	break;
      }
      else{
	fromB = false;
	break;	
      }
    }
    else fromB = isFromB(mom);
  }
  return fromB;
}
/*
void CheatHFHadronReplacer::visible
(
 reco::Candidate::PolarLorentzVector &v, const reco::Candidate &particle, int doCharge) const
{


  unsigned int npart=particle.numberOfDaughters();
  
  for (size_t i = 0; i < npart; ++i) {
    if(particle.daughter(i)->status()==1){
      int charge = particle.daughter(i)->charge();
      if(doCharge == 1 && charge ==0) continue;
      int pdgid = abs(particle.daughter(i)->pdgId());
      if(doCharge == 0 && charge !=0) continue;
      if(pdgid == 12 || pdgid == 14 || pdgid == 16) continue;     

      reco::Candidate::PolarLorentzVector vTemp(0.,0.,0.,0.);
      vTemp.SetPt(particle.daughter(i)->pt());
      vTemp.SetEta(particle.daughter(i)->eta());
      vTemp.SetPhi(particle.daughter(i)->phi());
      vTemp.SetM(particle.daughter(i)->mass());
      v+=vTemp;
      //std::cout<<" adding a particle:  pt= "<<particle.daughter(i)->pt()<<" pdg id "<<particle.daughter(i)->pdgId()<<std::endl;
    }
    else{
      visible(v,*particle.daughter(i),doCharge);
    }
  }
  
}
*/
std::unique_ptr<reco::GenParticleCollection> CheatHFHadronReplacer::addDaughters(const reco::Candidate &particle, int bCode) const {
  std::cout << "Found a B!" << std::endl;
  // Create collection to return
  static reco::GenParticleCollection daughterCollection;
  //static auto daughterCollection = std::make_unique<reco::GenParticleCollection>();

  // Get daughters
  unsigned int ndaught = particle.numberOfDaughters();

  for (size_t i = 0; i < ndaught; i++) {
    const reco::Candidate& daughter = *particle.daughter(i);
    daughterCollection.push_back(reco::GenParticle(daughter.charge(), daughter.p4(), daughter.vertex(), daughter.pdgId(), bCode, true));
    unsigned int ndaughtdaught = daughter.numberOfDaughters();
    if (ndaughtdaught > 0) {
      addDaughters(daughter, bCode);
    }
  }
  return std::make_unique<reco::GenParticleCollection>(daughterCollection);
}
			    

void CheatHFHadronReplacer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("genParticles", edm::InputTag("genParticles"));
  descriptions.add("CheatHFHadronReplacer", desc);
}


DEFINE_FWK_MODULE(CheatHFHadronReplacer);


