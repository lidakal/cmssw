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
  reco::GenParticleCollection addDaughters(const reco::Candidate &particle, reco::GenParticleCollection daughterCollection, int bCode) const;
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
  int fromBs = 0;
  for (const reco::GenParticle& genPart : *genParticles) {
    if(isFinalB(genPart)){
      // Add the B
      outputCollection->push_back(genPart);
      std::cout << "Found a B!" << std::endl;
      // Add the daughters to the output collection
      reco::GenParticleCollection daughterCollection = {};
      daughterCollection = addDaughters(genPart, daughterCollection, bCode);
      std::cout << "Daughters added to collection : " << daughterCollection.size() << std::endl;
      (*outputCollection).insert((*outputCollection).end(), daughterCollection.begin(), daughterCollection.end());
      bCode += 100;
    }
    if (genPart.status()!=1) continue;
    // Skip b daughters
    if (isFromB(genPart)) {
      std::cout << "Found from B!" << std::endl;
      fromBs++;
      continue;	
    }
    // Add the rest of the particles			 
    
    outputCollection->push_back(genPart);
  }
  std::cout << "particles fromB found (skipped) : " << fromBs << std::endl;
  std::cout << "Total particles in collection: " << (*outputCollection).size() << std::endl;
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
reco::GenParticleCollection CheatHFHadronReplacer::addDaughters(const reco::Candidate &particle, reco::GenParticleCollection collection, int bCode) const {

  reco::GenParticleCollection daughterCollection = collection;

  // Get daughters
  unsigned int ndaught = particle.numberOfDaughters();

  for (size_t i = 0; i < ndaught; i++) {
    
    const reco::Candidate& daughter = *particle.daughter(i);
    
    if (daughter.status() == 1) {
      std::cout << "Found a daughter!" << std::endl;
      daughterCollection.push_back(reco::GenParticle(daughter.charge(), daughter.p4(), daughter.vertex(), daughter.pdgId(), bCode, true));
      //daughterCollection.push_back(reco::GenParticle(daughter));
    }
    unsigned int ndaughtdaught = daughter.numberOfDaughters();
    if (ndaughtdaught > 0) {
      daughterCollection = addDaughters(daughter, daughterCollection, bCode);
    }
  }
  //std::cout << "Daughters found : " << daughterCollection.size() << std::endl;
  return daughterCollection;
}
			    

void CheatHFHadronReplacer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("genParticles", edm::InputTag("genParticles"));
  descriptions.add("CheatHFHadronReplacer", desc);
}


DEFINE_FWK_MODULE(CheatHFHadronReplacer);


