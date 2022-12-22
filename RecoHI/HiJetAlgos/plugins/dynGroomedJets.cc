
/*

  Perform dynamic grooming on jets
  Cannot be done in FastjetJetProducer, as this algorithm is not part of fastjet
  I borrowed liberally from CompoundJetProducer
  If I was a better programmer, I would just call that class

  -Matt Nguyen, 02-02-2022 (Groundhog's day) 

*/

#include <memory>
#include <vector>
#include <cmath>
#include <map>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/SoftDrop.hh"

#include "RecoJets/JetProducers/interface/JetSpecific.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/ShallowTagInfo.h"
#include "CommonTools/UtilAlgos/interface/DeltaR.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include <xgboost/c_api.h> 

using namespace edm;
using namespace std;
using namespace reco;
using namespace fastjet;

template <class T>
class dynGroomedJets : public global::EDProducer<> {
public:
  explicit dynGroomedJets(const ParameterSet&);
  ~dynGroomedJets() override = default;

  static void fillDescriptions(ConfigurationDescriptions&);

private:
  void produce(StreamID, Event&, const EventSetup&) const override;

  bool IterativeDeclustering(vector<PseudoJet> jetConstituents, PseudoJet *sub1, PseudoJet *sub2, vector<PseudoJet> &constit1, vector<PseudoJet> &constit2) const;
  BasicJet ConvertFJ2BasicJet(PseudoJet *fj, vector<PseudoJet>, Handle<PFCandidateCollection>, const EventSetup& iSetup) const;

  int trkInVector(reco::CandidatePtr trk, std::vector<reco::CandidatePtr> tracks) const;
  int trkGenPartMatch(reco::Jet::Constituent constituent, std::vector<pat::PackedGenParticle> genParticles, double ptCut) const;
  vector<PseudoJet> aggregateHF(const T& jet, float ptCut, std::vector<pat::PackedGenParticle> genParticles) const;
  vector<PseudoJet> aggregateHF(const T& jet, float ptCut, const std::string ipTagInfoLabel_, const std::string svTagInfoLabel_, const std::vector<pat::PackedGenParticle> genParticles) const;

  // ------------- member data ----------------------------
  EDGetTokenT<View<T> > jetSrc_;
  EDGetTokenT<PFCandidateCollection> constitSrc_;
  EDGetTokenT<edm::View<pat::PackedCandidate>>  packedConstitSrc_;
  const edm::EDGetTokenT<std::vector<reco::CandSecondaryVertexTagInfo>> svTagInfoSrc_;
  const edm::EDGetTokenT<std::vector<pat::PackedGenParticle>> genParticleSrc_;
  bool writeConstits_;
  bool doLateSD_;
  bool chargedOnly_;
  bool aggregateHF_;
  std::string ipTagInfoLabel_;
  std::string svTagInfoLabel_;
  double zcut_;
  double beta_;
  double dynktcut_;
  double rParam_;
  bool doGenJets_;
};

template <class T>
dynGroomedJets<T>::dynGroomedJets(const ParameterSet& iConfig)
  : jetSrc_(consumes<View<T> >(iConfig.getParameter<InputTag>("jetSrc"))),
    constitSrc_(mayConsume<PFCandidateCollection >(iConfig.getParameter<InputTag>("constitSrc"))),
    packedConstitSrc_(mayConsume<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("constitSrc"))),
    genParticleSrc_(mayConsume<std::vector<pat::PackedGenParticle>>(iConfig.getParameter<InputTag>("taggedGenParticles"))),
    writeConstits_(iConfig.getParameter<bool>("writeConstits")),
    doLateSD_(iConfig.getParameter<bool>("doLateSD")),
    chargedOnly_(iConfig.getParameter<bool>("chargedOnly")),
    aggregateHF_(iConfig.getParameter<bool>("aggregateHF")), // whatever is tagged in taggedGenParticles will be aggregated
    ipTagInfoLabel_(iConfig.getParameter<std::string>("ipTagInfoLabel")),
    svTagInfoLabel_(iConfig.getParameter<std::string>("svTagInfoLabel")),
    zcut_(iConfig.getParameter<double>("zcut")),
    beta_(iConfig.getParameter<double>("beta")),
    dynktcut_(iConfig.getParameter<double>("dynktcut")),
    rParam_(iConfig.getParameter<double>("rParam")),
    doGenJets_(iConfig.getParameter<bool>("doGenJets")) {
  string alias = (iConfig.getParameter<InputTag>("jetSrc")).label();
  produces<BasicJetCollection>().setBranchAlias(alias);
  produces<vector<BasicJet>>("SubJets").setBranchAlias(alias);
}

template <class T>
void dynGroomedJets<T>::produce(StreamID, Event& iEvent, const EventSetup& iSetup) const {
  // std::cout << "In dynGroomedJets" << std::endl;

  auto jetCollection = make_unique<reco::BasicJetCollection>();
  auto subjetCollection = make_unique<reco::BasicJetCollection>();

  // This will store the handle for the subjets after we write them
  OrphanHandle<reco::BasicJetCollection> subjetHandleAfterPut;
  // this is the mapping of subjet to hard jet
  vector<vector<int>> indices;
  // this is the list of hardjet 4-momenta
  vector<math::XYZTLorentzVector> p4_hardJets;
  // this is the hardjet areas
  vector<double> area_hardJets;
  vector<bool> isHardest;

  Handle<View<T> > jets;
  iEvent.getByToken(jetSrc_, jets);

  Handle<PFCandidateCollection> pfcands;
  bool isPF = iEvent.getByToken(constitSrc_, pfcands);

  Handle<edm::View<pat::PackedCandidate> > pfcandsPacked;
  bool isPackedPF = iEvent.getByToken(packedConstitSrc_, pfcandsPacked);

  // -- For aggregation -- //
  // Grab gen particles
  edm::Handle<std::vector<pat::PackedGenParticle>> genParticles;
  iEvent.getByToken(genParticleSrc_, genParticles);
  // std::cout << "nb of gen particles in collection: " << (*genParticles).size() << std::endl; 
  // --------------------- //

  indices.resize(jets->size());

  // ptCut to apply on tracks and gen particles
  float ptCut = 1.;

  int jetIndex = 0;
  for (const T& jet : *jets) { 
    // std::cout << "new jet with pt: " << jet.pt() << std::endl;
    const reco::GenJet *genJetTemp = jet.genJet();
    if (genJetTemp) 
      // std::cout << "\thas gen jet pt: " << genJetTemp->pt() << std::endl;

    p4_hardJets.push_back(math::XYZTLorentzVector(jet.px(), jet.py(), jet.pz(), jet.energy()));
    area_hardJets.push_back(jet.jetArea());

    std::vector<fastjet::PseudoJet> jetConstituents = {};

    fastjet::PseudoJet *subFJ1 = new fastjet::PseudoJet();
    fastjet::PseudoJet *subFJ2 = new fastjet::PseudoJet();
    std::vector<fastjet::PseudoJet> constit1;
    std::vector<fastjet::PseudoJet> constit2;

    // aggregate B based on jet type
    if (aggregateHF_) {
      if (doGenJets_) {
		    // std::cout << "\tAggregating HF for gen jet" << std::endl;
        const reco::GenJet *genJet = jet.genJet();
        if (genJet) jetConstituents = aggregateHF(*genJet, ptCut, *genParticles);
      } else {
		    // std::cout << "\tAggregating HF for reco jet" << std::endl;        
        jetConstituents = aggregateHF(jet, ptCut, ipTagInfoLabel_, svTagInfoLabel_, *genParticles);
      } 
    } else {
      // std::cout << "\tNot aggregating" << std::endl;
      std::vector<reco::CandidatePtr> constituents = {}; 
      if (doGenJets_) {
        const reco::GenJet *genJet = jet.genJet();
        if (genJet) constituents = genJet->getJetConstituents();
      } else {
        constituents = jet.getJetConstituents();
      }

      for (reco::CandidatePtr constituent : constituents) {
        if ((chargedOnly_) && (constituent->charge() == 0)) continue;
        if (constituent->pt() < ptCut) continue;
        jetConstituents.push_back(fastjet::PseudoJet(constituent->px(), constituent->py(), constituent->pz(), constituent->energy()));
      }
      // std::cout << "\tjet has " << jetConstituents.size() << " constituents" << std::endl;
    }
    //std::cout << "nb of jet constituents: " << jetConstituents.size() << std::endl;
    //for (auto constituent : jetConstituents) {
    //std::cout << "eta: " << constituent.eta() << std::endl;
    //std::cout << "phi: " << constituent.phi() << std::endl;
    //}


    // Iterative declustering
    isHardest.push_back(IterativeDeclustering(jetConstituents, subFJ1, subFJ2, constit1, constit2));

    //std::cout << "sjt1 pt: " << subFJ1->pt() << std::endl;
    //std::cout << "sjt2 pt: " << subFJ2->pt() << std::endl;

    // Convert fastjets to basicjets    
    BasicJet subjet1 = ConvertFJ2BasicJet(subFJ1, constit1, pfcands, iSetup);
    BasicJet subjet2 = ConvertFJ2BasicJet(subFJ2, constit2, pfcands, iSetup);

    if(subjet1.pt() > 1.0e-3){
      indices[jetIndex].push_back(subjetCollection->size());
      if (subFJ1->has_area()) subjet1.setJetArea(subFJ1->area());
      subjetCollection->push_back(subjet1);
    }
    if(subjet2.pt() > 1.0e-3){
      indices[jetIndex].push_back(subjetCollection->size());
      if(subFJ2->has_area())subjet2.setJetArea(subFJ2->area());
      subjetCollection->push_back(subjet2);
    }

    jetIndex++;
  }

  // put subjets into event record
  subjetHandleAfterPut = iEvent.put(move(subjetCollection), "SubJets"); 

  // Now create the hard jets with ptr's to the subjets as constituents
  vector<math::XYZTLorentzVector>::const_iterator ip4 = p4_hardJets.begin(), ip4Begin = p4_hardJets.begin(), ip4End = p4_hardJets.end();
  
  for (; ip4 != ip4End; ++ip4) {
    int p4_index = ip4 - ip4Begin;
    vector<int>& ind = indices[p4_index];
    vector<CandidatePtr> i_hardJetConstituents;

    // Add the subjets to the hard jet
    for (vector<int>::const_iterator isub = ind.begin(); isub != ind.end(); ++isub) {
      CandidatePtr candPtr(subjetHandleAfterPut, *isub, false);
      i_hardJetConstituents.push_back(candPtr);
    }
    
    Particle::Point point(0, 0, 0);
    //cout<<" size of i_hardJetConstituents "<<i_hardJetConstituents.size()<<endl;
    BasicJet toput(*ip4, point, i_hardJetConstituents);
    if(isHardest[ip4-ip4Begin]) toput.setJetArea(0.8);
    else toput.setJetArea(0.4);
    //toput.setIsWeighted(isHardest[ip4 - ip4Begin]);  // Not yet in CMSSW_9_X_X
    //cout << "hardjet added to collection" << endl;
    jetCollection->push_back(toput);
  }

  //cout << "jetCollection size: " << jetCollection->size() << endl; 
  iEvent.put(move(jetCollection));
}

template <class T>
BasicJet dynGroomedJets<T>::ConvertFJ2BasicJet(fastjet::PseudoJet *fj, vector<fastjet::PseudoJet> constit, Handle<PFCandidateCollection> pfcands, const EventSetup& iSetup) const
{

  math::XYZTLorentzVector p4(fj->px(), fj->py(), fj->pz(), fj->e());  
  Particle::Point point(0, 0, 0);
  vector<CandidatePtr> constituents;
  if(writeConstits_){
    for(uint j=0;j<constit.size();j++){
      double constitPt = constit[j].pt();
      double constitEta = constit[j].eta();
      
      bool foundMatch = false;
      int iCand = -1;
      if (foundMatch) {} // so the compiler does not complain that it is not used
      for (const PFCandidate& pfcand : *pfcands) {
	iCand++;
	
	if(fabs(constitPt - pfcand.pt()) < 0.0001 && fabs(constitEta - pfcand.eta()) <0.0001 ){
	  foundMatch = true;
	  constituents.push_back(CandidatePtr(pfcands, iCand));
	  break;
	}
      }
      //if(!foundMatch) cout<<" couldn't find match "<<endl;
    }
  }

  BasicJet basicjet;
  writeSpecific(basicjet, p4, point, constituents, iSetup);

  return basicjet;

}

/*
template <class T>
bool dynGroomedJets<T>::IterativeDeclusteringReco(const T& jet, std::vector<reco::CandidatePtr>& ipTracks, std::vector<reco::CandidatePtr>& svTracks, float ptCut, fastjet::PseudoJet *sub1, fastjet::PseudoJet *sub2, std::vector<fastjet::PseudoJet> &constit1, std::vector<fastjet::PseudoJet> &constit2) const
{
  
  //  Iterative declustering with B aggregation 
  //  using SV information for reco jets
  
  bool flagSubjet=false;
  bool isHardest=false;
  double kt1=-1;
  double nsplit=0;
  double nsel=0;
  double nsdin=-1;
  double jet_radius_ca = 1.0;
  JetDefinition jet_def(genkt_algorithm,jet_radius_ca,0,static_cast<RecombinationScheme>(0), Best);

  // Reclustering jet constituents with new algorithm                                                                                          
  try
    {      
      // Input particle collection 
      std::vector<reco::CandidatePtr> daughters = jet.getJetConstituents();

      // Particle collection to use in the declustering
      std::vector<fastjet::PseudoJet> particles;
      int nCharged = 0;

      // Particle collection to aggregate into a pseudo-B 
      std::vector<reco::CandidatePtr> hfProducts;
      
      
      //cout << "nb of tracks in jet: " << ipTracks.size() << std::endl;
      //cout << "nb of secondary vertices in jet: " << svTagInfo.nVertices() << std::endl;
      //cout << "nb of tracks associated to SV's: " << svTagInfo.nVertexTracks() << std::endl;
      //cout << "nb of tracks associated to 1st sv: " << svTagInfo.nVertexTracks(0) << std::endl;
      //cout << "nb of tracks associated to 2nd sv: " << svTagInfo.nVertexTracks(1) << std::endl;

      //cout << "nb of tracks associated to SV's ? : "<< svTagInfo.vertexTracks().size() << std::endl;
      
      
      //cout << "new jet" << endl;

      for(reco::CandidatePtr daughter : daughters) {
        // Charge and pt cut
        if(chargedOnly_ && daughter->charge() == 0) continue;
        if(daughter->pt() < ptCut) continue;

        //float eps = 0.001;

        //cout << "new daughter" << endl;

        // Look for track in ipTracks, otherwise discard
        bool foundTrack = trkInVector(daughter, ipTracks);
        if(!foundTrack) continue;

        // Look for track in svTracks
        bool foundTrackInSV = trkInVector(daughter, svTracks);
        if (foundTrackInSV) { 
          hfProducts.push_back(daughter);
        } else {
          particles.push_back(PseudoJet(daughter->px(), daughter->py(), daughter->pz(), daughter->energy()));
        }
        nCharged++;
      }

      //cout << "all non B decay products added to particles" << endl;
      //cout << "tracks associated to any sv : " << hfProducts.size() << endl;

      // Aggregate all hfProducts into a single pseudo-B
      if (hfProducts.size() > 0) {
        double px = 0;
        double py = 0;
        double pz = 0;
        double energy = 0;
        for (reco::CandidatePtr bProduct : hfProducts) {
          px += bProduct->px();
          py += bProduct->py();
          pz += bProduct->pz();
          energy += bProduct->energy();
        }
        //cout << "Added B to particles with " << hfProducts.size() << " constituents" << endl;
        particles.push_back(PseudoJet(px, py, pz, energy));
        nCharged++;
      }

      if(nCharged == 0)	return false;

      ClusterSequence csiter(particles, jet_def);
      vector<PseudoJet> output_jets = csiter.inclusive_jets(0);
      output_jets = sorted_by_pt(output_jets);
      
      
      PseudoJet jj = output_jets[0];
      PseudoJet j1;
      PseudoJet j2;                                                                                                               
      PseudoJet j1first;
      PseudoJet j2first;
      
      while(jj.has_parents(j1,j2)) {
        if(j1.perp() < j2.perp()) swap(j1,j2);
        
        double delta_R = j1.delta_R(j2);
        double cut=zcut_*pow(delta_R/rParam_,beta_);
        double z=j2.perp()/(j1.perp()+j2.perp());
        if(z > cut && (doLateSD_ || !flagSubjet) ){
          flagSubjet=true;
          j1first =j1;
          j2first =j2;
          *sub1 = j1first;
          *sub2 = j2first;
          nsdin=nsplit;
        }
        double dyn=z*(1-zreco::GenJet)*j2.perp()*pow(delta_R/rParam_,dynktcut_);
        
        if(dyn>kt1){
          nsel=nsplit;
          kt1=dyn;
        }
        nsplit=nsplit+1;
        jj=j1;
      }
          
      if(!flagSubjet) *sub1 = output_jets[0];

      if(sub1->has_constituents()) constit1 = sub1->constituents(); 
      if(sub2->has_constituents()) constit2 = sub2->constituents(); 
      if(nsel==nsdin)isHardest = true; 
    } catch (fastjet::Error) {} // { return -1; }

  //cout << "sub1 pt: " << sub1->pt() << endl;
  //cout << "sub2 pt: " << sub2->pt() << endl;
  return isHardest;
}
*/

/*
    // IterativeDeclustering prototype

template <class T>
bool dynGroomedJets<T>::IterativeDeclustering(const T& jet, PseudoJet *sub1, PseudoJet *sub2, vector<PseudoJet> &constit1, vector<PseudoJet> &constit2) const
{


  bool flagSubjet=false;
  bool isHardest=false;
  double kt1=-1;
  double nsplit=0;
  double nsel=0;
  double nsdin=-1;
  double jet_radius_ca = 1.0;
  JetDefinition jet_def(genkt_algorithm,jet_radius_ca,0,static_cast<RecombinationScheme>(0), Best);

  // Reclustering jet constituents with new algorithm                                                                                          
  
  try
    {
      vector<PseudoJet> particles;

      int nCharged = 0;
      auto daughters = jet.getJetConstituents();
      //std::cout<<" nConst "<<daughters.size()<<std::endl;
	
      for(auto it = daughters.begin(); it!=daughters.end(); ++it){
	if(chargedOnly_ && (**it).charge() == 0) continue;
	if((**it).pt()<1.0) continue;

	particles.push_back(PseudoJet((**it).px(), (**it).py(), (**it).pz(), (**it).energy()));
	nCharged++;
      }
      //std::cout<<" nCharged =  "<<nCharged<<std::endl;
      if(nCharged == 0) return false;
      ClusterSequence csiter(particles, jet_def);
      vector<PseudoJet> output_jets = csiter.inclusive_jets(0);
      output_jets = sorted_by_pt(output_jets);
      
      
      PseudoJet jj = output_jets[0];
      PseudoJet j1;
      PseudoJet j2;                                                                                                               
      PseudoJet j1first;
      PseudoJet j2first;
      
      while(jj.has_parents(j1,j2))
	{
	  
	  if(j1.perp() < j2.perp()) swap(j1,j2);
	  
	  double delta_R = j1.delta_R(j2);
	  double cut=zcut_*pow(delta_R/rParam_,beta_);
	  double z=j2.perp()/(j1.perp()+j2.perp());
	  if(z > cut && (doLateSD_ || !flagSubjet) ){
	    flagSubjet=true;
	    j1first =j1;
	    j2first =j2;
	    *sub1 = j1first;
	    *sub2 = j2first;
	    nsdin=nsplit;
	  }
	  double dyn=z*(1-z)*j2.perp()*pow(delta_R/rParam_,dynktcut_);
	  
	  if(dyn>kt1){
	    nsel=nsplit;
	    kt1=dyn;
	  }
	  nsplit=nsplit+1;
	  jj=j1;

	  
	}
      
      if(!flagSubjet) *sub1 = output_jets[0];

      if(sub1->has_constituents()) constit1 = sub1->constituents(); 
      if(sub2->has_constituents()) constit2 = sub2->constituents(); 
      if(nsel==nsdin)isHardest = true; 
    } catch (fastjet::Error) { }//return -1; }

  return isHardest;
}
*/

/*
template <class T>
bool dynGroomedJets<T>::IterativeDeclusteringGen(const T& jet, std::vector<pat::PackedGenParticle> genParticles, float ptCut, fastjet::PseudoJet *sub1, fastjet::PseudoJet *sub2, std::vector<fastjet::PseudoJet>& constit1, std::vector<fastjet::PseudoJet>& constit2) const
{
   
  //  Iterative declustering with aggregation
  //  of B hadron based on gen particle matching
  
  bool flagSubjet = false;
  bool isHardest = false;
  double kt1 = -1;
  double nsplit = 0;
  double nsel = 0;
  double nsdin = -1;
  double jet_radius_ca = 1.0;
  JetDefinition jet_def(genkt_algorithm, jet_radius_ca, 0, static_cast<RecombinationScheme>(0), Best);

  try {
    // Input and output particle collections
    std::vector<reco::CandidatePtr> inputConstituents = jet.getJetConstituents();
    std::vector<fastjet::PseudoJet> outputConstituents;
    std::map<int, std::vector<reco::CandidatePtr>> hfConstituentsMap;

    // Keep track of charged particles
    int nCharged = 0;
  

    for (reco::CandidatePtr constituent : inputConstituents) {
      // Charge and pt cut
      if(chargedOnly_ && constituent->charge() == 0) continue;
      if(constituent->pt() < ptCut) continue;

      // Match with gen particle -- should be exact match
      int genStatus = trkGenPartMatch(constituent, genParticles, ptCut);
    
      // Collect particles coming from B decays, add the rest directly to the collection
      if (genStatus >= 100) {
	      hfConstituentsMap[genStatus].push_back(constituent);
      } else if (genStatus == 1) {
        nCharged++;
        outputConstituents.push_back(PseudoJet(constituent->px(), constituent->py(), constituent->pz(), constituent->energy()));
      }
    }

    // Aggregate particles coming from B decays into pseudo-B's 
    // and add them to the collection
    for (auto it = hfConstituentsMap.begin(); it != hfConstituentsMap.end(); ++it) {
      double px = 0;
      double py = 0;
      double pz = 0;
      double energy = 0;
      for (reco::CandidatePtr hfConstituent : it->second) {
        px += hfConstituent->px();
        py += hfConstituent->py();
        pz += hfConstituent->pz();
        energy += hfConstituent->energy();
      }
      //cout << "Added B to particles with " << (it->second).size() << "constituents" << endl;
      nCharged++;
      outputConstituents.push_back(PseudoJet(px, py, pz, energy));
    }

    //cout << "nCharged = " << nCharged << endl;
    if (nCharged == 0) return false;

    // Perform the declustering 
    ClusterSequence csiter(outputConstituents, jet_def);
    vector<PseudoJet> output_jets = csiter.inclusive_jets(0);
    output_jets = sorted_by_pt(output_jets);
    
    PseudoJet jj = output_jets[0];
    PseudoJet j1;
    PseudoJet j2;                                                                                                             
    PseudoJet j1first;
    PseudoJet j2first;
   
    while(jj.has_parents(j1, j2)) {
      if(j1.perp() < j2.perp()) swap(j1,j2);
	  
      double delta_R = j1.delta_R(j2);
      double cut = zcut_ * pow(delta_R / rParam_,beta_);
      double z = j2.perp() / (j1.perp() + j2.perp());
      if(z > cut && (doLateSD_ || !flagSubjet) ) {
        flagSubjet=true;
        j1first = j1;
        j2first = j2;
        *sub1 = j1first;
        *sub2 = j2first;
        nsdin = nsplit;
      }
      double dyn = z * (1-z) * j2.perp() * pow(delta_R / rParam_, dynktcut_);
	 
      if(dyn > kt1) {
        nsel = nsplit;
        kt1 = dyn;
        }
        nsplit = nsplit + 1;
        jj = j1;
    } 

    if(!flagSubjet) *sub1 = output_jets[0];

    if (sub1->has_constituents()) constit1 = sub1->constituents(); 
    if (sub2->has_constituents()) constit2 = sub2->constituents(); 
    if (nsel == nsdin) isHardest = true; 
  } catch (fastjet::Error) {} // {return -1;}
  //cout << "sub1 pt: " << sub1->pt() << endl;
  //cout << "sub2 pt: " << sub2->pt() << endl;
  return isHardest;
}
*/

template <class T>
bool dynGroomedJets<T>::IterativeDeclustering(vector<PseudoJet> jetConstituents, 
                                              PseudoJet *sub1, PseudoJet *sub2, 
                                              vector<PseudoJet> &constit1, vector<PseudoJet> &constit2) const
{
   
  // Iterative declustering for any type of jet
  // given its constituents

  bool flagSubjet=false;
  bool isHardest=false;
  double kt1=-1;
  double nsplit=0;
  double nsel=0;
  double nsdin=-1;
  double jet_radius_ca = 1.0;
  JetDefinition jet_def(genkt_algorithm,jet_radius_ca,0,static_cast<RecombinationScheme>(0), Best);

  // Return false if no constituents
  if (jetConstituents.size() == 0) return false;

  // Reclustering jet constituents with new algorithm                                                                                          
  try
    {
      ClusterSequence csiter(jetConstituents, jet_def);
      vector<PseudoJet> output_jets = csiter.inclusive_jets(0);
      output_jets = sorted_by_pt(output_jets);
      
      
      PseudoJet jj = output_jets[0];
      PseudoJet j1;
      PseudoJet j2;                                                                                                               
      PseudoJet j1first;
      PseudoJet j2first;
      
      while(jj.has_parents(j1,j2)) {
	  
        if(j1.perp() < j2.perp()) swap(j1,j2);
        
        double delta_R = j1.delta_R(j2);
        double cut=zcut_*pow(delta_R/rParam_,beta_);
        double z=j2.perp()/(j1.perp()+j2.perp());
        if(z > cut && (doLateSD_ || !flagSubjet) ) {
          flagSubjet=true;
          j1first =j1;
          j2first =j2;
          *sub1 = j1first;
          *sub2 = j2first;
          nsdin=nsplit;
        }
        double dyn=z*(1-z)*j2.perp()*pow(delta_R/rParam_,dynktcut_);
        
        if(dyn>kt1){
          nsel=nsplit;
          kt1=dyn;
        }
        nsplit=nsplit+1;
        jj=j1;
      }
      
      if(!flagSubjet) *sub1 = output_jets[0];

      if (sub1->has_constituents()) constit1 = sub1->constituents(); 
      if (sub2->has_constituents()) constit2 = sub2->constituents(); 
      if (nsel==nsdin) isHardest = true; 
    } catch (fastjet::Error) { }//return -1; }
  return isHardest;
}


template <class T>
void dynGroomedJets<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setComment("Dynamically groomed jets");
  desc.add<bool>("writeConstits", false);
  desc.add<bool>("doLateSD", false);
  desc.add<bool>("chargedOnly", false);
  desc.add<bool>("aggregateHF", false);
  desc.add<std::string>("ipTagInfoLabel", "pfImpactParameter");
  desc.add<std::string>("svTagInfoLabel", "pfInclusiveSecondaryVertexFinder");
  desc.add<edm::InputTag>("taggedGenParticles", edm::InputTag("HFdecayProductTagger"));
  desc.add<double>("zcut", 0.1);
  desc.add<double>("beta", 0.0);
  desc.add<double>("dynktcut", 1.0);
  desc.add<double>("rParam", 0.4);
  desc.add<bool>("doGenJets", false);
  desc.add<edm::InputTag>("jetSrc", edm::InputTag("slimmedJets"));
  desc.add<edm::InputTag>("constitSrc", edm::InputTag("packedPFCandidates"));
  descriptions.add("dynGroomedPATJets", desc);
  // if (typeid(T) == typeid(pat::Jet)) {
  //   desc.add<edm::InputTag>("jetSrc", edm::InputTag("slimmedJets"));
  //   desc.add<edm::InputTag>("constitSrc", edm::InputTag("packedPFCandidates"));
  //   descriptions.add("dynGroomedPFJets", desc);
  // }
  // else if (typeid(T) == typeid(reco::GenJet)) {
  //   desc.add<edm::InputTag>("jetSrc", edm::InputTag("slimmedGenJets"));
  //   desc.add<edm::InputTag>("constitSrc", edm::InputTag("packedGenParticles"));
  //   descriptions.add("dynGroomedGenJets", desc);
  // }

}

// Function to match reconstructed tracks to gen particles geometrically
template <class T>
int dynGroomedJets<T>::trkGenPartMatch(reco::Jet::Constituent constituent, std::vector<pat::PackedGenParticle> genParticles, double ptCut) const
{
  // return status of matched particle
  int status = 1;
  
  double dRmin = 100.;
  double dR = dRmin;

  double trkPt = constituent->pt();
  double trkEta = constituent->eta();
  double trkPhi = constituent->phi();

  for (pat::PackedGenParticle genParticle : genParticles) {
    double partPt = genParticle.pt();
    double relPt = trkPt / partPt;
    //cout << "genParticle status: " << genParticle.status() << endl;
    if (partPt < 0.5 * ptCut) continue;
	  if(chargedOnly_ && genParticle.charge() == 0) continue;
	
    double partEta = genParticle.eta();
    double partPhi = genParticle.phi();

    double dEta = trkEta - partEta;
    double dPhi = std::acos(std::cos(trkPhi - partPhi));

    dR = std::sqrt((dEta * dEta) + (dPhi * dPhi));

    // edit for conservative matching:
    if (dR > 0.02) continue;
    if ((relPt < 0.8) || (relPt > 1.2)) continue;

    if (dR < dRmin) {
      dRmin = dR;
      status = genParticle.status();
    }
  }
  //cout << "dRmin: " << dRmin << endl;
  return status;  
}

template <class T>
int dynGroomedJets<T>::trkInVector(reco::CandidatePtr trk, std::vector<reco::CandidatePtr> tracks) const {
  int indexInVector = -1;
  auto it = find(tracks.begin(), tracks.end(), trk);
  if (it != tracks.end()) {
	indexInVector = it - tracks.begin();
  }
  return indexInVector;
}


// Function to aggregate reconstructed tracks into pseudo-B's 
template <class T>
std::vector<fastjet::PseudoJet> dynGroomedJets<T>::aggregateHF(const T& jet, float ptCut,
                                                               const std::string ipTagInfoLabel_, 
                                                               const std::string svTagInfoLabel_,
                                                               const std::vector<pat::PackedGenParticle> genParticles) const
{
  // std::cout << "Aggregating Bs in reco jet" << std::endl;
  // Input particle collection 
  std::vector<reco::CandidatePtr> inputJetConstituents = jet.getJetConstituents();
  // std::cout << "Got the jet constituents" << std::endl;
  // Grab the IP and SV tag info 
  const reco::CandIPTagInfo *ipTagInfo = jet.tagInfoCandIP(ipTagInfoLabel_.c_str());
  const reco::CandSecondaryVertexTagInfo *svTagInfo = jet.tagInfoCandSecondaryVertex(svTagInfoLabel_.c_str());

  bool withBDT = false;
  if (withBDT) {
    // Load the BDT 
    std::string model_path = "~/gitRepos/jetAnalysis/bdt/saved_models/xgboost_bdt_cmssw.model";
    // std::vector<std::string> features = {"f0",  "f1",  "f2",  "f3"};

    std::cout << "Loading the trained BDT from " << model_path << std::endl;
    // Declare Object
    BoosterHandle booster_;
    // Allocate memory in C style
    XGBoosterCreate(NULL, 0, &booster_);
    // Load Model
    XGBoosterLoadModel(booster_, model_path.c_str()); // second argument should be a const char *.

    const int nFeatures = 4;
    const int nEntries = 1;
    const float missing_value = -1000000.;
    float TestData[nEntries][nFeatures] = {{1, 2., 3., 4.}}; 
    // Suppose nEntries data points, each data point has nFeatures dimension
    
    // Assign data to the "TestData" 2d array ... 
    // Declare object
    DMatrixHandle data_;
    // Allocate memory and use external float array to initialize
    XGDMatrixCreateFromMat((float *)TestData, nEntries, nFeatures, missing_value, &data_); 
    // The first argument takes in float * namely 1d float array only, 
    // 2nd & 3rd: shape of input, 4th: value to replace missing ones

    bst_ulong out_len; // bst_ulong is a typedef of unsigned long
    const float *f; // array to store predictions
    XGBoosterPredict(booster_, data_, 0, 0, &out_len, &f);// lower version API
    
  }

  // Particle collection to use in the declustering
  std::vector<fastjet::PseudoJet> outputJetConstituents = {};

  // Particle collection to aggregate into a pseudo-B 
  //std::vector<reco::CandidatePtr> hfProducts;
  std::map<int, std::vector<reco::CandidatePtr>> hfConstituentsMap;

  // Get ipTracks and svTracks for jet
  std::vector<reco::CandidatePtr> ipTracks = ipTagInfo->selectedTracks();
  std::vector<reco::CandidatePtr> svTracks = svTagInfo->vertexTracks();

  for(reco::CandidatePtr constituent : inputJetConstituents) {
    // Charge and pt cut
    if(chargedOnly_ && constituent->charge() == 0) continue;
    if(constituent->pt() < ptCut) continue;

    // Look for track in ipTracks, otherwise discard
    int whichTrack = trkInVector(constituent, ipTracks);
    if (whichTrack < 0) continue;
    // std::cout << "New track!" << std::endl;

    // Look for track in svTracks
    // int whichTrackInSV = trkInVector(constituent, svTracks);
    // bool inSV = (whichTrackInSV >= 0);

    // if (inSV) {
    //   Measurement1D m1D = svTagInfo.flightDistance(0);
    //   std::cout << "svtxdls of sv0: " << m1D.significance() << std::endl;
    // }

    // get ip3dSig
    // reco::btag::TrackIPData IPdata = ipTagInfo.impactParameterData()[whichTrack];
    // Float_t ip3dSig = IPdata.ip3d.significance();
    // std::cout << "constituent pt: " << constituent->pt() << ", ip data pt: " << ipTracks[whichTrack]->pt() << std::endl;
    // std::cout << "new constituent:" << std::endl;
    // std::cout << "ip3dSig = " << ip3dSig << std::endl;
    // std::cout << "is in SV = " << inSV << std::endl;

    // bool isHFproduct = (inSV && (std::abs(ip3dSig) > 3.)) || (std::abs(ip3dSig) > 9.);
    // ---- End Aggregating with BDT ---- //
    

	  // --- Aggregating with truth info ---- //
    int genStatus = trkGenPartMatch(constituent, genParticles, ptCut);
    // std::cout << "Truth status: " << genStatus << std::endl;
    
    // Collect particles coming from HF decays, add the rest directly to the collection
    if (genStatus >= 100) {
      hfConstituentsMap[genStatus].push_back(constituent);
	  	// std::cout << "Found HF constituent" << std::endl;
    } else if (genStatus == 1) {
      outputJetConstituents.push_back(PseudoJet(constituent->px(), constituent->py(), constituent->pz(), constituent->energy()));
		  // std::cout << "Not a B constituent, added to collection" << std::endl;
    }
  } // end loop over input jet constituents

  // Aggregate particles coming from HF decays into pseudo-B/C's 
  // and add them to the collection
  for (auto it = hfConstituentsMap.begin(); it != hfConstituentsMap.end(); ++it) {
	  reco::Candidate::PolarLorentzVector pseudoHF(0., 0., 0., 0.);
    for (reco::CandidatePtr hfConstituent : it->second) {
      reco::Candidate::PolarLorentzVector productLorentzVector(0., 0., 0., 0.);
      productLorentzVector.SetPt(hfConstituent->pt());
      productLorentzVector.SetEta(hfConstituent->eta());
      productLorentzVector.SetPhi(hfConstituent->phi());
      productLorentzVector.SetM(hfConstituent->mass());
      pseudoHF += productLorentzVector;
	  }
    outputJetConstituents.push_back(PseudoJet(pseudoHF.Px(), pseudoHF.Py(), pseudoHF.Pz(), pseudoHF.E()));
  } // end map loop

	//std::cout << "daughter pt: " << constituent->pt() << std::endl;
	//std::cout << "daughter inSV: " << foundTrackInSV << std::endl;
	//std::cout << "daughter ip3dSig: " << ip3dSig << std::endl;

  //if (isHFproduct) { 
    //std::cout << "is b product" << std::endl;
  //  hfProducts.push_back(constituent);
	//} else {
	    //std::cout << "is not b product, added to outputJetConstituents" << std::endl;
  //    outputJetConstituents.push_back(PseudoJet(constituent->px(), constituent->py(), constituent->pz(), constituent->energy()));
	//}
  

  // Aggregate all hfProducts into a single pseudo-B
  //if (hfProducts.size() > 0) {
  //  reco::Candidate::PolarLorentzVector pseudoHF(0., 0., 0., 0.);
  //    for (reco::CandidatePtr product : hfProducts) {
  //      reco::Candidate::PolarLorentzVector productLorentzVector(0., 0., 0., 0.);
  //      productLorentzVector.SetPt(product->pt());
  //      productLorentzVector.SetEta(product->eta());
  //      productLorentzVector.SetPhi(product->phi());
  //      productLorentzVector.SetM(product->mass());
  //      pseudoHF += productLorentzVector;
  //      //std::cout << "b product lorentz vector: " << productLorentzVector << std::endl;
  //  }
	  //std::cout << "pseudo-B added to consituents" << std::endl;
  //  outputJetConstituents.push_back(PseudoJet(pseudoHF.Px(), pseudoHF.Py(), pseudoHF.Pz(), pseudoHF.E()));
  //}
  return outputJetConstituents;
} // end aggregateHF()

// Function to aggregate gen particles coming from B decays into pseudo-B's 
template <class T>
std::vector<fastjet::PseudoJet> dynGroomedJets<T>::aggregateHF(const T& jet, float ptCut,
                                                              std::vector<pat::PackedGenParticle> genParticles) const
{
    // std::cout << "Aggregating HF at gen level" << std::endl;
  
    // Input and output particle collections
    std::vector<reco::CandidatePtr> inputConstituents = jet.getJetConstituents();
    std::vector<fastjet::PseudoJet> outputConstituents;
    std::map<int, std::vector<reco::CandidatePtr>> hfConstituentsMap;

    for (reco::CandidatePtr constituent : inputConstituents) {
      // Charge and pt cut
      if(chargedOnly_ && constituent->charge() == 0) continue;
      if(constituent->pt() < ptCut) continue;
	  
      // Match with gen particle -- should be exact match
      int genStatus = trkGenPartMatch(constituent, genParticles, ptCut);
	  //std::cout << "New jet constituent with pt = " << constituent->pt() << " and status = " << genStatus << std::endl;

      // Collect particles coming from HF decays, add the rest directly to the collection
      if (genStatus >= 100) {
	      hfConstituentsMap[genStatus].push_back(constituent);
		    // std::cout << "Found HF constituent" << std::endl;
      } else if (genStatus == 1) {
        outputConstituents.push_back(PseudoJet(constituent->px(), constituent->py(), constituent->pz(), constituent->energy()));
		    // std::cout << "Not a HF constituent, added to collection" << std::endl;
      }
    }

    // Aggregate particles coming from HF decays into pseudo-B/D's 
    // and add them to the collection
    for (auto it = hfConstituentsMap.begin(); it != hfConstituentsMap.end(); ++it) {
        reco::Candidate::PolarLorentzVector pseudoHF(0., 0., 0., 0.);
        for (reco::CandidatePtr hfConstituent : it->second) {
          reco::Candidate::PolarLorentzVector productLorentzVector(0., 0., 0., 0.);
          productLorentzVector.SetPt(hfConstituent->pt());
          productLorentzVector.SetEta(hfConstituent->eta());
          productLorentzVector.SetPhi(hfConstituent->phi());
          productLorentzVector.SetM(hfConstituent->mass());
          pseudoHF += productLorentzVector;
          //std::cout << "hf product lorentz vector: " << productLorentzVector << std::endl;
      }
      // std::cout << "added HF to collection" << std::endl;
      outputConstituents.push_back(PseudoJet(pseudoHF.Px(), pseudoHF.Py(), pseudoHF.Pz(), pseudoHF.E()));
    }
  return outputConstituents;
}
  
// using dynGroomedGenJets = dynGroomedJets<reco::GenJet>;
// using dynGroomedPFJets = dynGroomedJets<pat::Jet>;
using dynGroomedPATJets = dynGroomedJets<pat::Jet>;

// define this as a plug-in
// DEFINE_FWK_MODULE(dynGroomedPFJets);
// DEFINE_FWK_MODULE(dynGroomedGenJets);
DEFINE_FWK_MODULE(dynGroomedPATJets);
