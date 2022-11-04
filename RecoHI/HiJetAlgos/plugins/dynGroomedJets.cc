
/*

  Perform dynamic grooming on jets
  Cannot be done in FastjetJetProducer, as this algorithm is not part of fastjet
  I borrowed liberally from CompoundJetProducer
  If I was a better programmer, I would just call that class

  -Matt Nguyen, 02-02-2022 (Groundhog's day) 

*/

#include <dirent.h> // DEBUG

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

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"
#include <xgboost/c_api.h> 

template <class T>
class dynGroomedJets : public edm::global::EDProducer<> {
public:
  explicit dynGroomedJets(const edm::ParameterSet&);
  ~dynGroomedJets() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  bool IterativeDeclustering(std::vector<fastjet::PseudoJet>, 
                             fastjet::PseudoJet *, fastjet::PseudoJet *, 
                             std::vector<fastjet::PseudoJet> &, std::vector<fastjet::PseudoJet> &) const;

  reco::BasicJet ConvertFJ2BasicJet(fastjet::PseudoJet *, std::vector<fastjet::PseudoJet>, 
                              edm::Handle<reco::PFCandidateCollection>, const edm::EventSetup&) const;

  int trkInVector(reco::CandidatePtr, std::vector<reco::CandidatePtr>) const;
  int trkGenPartMatch(reco::Jet::Constituent, reco::GenParticleCollection, double) const;

  std::vector<fastjet::PseudoJet> aggregateHFgen(const T&, float, reco::GenParticleCollection) const;
  std::vector<fastjet::PseudoJet> aggregateHFreco(const T&, float, const reco::CandIPTagInfo, 
                                                  const reco::CandSecondaryVertexTagInfo&, 
                                                  const reco::GenParticleCollection, 
                                                  BoosterHandle) const;

  // ------------- member data ----------------------------
  edm::EDGetTokenT<edm::View<T> > jetSrc_;
  edm::EDGetTokenT<reco::PFCandidateCollection> constitSrc_;
  edm::EDGetTokenT<edm::View<reco::BaseTagInfo>> tagInfoSrc_;
  const edm::EDGetTokenT<std::vector<reco::CandSecondaryVertexTagInfo>> svTagInfos_;
  const edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticleSrc_;
  bool writeConstits_;
  bool doLateSD_;
  bool chargedOnly_;
  bool aggregateHF_;
  edm::FileInPath model_path;
  double zcut_;
  double beta_;
  double dynktcut_;
  double rParam_;
};

template <class T>
dynGroomedJets<T>::dynGroomedJets(const edm::ParameterSet& iConfig)
  : jetSrc_(consumes<edm::View<T> >(iConfig.getParameter<edm::InputTag>("jetSrc"))),
    constitSrc_(mayConsume<reco::PFCandidateCollection >(iConfig.getParameter<edm::InputTag>("constitSrc"))),
    tagInfoSrc_(mayConsume<edm::View<reco::BaseTagInfo>>(iConfig.getParameter<edm::InputTag>("ak4PFPfImpactParameterTagInfos"))),
    svTagInfos_(mayConsume<std::vector<reco::CandSecondaryVertexTagInfo>>(iConfig.getParameter<edm::InputTag>("ak4PFPfInclusiveSecondaryVertexFinderTagInfos"))),
    genParticleSrc_(mayConsume<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("HFdecayProductTagger"))),
    writeConstits_(iConfig.getParameter<bool>("writeConstits")),
    doLateSD_(iConfig.getParameter<bool>("doLateSD")),
    chargedOnly_(iConfig.getParameter<bool>("chargedOnly")),
    aggregateHF_(iConfig.getParameter<bool>("aggregateHF")), // whatever is tagged in HFdecayProductTagger will be aggregated at truth
    model_path(iConfig.getParameter<edm::FileInPath>("model_path")),
    zcut_(iConfig.getParameter<double>("zcut")),
    beta_(iConfig.getParameter<double>("beta")),
    dynktcut_(iConfig.getParameter<double>("dynktcut")),
    rParam_(iConfig.getParameter<double>("rParam")){
  std::string alias = (iConfig.getParameter<edm::InputTag>("jetSrc")).label();
  produces<reco::BasicJetCollection>().setBranchAlias(alias);
  produces<std::vector<reco::BasicJet>>("SubJets").setBranchAlias(alias);
}

template <class T>
void dynGroomedJets<T>::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  //std::cout << "In dynGroomedJets::produce" << std::endl;

  // ptCut to apply on tracks and gen particles
  float ptCut = 1.;

  auto jetCollection = std::make_unique<reco::BasicJetCollection>();
  auto subjetCollection = std::make_unique<reco::BasicJetCollection>();

  // This will store the handle for the subjets after we write them
  edm::OrphanHandle<reco::BasicJetCollection> subjetHandleAfterPut;
  // this is the mapping of subjet to hard jet
  std::vector<std::vector<int>> indices;
  // this is the list of hardjet 4-momenta
  std::vector<math::XYZTLorentzVector> p4_hardJets;
  // this is the hardjet areas
  std::vector<double> area_hardJets;
  std::vector<bool> isHardest;

  edm::Handle<edm::View<T> > pfjets;
  iEvent.getByToken(jetSrc_, pfjets);
  indices.resize(pfjets->size());

  edm::Handle<reco::PFCandidateCollection> pfcands;
  iEvent.getByToken(constitSrc_, pfcands);

  // Grab IP and SV tag info
  edm::Handle<edm::View<reco::BaseTagInfo>> ipTagInfoHandle;
  edm::Handle<std::vector<reco::CandSecondaryVertexTagInfo>> svTagInfoHandle;

  if (typeid(T) == typeid(reco::PFJet)) {
    iEvent.getByToken(tagInfoSrc_, ipTagInfoHandle);
    iEvent.getByToken(svTagInfos_, svTagInfoHandle);
  }

  // Grab gen particles
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleSrc_, genParticles);
  //std::cout << "nb of gen particles in collection: " << (*genParticles).size() << std::endl; 

  // Load trained booster for aggregation
  BoosterHandle booster_;
  XGBoosterCreate(NULL, 0, &booster_);
  if (aggregateHF_) {
    // std::cout << "Loading the trained BDT from " << model_path << std::endl;
    auto retModel = XGBoosterLoadModel(booster_, model_path.fullPath().c_str()); // second argument should be a const char *.
    // if (retModel == 0) std::cout << "Model loaded successfully" << std::endl;
    // else std::cout << "Model was not loaded successfully, return value " << retModel << std::endl;
  }

  // Loop over PFJets
  for (std::size_t jetIndex = 0; jetIndex < (*pfjets).size(); jetIndex++) { 
    const T& pfjet = (*pfjets)[jetIndex];
    // std::cout << "new pfjet with pt: " << pfjet.pt() << std::endl;

    p4_hardJets.push_back(math::XYZTLorentzVector(pfjet.px(), pfjet.py(), pfjet.pz(), pfjet.energy()));
    area_hardJets.push_back(pfjet.jetArea());

    std::vector<fastjet::PseudoJet> jetConstituents;
    fastjet::PseudoJet *subFJ1 = new fastjet::PseudoJet();
    fastjet::PseudoJet *subFJ2 = new fastjet::PseudoJet();
    std::vector<fastjet::PseudoJet> constit1;
    std::vector<fastjet::PseudoJet> constit2;
	
    if (aggregateHF_) {
      if (typeid(T) == typeid(reco::PFJet)) {
		    // std::cout << "Aggregating HF for PFJet" << std::endl;
        // Grab IP and SV tag info for jet
        const std::vector<reco::CandSecondaryVertexTagInfo>& svTagInfos = *svTagInfoHandle;
        const reco::CandSecondaryVertexTagInfo& svTagInfo = svTagInfos[jetIndex];
        const reco::CandIPTagInfo ipTagInfo = *(svTagInfo.trackIPTagInfoRef().get());
        //std::vector<reco::CandidatePtr> ipTracks = ipTagInfo.selectedTracks();
        //std::vector<reco::CandidatePtr> svTracks = svTagInfo.vertexTracks();

        jetConstituents = aggregateHFreco(pfjet, ptCut, ipTagInfo, svTagInfo, *genParticles, booster_);
      } else if (typeid(T) == typeid(reco::GenJet)) {
		    // std::cout << "Aggregating HF for GenJet" << std::endl;
        jetConstituents = aggregateHFgen(pfjet, ptCut, *genParticles);
      }
    } else {
      for (reco::CandidatePtr constituent : pfjet.getJetConstituents()) {
        //std::cout << "Not aggregating" << std::endl;
        //std::cout << "const charge: " << constituent->charge() << std::endl;
        //std::cout << "const pt: " << constituent->pt() << std::endl;
        // charge and pt cut
        if ((chargedOnly_) && (constituent->charge() == 0)) continue;
        if (constituent->pt() < ptCut) continue;
		    //std::cout << "const kept" << std::endl;
        jetConstituents.push_back(fastjet::PseudoJet(constituent->px(), constituent->py(), constituent->pz(), constituent->energy()));
      }
    }
    // std::cout << "nb of jet constituents: " << jetConstituents.size() << std::endl;
    //for (auto constituent : jetConstituents) {
    //std::cout << "eta: " << constituent.eta() << std::endl;
    //std::cout << "phi: " << constituent.phi() << std::endl;
    //}


    // Iterative declustering
    isHardest.push_back(IterativeDeclustering(jetConstituents, subFJ1, subFJ2, constit1, constit2));

    //std::cout << "sjt1 pt: " << subFJ1->pt() << std::endl;
    //std::cout << "sjt2 pt: " << subFJ2->pt() << std::endl;

    // Convert fastjets to basicjets    
    reco::BasicJet subjet1 = ConvertFJ2BasicJet(subFJ1, constit1, pfcands, iSetup);
    reco::BasicJet subjet2 = ConvertFJ2BasicJet(subFJ2, constit2, pfcands, iSetup);

    if (subjet1.pt() > 1.0e-3){
      indices[jetIndex].push_back(subjetCollection->size());
      if (subFJ1->has_area()) subjet1.setJetArea(subFJ1->area());
      subjetCollection->push_back(subjet1);
    }
    if (subjet2.pt() > 1.0e-3){
      indices[jetIndex].push_back(subjetCollection->size());
      if(subFJ2->has_area())subjet2.setJetArea(subFJ2->area());
      subjetCollection->push_back(subjet2);
    }
  }

  // put subjets into event record
  subjetHandleAfterPut = iEvent.put(move(subjetCollection), "SubJets"); 

  // Now create the hard jets with ptr's to the subjets as constituents
  std::vector<math::XYZTLorentzVector>::const_iterator ip4 = p4_hardJets.begin(), 
                                                       ip4Begin = p4_hardJets.begin(), 
                                                       ip4End = p4_hardJets.end();
  
  for (; ip4 != ip4End; ++ip4) {
    int p4_index = ip4 - ip4Begin;
    std::vector<int>& ind = indices[p4_index];
    std::vector<reco::CandidatePtr> i_hardJetConstituents;

    // Add the subjets to the hard jet
    for (std::vector<int>::const_iterator isub = ind.begin(); isub != ind.end(); ++isub) {
      reco::CandidatePtr candPtr(subjetHandleAfterPut, *isub, false);
      i_hardJetConstituents.push_back(candPtr);
    }
    
    reco::Particle::Point point(0, 0, 0);
    //cout<<" size of i_hardJetConstituents "<<i_hardJetConstituents.size()<<endl;
    reco::BasicJet toput(*ip4, point, i_hardJetConstituents);
    if(isHardest[ip4 - ip4Begin]) toput.setJetArea(0.8);
    else toput.setJetArea(0.4);
    //toput.setIsWeighted(isHardest[ip4 - ip4Begin]);  // Not yet in CMSSW_9_X_X
    //cout << "hardjet added to collection" << endl;
    jetCollection->push_back(toput);
  }

  //cout << "jetCollection size: " << jetCollection->size() << endl; 
  iEvent.put(move(jetCollection));
}

template <class T>
reco::BasicJet dynGroomedJets<T>::ConvertFJ2BasicJet(fastjet::PseudoJet *fj, 
                                                     std::vector<fastjet::PseudoJet> constit, 
                                                     edm::Handle<reco::PFCandidateCollection> pfcands, 
                                                     const edm::EventSetup& iSetup) const
{

  math::XYZTLorentzVector p4(fj->px(), fj->py(), fj->pz(), fj->e());  
  reco::Particle::Point point(0, 0, 0);
  std::vector<reco::CandidatePtr> constituents;
  if (writeConstits_){
    for (uint j=0;j<constit.size();j++) {
      double constitPt = constit[j].pt();
      double constitEta = constit[j].eta();
      
      bool foundMatch = false;
      int iCand = -1;
      if (foundMatch) {} // so the compiler does not complain that it is not used
      for (const reco::PFCandidate& pfcand : *pfcands) {
        iCand++;
        
        if(fabs(constitPt - pfcand.pt()) < 0.0001 && fabs(constitEta - pfcand.eta()) <0.0001 ){
          foundMatch = true;
          constituents.push_back(reco::CandidatePtr(pfcands, iCand));
          break;
        }
      }
      //if(!foundMatch) cout<<" couldn't find match "<<endl;
    }
  }

  reco::BasicJet basicjet;
  reco::writeSpecific(basicjet, p4, point, constituents, iSetup);

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
bool dynGroomedJets<T>::IterativeDeclusteringGen(const T& jet, reco::GenParticleCollection genParticles, float ptCut, fastjet::PseudoJet *sub1, fastjet::PseudoJet *sub2, std::vector<fastjet::PseudoJet>& constit1, std::vector<fastjet::PseudoJet>& constit2) const
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
bool dynGroomedJets<T>::IterativeDeclustering(std::vector<fastjet::PseudoJet> jetConstituents, 
                                              fastjet::PseudoJet *sub1, fastjet::PseudoJet *sub2, 
                                              std::vector<fastjet::PseudoJet> &constit1, 
                                              std::vector<fastjet::PseudoJet> &constit2) const
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
  fastjet::JetDefinition jet_def(fastjet::genkt_algorithm, jet_radius_ca, 0, 
                                 static_cast<fastjet::RecombinationScheme>(0), 
                                 fastjet::Best);

  // Return false if no constituents
  if (jetConstituents.size() == 0) return false;

  // Reclustering jet constituents with new algorithm                                                                                          
  try
    {
      fastjet::ClusterSequence csiter(jetConstituents, jet_def);
      std::vector<fastjet::PseudoJet> output_jets = csiter.inclusive_jets(0);
      output_jets = sorted_by_pt(output_jets);
      
      
      fastjet::PseudoJet jj = output_jets[0];
      fastjet::PseudoJet j1;
      fastjet::PseudoJet j2;                                                                                                               
      fastjet::PseudoJet j1first;
      fastjet::PseudoJet j2first;
      
      while (jj.has_parents(j1,j2)) {
	  
        if (j1.perp() < j2.perp()) std::swap(j1, j2);
        
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
        
        if(dyn > kt1){
          nsel = nsplit;
          kt1 = dyn;
        }
        nsplit = nsplit + 1;
        jj = j1;
      }
      
      if (!flagSubjet) *sub1 = output_jets[0];

      if (sub1->has_constituents()) constit1 = sub1->constituents(); 
      if (sub2->has_constituents()) constit2 = sub2->constituents(); 
      if (nsel == nsdin) isHardest = true; 
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
  desc.add<edm::FileInPath>("model_path");
  desc.add<double>("zcut", 0.1);
  desc.add<double>("beta", 0.0);
  desc.add<double>("dynktcut", 1.0);
  desc.add<double>("rParam", 0.4);
  if (typeid(T) == typeid(reco::PFJet)) {
    desc.add<edm::InputTag>("jetSrc", edm::InputTag("ak4PFJets"));
    desc.add<edm::InputTag>("constitSrc", edm::InputTag("particleFlow"));
    desc.add<edm::InputTag>("ak4PFPfImpactParameterTagInfos", edm::InputTag("ak4PFPfImpactParameterTagInfos"));
    desc.add<edm::InputTag>("ak4PFPfInclusiveSecondaryVertexFinderTagInfos", edm::InputTag("ak4PFPfInclusiveSecondaryVertexFinderTagInfos"));
    desc.add<edm::InputTag>("HFdecayProductTagger", edm::InputTag("HFdecayProductTagger"));
    descriptions.add("dynGroomedPFJets", desc);
  }
  else if (typeid(T) == typeid(reco::GenJet)) {
    desc.add<edm::InputTag>("jetSrc", edm::InputTag("ak4GenJets"));
    desc.add<edm::InputTag>("constitSrc", edm::InputTag("genParticles"));
    desc.add<edm::InputTag>("ak4PFPfImpactParameterTagInfos", edm::InputTag("ak4GenJets"));
    desc.add<edm::InputTag>("ak4PFPfInclusiveSecondaryVertexFinderTagInfos", edm::InputTag("ak4GenJets"));
    desc.add<edm::InputTag>("HFdecayProductTagger", edm::InputTag("HFdecayProductTagger"));
    descriptions.add("dynGroomedGenJets", desc);
  }

}

// Function to match reconstructed tracks to gen particles geometrically
template <class T>
int dynGroomedJets<T>::trkGenPartMatch(reco::Jet::Constituent constituent, 
                                       reco::GenParticleCollection genParticles, 
                                       double ptCut) const
{
  // return status of matched particle
  int status = -1;
  
  double dRmin = 100.;
  double dR = dRmin;

  double trkEta = constituent->eta();
  double trkPhi = constituent->phi();

  for (reco::GenParticle genParticle : genParticles) {
    double partPt = genParticle.pt();
    //cout << "genParticle status: " << genParticle.status() << endl;
    if (partPt < ptCut) continue;
	  if (chargedOnly_ && genParticle.charge() == 0) continue;
	
    double partEta = genParticle.eta();
    double partPhi = genParticle.phi();

    double dEta = trkEta - partEta;
    double dPhi = std::acos(std::cos(trkPhi - partPhi));

    dR = std::sqrt((dEta * dEta) + (dPhi * dPhi));

    if (dR < dRmin) {
      dRmin = dR;
      status = genParticle.status();
    }
  }
  //cout << "dRmin: " << dRmin << endl;
  return status;  
}

template <class T>
int dynGroomedJets<T>::trkInVector(reco::CandidatePtr trk, 
                                   std::vector<reco::CandidatePtr> tracks) const {
  int indexInVector = -1;
  auto it = find(tracks.begin(), tracks.end(), trk);
  if (it != tracks.end()) {
	indexInVector = it - tracks.begin();
  }
  return indexInVector;
}


// Function to aggregate reconstructed tracks into pseudo-B's 
template <class T>
std::vector<fastjet::PseudoJet> dynGroomedJets<T>::aggregateHFreco(const T& jet, float ptCut,
                                                                   const reco::CandIPTagInfo ipTagInfo, 
                                                                   const reco::CandSecondaryVertexTagInfo& svTagInfo,
                                                                   const reco::GenParticleCollection genParticles,
                                                                   BoosterHandle booster_) const
{
  // std::cout << "Aggregating Bs in reco jet" << std::endl;
  // Input particle collection 
  std::vector<reco::CandidatePtr> inputJetConstituents = jet.getJetConstituents();
  // Particle collection to use in the declustering
  std::vector<fastjet::PseudoJet> outputJetConstituents;
  // Particle collection to aggregate into a pseudo-B 
  std::map<int, std::vector<reco::CandidatePtr>> hfConstituentsMap;

  // [DEBUG]
  // std::string path = "../../../RecoHI/HiJetAlgos/data/";
  // DIR *dir;
  // struct dirent *ent;
  // if ((dir = opendir (path.c_str())) != NULL) {
  //   /* print all the files and directories within directory */
  //   while ((ent = readdir (dir)) != NULL) {
  //     std::cout << "**DEBUG** " << ent->d_name << "\n" << std::endl;
  //   }
  //   closedir (dir);
  // } else {
  //   /* could not open directory */
  //   std::cout << "**DEBUG** Could not open directory " << path << std::endl;
  // }


  // [DEBUG]: Prepare data to make predictions
  // const int nFeatures_ = 4;
  // const int nEntries_ = 3;
  // const float missing_value_ = -1000000.;
  // float TestData_[nEntries_][nFeatures_] = {{1., 2., 3., 5.},
  //                                        {0., missing_value_, missing_value_, 1.},
  //                                        {0., missing_value_, missing_value_, 10.}}; 
  // Assign data to the "TestData" 2d array ... 
  // Declare object & allocate memory & initialize
  // XGDMatrixCreateFromMat((float *)TestData_, nEntries_, nFeatures_, missing_value_, &data_); 
  
  // [DEBUG]: Make predictions on TestData_
  
  // std::cout << "**DEBUG** out_len before: " << out_len << endl;
  // auto ret = XGBoosterPredict(booster_, data_, 0, 0, &out_len, &f);// lower version API
  // std::cout << "**DEBUG** ret: " << ret << std::endl;
  // std::cout << "**DEBUG** ret=0 success; ret=-1 failure " << std::endl;
  // std::cout << "**DEBUG** out_len after: " << out_len << endl;
  // std::cout << "**DEBUG** Model prediction for:";
  // for (int i = 0; i < nEntries_; i++) {
  //   std::cout << "\n{ ";
  //   for (int j = 0; j < nFeatures_; j++) {
  //     std::cout << TestData_[i][j] << " ";
  //   }
  //   std::cout << "} is: " << f[i] << std::endl;
  // }

  // Prepare the data array to use for the tracks
  bst_ulong out_len = 0; // bst_ulong is a typedef of unsigned long
  const float *f; // array to store predictions
  DMatrixHandle data_;

  const int nFeatures = 4;
  const int nEntries = 1;
  const double missing_value = -1000000.;
  float TestData[nEntries][nFeatures];
  float threshold = 0.71;

  // Get ipTracks and svTracks for jet
  std::vector<reco::CandidatePtr> ipTracks = ipTagInfo.selectedTracks();
  std::vector<reco::CandidatePtr> svTracks = svTagInfo.vertexTracks();

  for(reco::CandidatePtr constituent : inputJetConstituents) {
    // Charge and pt cut
    if(chargedOnly_ && constituent->charge() == 0) continue;
    if(constituent->pt() < ptCut) continue;

    // Look for track in ipTracks, otherwise discard
    int whichTrack = trkInVector(constituent, ipTracks);
    if (whichTrack < 0) continue;
    // std::cout << "New track!" << std::endl;

    // Look for track in SVs 
    float inSV = 0.;
    float svtxdls = missing_value;
    float svtxm = missing_value;
    float ip3dSig = missing_value;

    for (uint ivtx = 0; ivtx < svTagInfo.nVertices(); ivtx++) {
      std::vector<reco::CandidatePtr> isvTracks = svTagInfo.vertexTracks(ivtx);
      int whichTrackInSV = trkInVector(constituent, isvTracks);
      if (whichTrackInSV < 0) continue;

      // std::cout << "Track found in vtx " << ivtx << std::endl;
      // if (inSV == 1) std::cout << "********Track found in 2 vertices***********" << std::endl;
      inSV = 1;
      
      Measurement1D m1D = svTagInfo.flightDistance(ivtx);
      svtxdls = m1D.significance();
      Measurement1D m2D = svTagInfo.flightDistance(ivtx, true);
      const reco::VertexCompositePtrCandidate svtx = svTagInfo.secondaryVertex(ivtx);
      double svtxM = svtx.p4().mass();
      svtxm = svtxM; 
      break;
    } // end vtx loop

    // get ip3dSig
    reco::btag::TrackIPData IPdata = ipTagInfo.impactParameterData()[whichTrack];
    ip3dSig = IPdata.ip3d.significance();

    // std::cout << "track pt = " << constituent->pt();
    // std::cout << ", inSV = " << inSV;
    // std::cout << ", svtxdls = " << svtxdls;
    // std::cout << ", svtxm = " << svtxm;
    // std::cout << ", ip3dSig = " << ip3dSig << std::endl;

    int hfIndex = 100;
    bool isHFproduct = false;

    // ***********************************
    // **** USING CUT BASED SELECTION ****
    // ***********************************

    // isHFproduct = (inSV && (std::abs(ip3dSig) > 3.)) || (std::abs(ip3dSig) > 9.);

    // ***********************************
    // **** USING THE TRAINED BOOSTER ****
    // ***********************************

    TestData[0][0] = inSV;
    TestData[0][1] = svtxdls;
    TestData[0][2] = svtxm;
    TestData[0][3] = ip3dSig;
    XGDMatrixCreateFromMat((float *)TestData, nEntries, nFeatures, missing_value, &data_);
    auto retPred = XGBoosterPredict(booster_, data_, 0, 0, &out_len, &f);
    // if (retPred == 0) std::cout << "Booster made a prediction successfully" << std::endl;
    // else std::cout << "Booster could not make a prediction, return value " << retPred << std::endl;

    // std::cout << "Booster predicted probability: " << f[0] << std::endl;
    isHFproduct = (f[0] > threshold);
    // std::cout << isHFproduct << std::endl;

    // ***********************************
    // ****      USING TRUTH INFO     ****
    // ***********************************

    // int genStatus = trkGenPartMatch(constituent, genParticles, ptCut);
    // std::cout << "Truth status = " << genStatus << std::endl;
    // if (genStatus >= 100) {
    //   hfIndex = genStatus;
    //   isHFproduct = true;
    // }

    if (isHFproduct) {
      // std::cout << "Found HF constituent" << std::endl;
      hfConstituentsMap[hfIndex].push_back(constituent);
    } else {
      // std::cout << "Not a HF constituent, added to collection" << std::endl;
      outputJetConstituents.push_back(fastjet::PseudoJet(constituent->px(), constituent->py(), 
                                                         constituent->pz(), constituent->energy()));
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
    // std::cout << "Adding HF to the collection" << std::endl;
    outputJetConstituents.push_back(fastjet::PseudoJet(pseudoHF.Px(), pseudoHF.Py(), pseudoHF.Pz(), pseudoHF.E()));
  } // end map loop
 
  return outputJetConstituents;
} // end aggregateHF()

// Function to aggregate gen particles coming from B decays into pseudo-B's 
template <class T>
std::vector<fastjet::PseudoJet> dynGroomedJets<T>::aggregateHFgen(const T& jet, float ptCut,
                                                                  reco::GenParticleCollection genParticles) const
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
		    //std::cout << "Found HF constituent" << std::endl;
      } else if (genStatus == 1) {
        outputConstituents.push_back(fastjet::PseudoJet(constituent->px(), constituent->py(), constituent->pz(), constituent->energy()));
		    //std::cout << "Not a HF constituent, added to collection" << std::endl;
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
      //std::cout << "added HF to collection" << std::endl;
      outputConstituents.push_back(fastjet::PseudoJet(pseudoHF.Px(), pseudoHF.Py(), pseudoHF.Pz(), pseudoHF.E()));
    }
  return outputConstituents;
}
  


using dynGroomedGenJets = dynGroomedJets<reco::GenJet>;
using dynGroomedPFJets = dynGroomedJets<reco::PFJet>;

// define this as a plug-in
DEFINE_FWK_MODULE(dynGroomedPFJets);
DEFINE_FWK_MODULE(dynGroomedGenJets);
