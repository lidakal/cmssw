
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
  // Still used for gen jets
  // PseudoJet is fastjet::PseudoJet
  bool IterativeDeclustering(const T&, PseudoJet *, PseudoJet *, vector<PseudoJet> &, vector<PseudoJet> &) const;
  // Used for reco jets (vector<CandidatePtr> are the tracks)
  bool IterativeDeclustering(const T&, vector<CandidatePtr>, reco::GenParticleCollection, PseudoJet *, PseudoJet *, vector<PseudoJet> &, vector<PseudoJet> &) const;

  BasicJet ConvertFJ2BasicJet(PseudoJet *fj, vector<PseudoJet>, Handle<PFCandidateCollection>, const EventSetup& iSetup) const;

  int trkGenPartMatch(reco::Jet::Constituent constituent, reco::GenParticleCollection genParticles, double ptCut) const;

  EDGetTokenT<View<T> > jetSrc_;
  EDGetTokenT<PFCandidateCollection> constitSrc_;
  EDGetTokenT<View<BaseTagInfo>> tagInfoSrc_;
  EDGetTokenT<vector<GenParticle>> genParticleSrc_;
  bool writeConstits_;
  bool doLateSD_;
  bool chargedOnly_;
  double zcut_;
  double beta_;
  double dynktcut_;
  double rParam_;
};

template <class T>
dynGroomedJets<T>::dynGroomedJets(const ParameterSet& iConfig)
  : jetSrc_(consumes<View<T> >(iConfig.getParameter<InputTag>("jetSrc"))),
    constitSrc_(mayConsume<PFCandidateCollection >(iConfig.getParameter<InputTag>("constitSrc"))),
    tagInfoSrc_(mayConsume<View<BaseTagInfo>>(iConfig.getParameter<InputTag>("ak4PFPfImpactParameterTagInfos"))),
    genParticleSrc_(mayConsume<vector<GenParticle>>(iConfig.getParameter<InputTag>("CheatHFHadronReplacer"))),
    writeConstits_(iConfig.getParameter<bool>("writeConstits")),
    doLateSD_(iConfig.getParameter<bool>("doLateSD")),
    chargedOnly_(iConfig.getParameter<bool>("chargedOnly")),
    zcut_(iConfig.getParameter<double>("zcut")),
    beta_(iConfig.getParameter<double>("beta")),
    dynktcut_(iConfig.getParameter<double>("dynktcut")),
    rParam_(iConfig.getParameter<double>("rParam")){
  string alias = (iConfig.getParameter<InputTag>("jetSrc")).label();
  produces<BasicJetCollection>().setBranchAlias(alias);
  produces<vector<BasicJet>>("SubJets").setBranchAlias(alias);
}

template <class T>
void dynGroomedJets<T>::produce(StreamID, Event& iEvent, const EventSetup& iSetup) const {

  auto jetCollection = make_unique<BasicJetCollection>();
  auto subjetCollection = make_unique<vector<BasicJet>>();

  // This will store the handle for the subjets after we write them
  OrphanHandle< vector<BasicJet> > subjetHandleAfterPut;
  // this is the mapping of subjet to hard jet
  vector< vector<int> > indices;
  // this is the list of hardjet 4-momenta
  vector<math::XYZTLorentzVector> p4_hardJets;
  // this is the hardjet areas
  vector<double> area_hardJets;
  vector<bool> isHardest;

  Handle<View<T> > pfjets;
  iEvent.getByToken(jetSrc_, pfjets);

  Handle<PFCandidateCollection> pfcands;
  iEvent.getByToken(constitSrc_, pfcands);
 
  Handle<View<BaseTagInfo>> ipTagInfoHandle;
  if (typeid(T) == typeid(reco::PFJet)) iEvent.getByToken(tagInfoSrc_, ipTagInfoHandle);

  // std::vector<reco::GenParticle> = reco::GenParticleCollection typedef in GenParticleFwd.h
  Handle<vector<GenParticle>> genParticles;
  if (typeid(T) == typeid(reco::PFJet)) iEvent.getByToken(genParticleSrc_, genParticles);


  indices.resize(pfjets->size());

  int jetIndex = 0;
  for (const T& pfjet : *pfjets) {
    p4_hardJets.push_back(math::XYZTLorentzVector(pfjet.px(), pfjet.py(), pfjet.pz(), pfjet.energy()));
    area_hardJets.push_back(pfjet.jetArea());

    PseudoJet *subFJ1 = new PseudoJet();
    PseudoJet *subFJ2 = new PseudoJet();
    vector<PseudoJet> constit1;
    vector<PseudoJet> constit2;

    if (typeid(T) == typeid(reco::PFJet)) { 

      const edm::View<reco::BaseTagInfo> &ipTagInfos = *ipTagInfoHandle;
      const edm::RefToBase<Jet>& ipJetRef = ipTagInfos[jetIndex].jet();

      const reco::IPTagInfo<std::vector<CandidatePtr>, JetTagInfo>* ipTagInfo =
	dynamic_cast<const reco::IPTagInfo<vector<CandidatePtr>, JetTagInfo>*>(&ipTagInfos[jetIndex]);

      vector<reco::CandidatePtr> ipTracks = ipTagInfo->selectedTracks();
      //cout<<" nsel tracks "<<ipTracks.size()<<endl;
      //cout<<" nsel tracks = "<<ipTagInfo->selectedTracks().size()<<endl;
      isHardest.push_back(IterativeDeclustering(pfjet,ipTracks,*genParticles,subFJ1,subFJ2,constit1,constit2));
      /*
	// to loop over gen particles

      for(auto genParticle : *genParticles) {
	cout<<" I'm a gen particle "<<genParticle.phi()<<endl;
      }
      */

    }
    else{
      isHardest.push_back(IterativeDeclustering(pfjet,subFJ1,subFJ2,constit1,constit2));
    }
    //cout<<" jet pT = "<<pfjet.pt()<<" sj1 pT = "<<subFJ1->pt()<<" sj2 pT = "<<subFJ2->pt()<<endl;
    
    BasicJet subjet1 = ConvertFJ2BasicJet(subFJ1, constit1, pfcands, iSetup);
    BasicJet subjet2 = ConvertFJ2BasicJet(subFJ2, constit2, pfcands, iSetup);

    if(subjet1.pt()>1.0e-3){
      indices[jetIndex].push_back(subjetCollection->size());
      if(subFJ1->has_area())subjet1.setJetArea(subFJ1->area());
      subjetCollection->push_back(subjet1);
    }
    if(subjet2.pt()>1.0e-3){
      indices[jetIndex].push_back(subjetCollection->size());
      if(subFJ2->has_area())subjet2.setJetArea(subFJ2->area());
      subjetCollection->push_back(subjet2);
    }
    // put subjets into event record
    jetIndex++;
  }
  
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
    jetCollection->push_back(toput);
  }

  
  iEvent.put(move(jetCollection));
}

template <class T>
BasicJet dynGroomedJets<T>::ConvertFJ2BasicJet(PseudoJet *fj, vector<PseudoJet> constit, Handle<PFCandidateCollection> pfcands, const EventSetup& iSetup) const
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


template <class T>
bool dynGroomedJets<T>::IterativeDeclustering(const T& jet, vector<CandidatePtr> ipTracks, reco::GenParticleCollection genParticles, PseudoJet *sub1, PseudoJet *sub2, vector<PseudoJet> &constit1, vector<PseudoJet> &constit2) const
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
      std::map<int, reco::Jet::Constituents> pseudoBs;

      int nCharged = 0;
      // getJetConstituents() returns std::vector<reco::Jet::Constituent> where reco::Jet::Constituent = edm::Ptr<Candidate>
      auto daughters = jet.getJetConstituents();
      //std::cout<< "Constituents in jet: " << daughters.size() << std::endl;

      //cout << "new jet" << endl;
      // daughters is a vector of pointers
      for(reco::Jet::Constituent daughter : daughters) {
	if(chargedOnly_ && daughter->charge() == 0) continue;
	if(daughter->pt()<1.0) continue;
	float eps = 0.001;
	bool foundTrack = false;
	//cout<<" daughter pT "<<daughter->pt()<<endl;
	for(auto trk : ipTracks){
	  if(fabs(trk->eta()-daughter->eta())>eps) continue;
	  if(acos(cos(trk->phi()-daughter->phi()))>eps) continue;
	  foundTrack = true;
	  break;
	}
	//cout<< "foundTrack = " << foundTrack << endl;
	if(!foundTrack) continue;
	// add particle matching, skip particles matched as B decay prods
	// collect all B constituents
	int trkStatus = trkGenPartMatch(daughter, genParticles, 1.);
	//cout << "trkStatus : " << trkStatus << endl;
	if (trkStatus < 0) continue; // skip tracks with no match - impossible at the moment
	if (trkStatus >= 100) {
	  pseudoBs[trkStatus].push_back(daughter);
	  //std::cout << "Found B product with status " << trkStatus << std::endl;
	  nCharged++;
	}
	if (trkStatus > 1) continue;
	particles.push_back(PseudoJet(daughter->px(), daughter->py(), daughter->pz(), daughter->energy()));
	nCharged++;
      }

      for (auto it = pseudoBs.begin(); it != pseudoBs.end(); ++it) {
	double px = 0;
	double py = 0;
	double pz = 0;
	double energy = 0;
	for (reco::Jet::Constituent bDaughter : it->second) {
	  px += bDaughter->px();
	  py += bDaughter->py();
	  pz += bDaughter->pz();
	  energy += bDaughter->energy();
	}
	//cout << "Added B to particles" << endl;
	particles.push_back(PseudoJet(px, py, pz, energy));
      }

      //std::cout<<" nCharged =  "<<nCharged<<std::endl;
      if(nCharged == 0) return false;
      //cout << "jet kept" << endl;
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
    } catch (fastjet::Error) { /*return -1;*/ }

  return isHardest;
}




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
    } catch (fastjet::Error) { /*return -1;*/ }

  return isHardest;
}




template <class T>
void dynGroomedJets<T>::fillDescriptions(ConfigurationDescriptions& descriptions) {
  ParameterSetDescription desc;
  desc.setComment("Dynamically groomed jets");
  desc.add<bool>("writeConstits", false);
  desc.add<bool>("doLateSD", false);
  desc.add<bool>("chargedOnly", false);
  desc.add<double>("zcut", 0.1);
  desc.add<double>("beta", 0.0);
  desc.add<double>("dynktcut", 1.0);
  desc.add<double>("rParam", 0.4);
  if (typeid(T) == typeid(reco::PFJet)) {
    desc.add<edm::InputTag>("jetSrc", edm::InputTag("ak4PFJets"));
    desc.add<edm::InputTag>("constitSrc", edm::InputTag("particleFlow"));
    desc.add<edm::InputTag>("ak4PFPfImpactParameterTagInfos", edm::InputTag("ak4PFPfImpactParameterTagInfos"));
    desc.add<edm::InputTag>("CheatHFHadronReplacer", edm::InputTag("CheatHFHadronReplacer"));
    descriptions.add("dynGroomedPFJets", desc);
  }
  else if (typeid(T) == typeid(reco::GenJet)) {
    desc.add<edm::InputTag>("jetSrc", edm::InputTag("ak4GenJets"));
    desc.add<edm::InputTag>("constitSrc", edm::InputTag("genParticles"));
    desc.add<edm::InputTag>("ak4PFPfImpactParameterTagInfos", edm::InputTag("ak4GenJets"));
    desc.add<edm::InputTag>("CheatHFHadronReplacer", edm::InputTag("ak4GenJets"));
    descriptions.add("dynGroomedGenJets", desc);
  }

}

//reco::Jet::Constituent aka edm::Ptr<reco::Candidate> is a pointer
template <class T>
int dynGroomedJets<T>::trkGenPartMatch(reco::Jet::Constituent constituent, reco::GenParticleCollection genParticles, double ptCut) const
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
  return status;  
}


using dynGroomedGenJets = dynGroomedJets<reco::PFJet>;
using dynGroomedPFJets = dynGroomedJets<reco::GenJet>;

// define this as a plug-in
DEFINE_FWK_MODULE(dynGroomedPFJets);
DEFINE_FWK_MODULE(dynGroomedGenJets);
