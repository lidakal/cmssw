
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
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/ShallowTagInfo.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "RecoJets/JetProducers/interface/JetSpecific.h"
#include "CommonTools/UtilAlgos/interface/DeltaR.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include <fastjet/AreaDefinition.hh>
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/SoftDrop.hh"


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

  bool IterativeDeclustering(vector<PseudoJet>, PseudoJet *, PseudoJet *, vector<PseudoJet> &, vector<PseudoJet> &) const;
  BasicJet ConvertFJ2BasicJet(PseudoJet *, vector<PseudoJet>, Handle<PFCandidateCollection>, const EventSetup&) const;
  BasicJet ConvertFJ2BasicJet(PseudoJet *, vector<PseudoJet>, Handle<edm::View<pat::PackedCandidate> >, const EventSetup&) const;

  int trkInVector(reco::CandidatePtr, std::vector<reco::CandidatePtr>) const;
  int trkGenPartMatch(reco::Jet::Constituent, reco::GenParticleCollection, double) const;
  vector<PseudoJet> aggregateHF(const T&, float , reco::GenParticleCollection) const;
  vector<PseudoJet> aggregateHF(const T&, float, const reco::CandIPTagInfo, const reco::CandSecondaryVertexTagInfo&, const reco::GenParticleCollection) const;

  EDGetTokenT<View<T> > jetSrc_;
  EDGetTokenT<PFCandidateCollection> constitSrc_;
  EDGetTokenT<edm::View<pat::PackedCandidate>>  packedConstitSrc_;
  EDGetTokenT<View<BaseTagInfo>> tagInfoSrc_;
  const edm::EDGetTokenT<std::vector<reco::CandSecondaryVertexTagInfo>> svTagInfos_;
  const edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticleSrc_;

  bool writeConstits_;
  bool doLateSD_;
  bool chargedOnly_;
  bool aggregateHF_;
  double constPtCut_;
  double zcut_;
  double beta_;
  double dynktcut_;
  double rParam_;
};

template <class T>
dynGroomedJets<T>::dynGroomedJets(const ParameterSet& iConfig)
  :  jetSrc_(consumes<View<T> >(iConfig.getParameter<InputTag>("jetSrc"))),
     constitSrc_(mayConsume<PFCandidateCollection >(iConfig.getParameter<InputTag>("constitSrc"))),    
     packedConstitSrc_(mayConsume<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("constitSrc"))),
     svTagInfos_(mayConsume<std::vector<reco::CandSecondaryVertexTagInfo>>(iConfig.getParameter<edm::InputTag>("pfInclusiveSecondaryVertexFinderTagInfos"))),
     genParticleSrc_(mayConsume<vector<GenParticle>>(iConfig.getParameter<InputTag>("HFdecayProductTagger"))),
  writeConstits_(iConfig.getParameter<bool>("writeConstits")),
  doLateSD_(iConfig.getParameter<bool>("doLateSD")),
  chargedOnly_(iConfig.getParameter<bool>("chargedOnly")),
  aggregateHF_(iConfig.getParameter<bool>("aggregateHF")), 
  constPtCut_(iConfig.getParameter<double>("constPtCut")),
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

  Handle<View<T> > jets;
  iEvent.getByToken(jetSrc_, jets);

  Handle<PFCandidateCollection> pfcands;
  bool isPF = iEvent.getByToken(constitSrc_, pfcands);
 
  Handle<edm::View<pat::PackedCandidate> > pfcandsPacked;
  bool isPackedPF = iEvent.getByToken(packedConstitSrc_, pfcandsPacked);

  // stuff for aggregation
  // Grab SV tag info
  edm::Handle<std::vector<reco::CandSecondaryVertexTagInfo>> svTagInfoHandle;
  if (typeid(T) == typeid(reco::PFJet)) iEvent.getByToken(svTagInfos_, svTagInfoHandle);
  // Grab gen particles
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleSrc_, genParticles);


  indices.resize(jets->size());


  int jetIndex = 0;
  for (const T& jet : *jets) {
    p4_hardJets.push_back(math::XYZTLorentzVector(jet.px(), jet.py(), jet.pz(), jet.energy()));
    area_hardJets.push_back(jet.jetArea());

    vector<PseudoJet> jetConstituents;
    PseudoJet *subFJ1 = new PseudoJet();
    PseudoJet *subFJ2 = new PseudoJet();
    vector<PseudoJet> constit1;
    vector<PseudoJet> constit2;

    // aggregate B based on jet type
    if (aggregateHF_) {
      if (typeid(T) == typeid(reco::PFJet)) {
        // Grab IP and SV tag info for jet
        const std::vector<reco::CandSecondaryVertexTagInfo>& svTagInfos = *svTagInfoHandle;
        const reco::CandSecondaryVertexTagInfo& svTagInfo = svTagInfos[jetIndex];
        const reco::CandIPTagInfo ipTagInfo = *(svTagInfo.trackIPTagInfoRef().get());
        jetConstituents = aggregateHF(jet, constPtCut_, ipTagInfo, svTagInfo, *genParticles);
      } else if (typeid(T) == typeid(reco::GenJet)) {
        jetConstituents = aggregateHF(jet, constPtCut_, *genParticles);
      }
    } else {
      for (reco::CandidatePtr constituent : jet.getJetConstituents()) {
        // charge and pt cut
        if ((chargedOnly_) && (constituent->charge() == 0)) continue;
        if (constituent->pt() < constPtCut_) continue;
        jetConstituents.push_back(PseudoJet(constituent->px(), constituent->py(), constituent->pz(), constituent->energy()));
      }
    }

    isHardest.push_back(IterativeDeclustering(jetConstituents,subFJ1,subFJ2,constit1,constit2));
    
    BasicJet subjet1, subjet2;

    if(isPF){
      subjet1 = ConvertFJ2BasicJet(subFJ1, constit1, pfcands, iSetup);
      subjet2 = ConvertFJ2BasicJet(subFJ2, constit2, pfcands, iSetup);
    }
    else if(isPackedPF){
      subjet1 = ConvertFJ2BasicJet(subFJ1, constit1, pfcandsPacked, iSetup);
      subjet2 = ConvertFJ2BasicJet(subFJ2, constit2, pfcandsPacked, iSetup);
    }
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
    BasicJet toput(*ip4, point, i_hardJetConstituents);
    if(isHardest[ip4-ip4Begin]) toput.setJetArea(0.8);  // super hacky
    else toput.setJetArea(0.4);
    //toput.setIsWeighted(isHardest[ip4 - ip4Begin]);  // Not yet in 10_6_X
    jetCollection->push_back(toput);
  }

  
  iEvent.put(move(jetCollection));
}

template <class T>
BasicJet dynGroomedJets<T>::ConvertFJ2BasicJet(PseudoJet *fj, vector<PseudoJet> constit, Handle<edm::View<pat::PackedCandidate> > pfcandsPacked, const EventSetup& iSetup) const
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
      
      for (const pat::PackedCandidate& pfcandPacked : *pfcandsPacked) {
	iCand++;
	
	if(fabs(constitPt - pfcandPacked.pt()) < 0.0001 && fabs(constitEta - pfcandPacked.eta()) <0.0001 ){
	  foundMatch = true;
	  constituents.push_back(reco::CandidatePtr(pfcandsPacked, iCand));
	  break;
	}
      }
      if(!foundMatch) cout<<" couldn't find match "<<endl;
    }
  }

  BasicJet basicjet;
  writeSpecific(basicjet, p4, point, constituents, iSetup);

  return basicjet;

}

template <class T>
BasicJet dynGroomedJets<T>::ConvertFJ2BasicJet(PseudoJet *fj, vector<PseudoJet> constit, Handle<PFCandidateCollection> pfcandsPacked, const EventSetup& iSetup) const
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
      
      for (const PFCandidate& pfcand : *pfcandsPacked) {
	iCand++;
	
	if(fabs(constitPt - pfcand.pt()) < 0.0001 && fabs(constitEta - pfcand.eta()) <0.0001 ){
	  foundMatch = true;
	  constituents.push_back(reco::CandidatePtr(pfcandsPacked, iCand));
	  break;
	}
      }
      if(!foundMatch) cout<<" couldn't find match "<<endl;
    }
  }

  BasicJet basicjet;
  writeSpecific(basicjet, p4, point, constituents, iSetup);

  return basicjet;

}



template <class T>
bool dynGroomedJets<T>::IterativeDeclustering(vector<PseudoJet> jetConstituents, 
                                              PseudoJet *sub1, PseudoJet *sub2, 
                                              vector<PseudoJet> &constit1, vector<PseudoJet> &constit2) const
{


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

// Function to match reconstructed tracks to gen particles geometrically
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
    if (partPt < ptCut) continue;
    if(chargedOnly_ && genParticle.charge() == 0) continue;
	
    double partEta = genParticle.eta();
    double partPhi = genParticle.phi();

    double dEta = trkEta - partEta;
    double dPhi = std::acos(std::cos(trkPhi - partPhi));

    //dR = std::sqrt((dEta * dEta) + (dPhi * dPhi));
    dR = (dEta * dEta) + (dPhi * dPhi);  //faster if you don't take the sqrt

    if (dR < dRmin) {
      dRmin = dR;
      status = genParticle.status();
    }
  }
  return status;  
}

template <class T>
int dynGroomedJets<T>::trkInVector(reco::CandidatePtr trk, std::vector<reco::CandidatePtr> tracks) const {
  int indexInVector = -1;
  auto it = find(tracks.begin(), tracks.end(), trk);
  if (it != tracks.end()) indexInVector = it - tracks.begin();

  return indexInVector;
}


// Function to aggregate reconstructed tracks into pseudo-B's 
template <class T>
std::vector<fastjet::PseudoJet> dynGroomedJets<T>::aggregateHF(const T& jet, float ptCut,
                                                              const reco::CandIPTagInfo ipTagInfo, 
                                                              const reco::CandSecondaryVertexTagInfo& svTagInfo,
                                                              const reco::GenParticleCollection genParticles) const
{
  // Input particle collection 
  std::vector<reco::CandidatePtr> inputJetConstituents = jet.getJetConstituents();

  // Particle collection to use in the declustering
  std::vector<fastjet::PseudoJet> outputJetConstituents;

  // Particle collection to aggregate into a pseudo-B 
  std::map<int, std::vector<reco::CandidatePtr>> hfConstituentsMap;

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

    /*
    // Look for track in svTracks
    int whichTrackInSV = trkInVector(constituent, svTracks);
    bool inSV = (whichTrackInSV >= 0);

    if (inSV) Measurement1D m1D = svTagInfo.flightDistance(0);
    
    // get ip3dSig
    reco::btag::TrackIPData IPdata = ipTagInfo.impactParameterData()[whichTrack];
    Float_t ip3dSig = IPdata.ip3d.significance();
    */
    // use truth info 
    int genStatus = trkGenPartMatch(constituent, genParticles, ptCut);
    // Collect particles coming from HF decays, add the rest directly to the collection
    if (genStatus >= 100) {
      hfConstituentsMap[genStatus].push_back(constituent);
      //   	//std::cout << "Found HF constituent" << std::endl;
    } else if (genStatus == 1) {
      outputJetConstituents.push_back(PseudoJet(constituent->px(), constituent->py(), constituent->pz(), constituent->energy()));
      // 	  //std::cout << "Not a B constituent, added to collection" << std::endl;
    }
  } // end loop over input jet constituents

  // Aggregate particles coming from HF decays into pseudo-B/C's 
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
                                                              reco::GenParticleCollection genParticles) const
{
  // Input and output particle collections
  std::vector<reco::CandidatePtr> inputConstituents = jet.getJetConstituents();
  std::vector<fastjet::PseudoJet> outputConstituents;
  std::map<int, std::vector<reco::CandidatePtr>> hfConstituentsMap;
  
  for (reco::CandidatePtr constituent : inputConstituents) {
    // Charge and pt cut
    if(!constituent.isNonnull()){
      cout<<" NULL constituent !!! "<<endl;
      continue;
    }
    if(chargedOnly_ && constituent->charge() == 0) continue;
    if(constituent->pt() < ptCut) continue;
    
    // Match with gen particle -- should be exact match
    int genStatus = trkGenPartMatch(constituent, genParticles, ptCut);
    // Collect particles coming from HF decays, add the rest directly to the collection
    if (genStatus >= 100) hfConstituentsMap[genStatus].push_back(constituent);
    else if (genStatus == 1) 
      outputConstituents.push_back(PseudoJet(constituent->px(), constituent->py(), constituent->pz(), constituent->energy()));
  }

  // Aggregate particles coming from HF decays into pseudo-B/D's 
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
    outputConstituents.push_back(PseudoJet(pseudoHF.Px(), pseudoHF.Py(), pseudoHF.Pz(), pseudoHF.E()));
  }
  return outputConstituents;
}
  



template <class T>
void dynGroomedJets<T>::fillDescriptions(ConfigurationDescriptions& descriptions) {
  ParameterSetDescription desc;
  desc.setComment("Dynamically groomed jets");
  desc.add<bool>("writeConstits", false);
  desc.add<bool>("doLateSD", false);
  desc.add<bool>("chargedOnly", false);
  desc.add<bool>("aggregateHF", true);
  desc.add<double>("constPtCut",1.0);
  desc.add<double>("zcut", 0.1);
  desc.add<double>("beta", 0.0);
  desc.add<double>("dynktcut", 1.0);
  desc.add<double>("rParam", 0.4);
  if (typeid(T) == typeid(pat::Jet)) {
    desc.add<edm::InputTag>("jetSrc", edm::InputTag("slimmedJets"));
    desc.add<edm::InputTag>("constitSrc", edm::InputTag("packedPFCandidates"));
    desc.add<edm::InputTag>("pfInclusiveSecondaryVertexFinderTagInfos", edm::InputTag("pfInclusiveSecondaryVertexFinderTagInfos"));
    desc.add<edm::InputTag>("HFdecayProductTagger", edm::InputTag("HFdecayProductTagger"));
    descriptions.add("dynGroomedPFJets", desc);
  }
  else if (typeid(T) == typeid(reco::GenJet)) {
    desc.add<edm::InputTag>("jetSrc", edm::InputTag("slimmedGenJets"));
    desc.add<edm::InputTag>("constitSrc", edm::InputTag("packedGenParticles"));
    desc.add<edm::InputTag>("pfInclusiveSecondaryVertexFinderTagInfos", edm::InputTag("slimmedGenJets"));  //dummy
    desc.add<edm::InputTag>("HFdecayProductTagger", edm::InputTag("HFdecayProductTagger"));
    descriptions.add("dynGroomedGenJets", desc);
  }

}


using dynGroomedGenJets = dynGroomedJets<pat::Jet>;
using dynGroomedPFJets = dynGroomedJets<reco::GenJet>;

// define this as a plug-in
DEFINE_FWK_MODULE(dynGroomedPFJets);
DEFINE_FWK_MODULE(dynGroomedGenJets);
