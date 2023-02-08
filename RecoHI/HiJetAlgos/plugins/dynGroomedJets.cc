
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
#include "AnalysisDataFormats/TrackInfo/interface/TrackToGenParticleMap.h"

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
                             std::vector<fastjet::PseudoJet>&, 
                             std::vector<fastjet::PseudoJet>&) const;
  reco::BasicJet ConvertFJ2BasicJet(fastjet::PseudoJet *, 
                                    std::vector<fastjet::PseudoJet>, 
                                    edm::Handle<std::vector<reco::PFCandidate>>, 
                                    const edm::EventSetup&) const;
  reco::BasicJet ConvertFJ2BasicJet(fastjet::PseudoJet *, 
                                    std::vector<fastjet::PseudoJet>, 
                                    edm::Handle<edm::View<pat::PackedCandidate>>, 
                                    const edm::EventSetup&) const;

  
  typedef std::pair<std::vector<fastjet::PseudoJet>, pat::PackedCandidate> jetConstituentsPseudoHFPair;
  jetConstituentsPseudoHFPair aggregateHFGen(const T&, 
                                             reco::TrackToGenParticleMap) const;
  jetConstituentsPseudoHFPair aggregateHFReco(const T&, 
                                              reco::TrackToGenParticleMap, 
                                              BoosterHandle) const;

  // ------------- member data ----------------------------
  edm::EDGetTokenT<std::vector<T>> jetSrc_;
  edm::EDGetTokenT<std::vector<reco::PFCandidate>> constitSrc_;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>>  packedConstitSrc_;
  edm::EDGetTokenT<reco::TrackToGenParticleMap> candToGenParticleMapToken_;

  bool writeConstits_;
  bool doLateSD_;
  bool chargedOnly_;
  bool doGenJets_;
  
  double zcut_;
  double beta_;
  double dynktcut_;
  double rParam_;
  double ptCut_;

  bool aggregateHF_;
  bool withTruthInfo_;
  bool withCuts_;
  bool withBDT_;
  edm::FileInPath model_path_;

  std::string ipTagInfoLabel_;
  std::string svTagInfoLabel_;
};

template <class T>
dynGroomedJets<T>::dynGroomedJets(const edm::ParameterSet& iConfig) {
  // Get configuration parameters
  writeConstits_ = iConfig.getParameter<bool>("writeConstits");
  doLateSD_ = iConfig.getParameter<bool>("doLateSD");
  chargedOnly_ = iConfig.getParameter<bool>("chargedOnly");
  aggregateHF_ = iConfig.getParameter<bool>("aggregateHF"); 
  doGenJets_ = iConfig.getParameter<bool>("doGenJets");
  
  ptCut_ = iConfig.getParameter<double>("ptCut");
  zcut_ = iConfig.getParameter<double>("zcut");
  beta_ = iConfig.getParameter<double>("beta");
  dynktcut_ = iConfig.getParameter<double>("dynktcut");
  rParam_ = iConfig.getParameter<double>("rParam");

  if (aggregateHF_) {
    withTruthInfo_ = iConfig.getParameter<bool>("aggregateWithTruthInfo");
    withCuts_ = iConfig.getParameter<bool>("aggregateWithCuts");
    withBDT_ = iConfig.getParameter<bool>("aggregateWithBDT");

    if (withBDT_) 
      model_path_ = iConfig.getParameter<edm::FileInPath>("model_path");
  }

  // Get labels
  ipTagInfoLabel_ = iConfig.getParameter<std::string>("ipTagInfoLabel");
  svTagInfoLabel_ = iConfig.getParameter<std::string>("svTagInfoLabel");
  
  // Get tokens
  // std::cout << "DEBUG: before consume" << std::endl;
  jetSrc_ = consumes<std::vector<T>>(iConfig.getParameter<edm::InputTag>("jetSrc"));
  constitSrc_ = mayConsume<std::vector<reco::PFCandidate>>(iConfig.getParameter<edm::InputTag>("constitSrc"));
  packedConstitSrc_ = mayConsume<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("constitSrc"));
  // std::cout << "DEBUG: constitSrc consumed" << std::endl;

  if (aggregateHF_) {
    candToGenParticleMapToken_ = mayConsume<reco::TrackToGenParticleMap>(iConfig.getParameter<edm::InputTag>("candToGenParticleMap"));
  }

  std::string alias = (iConfig.getParameter<edm::InputTag>("jetSrc")).label();
  produces<std::vector<reco::BasicJet>>().setBranchAlias(alias);
  produces<std::vector<reco::BasicJet>>("SubJets").setBranchAlias(alias);
  produces<std::vector<pat::PackedCandidate>>("pseudoHF");
}

template <class T>
void dynGroomedJets<T>::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  // std::cout << "In dynGroomedJets" << std::endl;

  auto jetCollection = std::make_unique<reco::BasicJetCollection>();
  auto subjetCollection = std::make_unique<reco::BasicJetCollection>();
  auto pseudoHFCollection = std::make_unique<std::vector<pat::PackedCandidate>>();

  // This will store the handle for the subjets after we write them
  edm::OrphanHandle<std::vector<reco::BasicJet>> subjetHandleAfterPut;
  // this is the mapping of subjet to hard jet
  std::vector<std::vector<int>> indices;
  // this is the list of hardjet 4-momenta
  std::vector<math::XYZTLorentzVector> p4_hardJets;
  // this is the hardjet areas
  std::vector<double> area_hardJets;
  std::vector<bool> isHardest;

  edm::Handle<std::vector<T>> jets;
  iEvent.getByToken(jetSrc_, jets);

  edm::Handle<std::vector<reco::PFCandidate>> pfcands;
  bool isPF = iEvent.getByToken(constitSrc_, pfcands);

  edm::Handle<edm::View<pat::PackedCandidate>> pfcandsPacked;
  bool isPackedPF = iEvent.getByToken(packedConstitSrc_, pfcandsPacked);
  // std::cout << "DEBUG: pfcandPacked ok" << std::endl;

  // -- For aggregation -- //
  edm::Handle<reco::TrackToGenParticleMap> candToGenParticleMap;
  BoosterHandle booster_;
  if (aggregateHF_) {
    iEvent.getByToken(candToGenParticleMapToken_, candToGenParticleMap);
    if (withBDT_) {
      // Following instructions here: https://cms-ml.github.io/documentation/inference/xgboost.html
      XGBoosterCreate(NULL, 0, &booster_);
      auto res = XGBoosterLoadModel(booster_, model_path_.fullPath().c_str());    
      // std::cout << "model loaded returned " << res << std::endl;
    }
  }
  // --------------------- //

  indices.resize(jets->size());

  int jetIndex = 0;
  for (const T& jet : *jets) { 
    // std::cout << "new jet with pt: " << jet.pt() << std::endl;

    p4_hardJets.push_back(math::XYZTLorentzVector(jet.px(), jet.py(), jet.pz(), jet.energy()));
    area_hardJets.push_back(jet.jetArea());

    std::vector<fastjet::PseudoJet> jetConstituents = {};

    fastjet::PseudoJet *subFJ1 = new fastjet::PseudoJet();
    fastjet::PseudoJet *subFJ2 = new fastjet::PseudoJet();
    std::vector<fastjet::PseudoJet> constit1;
    std::vector<fastjet::PseudoJet> constit2;

    // aggregate B based on jet type
    if (aggregateHF_) {
      pat::PackedCandidate pseudoHF;
      if (doGenJets_) {
		    // std::cout << "\n***Aggregating HF for gen jet" << std::endl;
        const reco::GenJet *genJet = jet.genJet();
        if (genJet) {
          auto tempPair = aggregateHFGen(*genJet, *candToGenParticleMap);
          jetConstituents = tempPair.first;
          pseudoHF = tempPair.second;
        }
        // std::cout << "genjet i = " << jetIndex << " pt " << jet.pt() << " mb " << pseudoHF.mass() << std::endl;
      } else {
		    // std::cout << "\n***Aggregating HF for reco jet" << std::endl;  
        auto tempPair = aggregateHFReco(jet, *candToGenParticleMap, booster_);
        jetConstituents = tempPair.first;
        pseudoHF = tempPair.second;
        // std::cout << "reco jet i = " << jetIndex << " pt " << jet.pt() << " mb " << pseudoHF.mass() << std::endl;
      } 

      pseudoHFCollection->push_back(pseudoHF);
    } else {
      // std::cout << "\tNot aggregating" << std::endl;
      std::vector<edm::Ptr<reco::Candidate>> constituents = {}; 
      if (doGenJets_) {
        const reco::GenJet *genJet = jet.genJet();
        if (genJet) constituents = genJet->getJetConstituents();
      } else {
        constituents = jet.getJetConstituents();
      }

      for (edm::Ptr<reco::Candidate> constituent : constituents) {
        if ((chargedOnly_) && (constituent->charge() == 0)) continue;
        if (constituent->pt() < ptCut_) continue;
        jetConstituents.push_back(fastjet::PseudoJet(constituent->px(), constituent->py(), constituent->pz(), constituent->energy()));
      }
    }


    // Iterative declustering
    isHardest.push_back(IterativeDeclustering(jetConstituents, subFJ1, subFJ2, constit1, constit2));

    // Convert fastjets to basicjets 
    reco::BasicJet subjet1, subjet2;
    if(isPF){
      subjet1 = ConvertFJ2BasicJet(subFJ1, constit1, pfcands, iSetup);
      subjet2 = ConvertFJ2BasicJet(subFJ2, constit2, pfcands, iSetup);
    }
    else{
      subjet1 = ConvertFJ2BasicJet(subFJ1, constit1, pfcandsPacked, iSetup);
      subjet2 = ConvertFJ2BasicJet(subFJ2, constit2, pfcandsPacked, iSetup);
    }

    if (subjet1.pt() > 1.0e-3) {
      indices[jetIndex].push_back(subjetCollection->size());
      if (subFJ1->has_area()) subjet1.setJetArea(subFJ1->area());
      subjetCollection->push_back(subjet1);
    }
    if (subjet2.pt() > 1.0e-3) {
      indices[jetIndex].push_back(subjetCollection->size());
      if (subFJ2->has_area()) subjet2.setJetArea(subFJ2->area());
      subjetCollection->push_back(subjet2);
    }

    jetIndex++;
  } // end jet loop

  // put subjets into event record
  subjetHandleAfterPut = iEvent.put(move(subjetCollection), "SubJets"); 

  // Now create the hard jets with ptr's to the subjets as constituents
  std::vector<math::XYZTLorentzVector>::const_iterator ip4 = p4_hardJets.begin(), ip4Begin = p4_hardJets.begin(), ip4End = p4_hardJets.end();
  
  for (; ip4 != ip4End; ++ip4) {
    int p4_index = ip4 - ip4Begin;
    std::vector<int>& ind = indices[p4_index];
    std::vector<edm::Ptr<reco::Candidate>> i_hardJetConstituents;

    // Add the subjets to the hard jet
    for (std::vector<int>::const_iterator isub = ind.begin(); isub != ind.end(); ++isub) {
      edm::Ptr<reco::Candidate> candPtr(subjetHandleAfterPut, *isub, false);
      i_hardJetConstituents.push_back(candPtr);
    }
    
    reco::Particle::Point point(0, 0, 0);
    //cout<<" size of i_hardJetConstituents "<<i_hardJetConstituents.size()<<endl;
    reco::BasicJet toput(*ip4, point, i_hardJetConstituents);
    if (isHardest[ip4-ip4Begin]) toput.setJetArea(0.8);
    else toput.setJetArea(0.4);
    jetCollection->push_back(toput);
  }

  //cout << "jetCollection size: " << jetCollection->size() << endl; 
  iEvent.put(move(jetCollection));
  iEvent.put(move(pseudoHFCollection), "pseudoHF");
}

template <class T>
reco::BasicJet dynGroomedJets<T>::ConvertFJ2BasicJet(fastjet::PseudoJet *fj, 
                                                     std::vector<fastjet::PseudoJet> constit, 
                                                     edm::Handle<std::vector<reco::PFCandidate>> pfcands, 
                                                     const edm::EventSetup& iSetup) const
{
  math::XYZTLorentzVector p4(fj->px(), fj->py(), fj->pz(), fj->e());  
  reco::Particle::Point point(0, 0, 0);
  std::vector<edm::Ptr<reco::Candidate>> constituents;
  if (writeConstits_) {
    for(uint j = 0; j < constit.size(); j++){
      double constitPt = constit[j].pt();
      double constitEta = constit[j].eta();
      
      int iCand = -1;
      for (const reco::PFCandidate& pfcand : *pfcands) {
        iCand++;
        
        if (std::fabs(constitPt - pfcand.pt()) < 0.0001 && std::fabs(constitEta - pfcand.eta()) < 0.0001 ){
          constituents.push_back(edm::Ptr<reco::Candidate>(pfcands, iCand));
          break;
        }
      }
    }
  }
  reco::BasicJet basicjet;
  writeSpecific(basicjet, p4, point, constituents, iSetup);

  return basicjet;
}

template <class T>
reco::BasicJet dynGroomedJets<T>::ConvertFJ2BasicJet(fastjet::PseudoJet *fj, 
                                                     std::vector<fastjet::PseudoJet> constit, 
                                                     edm::Handle<edm::View<pat::PackedCandidate>> pfcands, 
                                                     const edm::EventSetup& iSetup) const
{
  math::XYZTLorentzVector p4(fj->px(), fj->py(), fj->pz(), fj->e());  
  reco::Particle::Point point(0, 0, 0);
  std::vector<edm::Ptr<reco::Candidate>> constituents;
  if (writeConstits_) {
    for(uint j = 0; j < constit.size(); j++){
      double constitPt = constit[j].pt();
      double constitEta = constit[j].eta();
      
      int iCand = -1;
      for (const pat::PackedCandidate& pfcand : *pfcands) {
        iCand++;
        
        if (std::fabs(constitPt - pfcand.pt()) < 0.0001 && std::fabs(constitEta - pfcand.eta()) < 0.0001 ){
          constituents.push_back(edm::Ptr<reco::Candidate>(pfcands, iCand));
          break;
        }
      }
    }
  }
  reco::BasicJet basicjet;
  writeSpecific(basicjet, p4, point, constituents, iSetup);

  return basicjet;
}

template <class T>
bool dynGroomedJets<T>::IterativeDeclustering(std::vector<fastjet::PseudoJet> jetConstituents, 
                                              fastjet::PseudoJet *sub1, fastjet::PseudoJet *sub2, 
                                              std::vector<fastjet::PseudoJet> &constit1, std::vector<fastjet::PseudoJet> &constit2) const
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
  fastjet::JetDefinition jet_def(fastjet::genkt_algorithm, jet_radius_ca, 0, static_cast <fastjet::RecombinationScheme>(0), fastjet::Best);

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
      
      while (jj.has_parents(j1, j2)) {
	  
        if (j1.perp() < j2.perp()) std::swap(j1,j2);
        
        double delta_R = j1.delta_R(j2);
        double cut = zcut_ * pow(delta_R / rParam_, beta_);
        double z = j2.perp() / (j1.perp() + j2.perp());
        if (z > cut && (doLateSD_ || !flagSubjet) ) {
          flagSubjet = true;
          j1first = j1;
          j2first = j2;
          *sub1 = j1first;
          *sub2 = j2first;
          nsdin = nsplit;
        }
        double dyn= z * (1-z) * j2.perp() * pow(delta_R / rParam_, dynktcut_);
        
        if (dyn > kt1) {
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
    } catch (fastjet::Error) {} //return -1; }
  return isHardest;
}

// Function to aggregate reconstructed tracks into pseudo-B's 
template <class T>
typename dynGroomedJets<T>::jetConstituentsPseudoHFPair dynGroomedJets<T>::aggregateHFReco(const T& jet, 
                                                               reco::TrackToGenParticleMap candToGenParticleMap, 
                                                               BoosterHandle booster_) const
{
  // std::cout << "Aggregating Bs in reco jet" << std::endl;

  // Input and output particle collections
  std::vector<edm::Ptr<reco::Candidate>> inputJetConstituents = jet.getJetConstituents();
  std::vector<fastjet::PseudoJet> outputJetConstituents = {};
  pat::PackedCandidate outputPseudoHF;

  // Particle collection to aggregate into pseudo-Bs
  std::map<int, std::vector<edm::Ptr<reco::Candidate>>> hfConstituentsMap;
  reco::Candidate::PolarLorentzVector totalPseudoHF(0., 0., 0., 0.);

  // Grab the IP and SV tag info from the jet
  const reco::CandIPTagInfo *ipTagInfo = jet.tagInfoCandIP(ipTagInfoLabel_.c_str());
  const std::vector<reco::btag::TrackIPData> ipData = ipTagInfo->impactParameterData();
  const std::vector<edm::Ptr<reco::Candidate>> ipTracks = ipTagInfo->selectedTracks();

  const reco::CandSecondaryVertexTagInfo *svTagInfo = jet.tagInfoCandSecondaryVertex(svTagInfoLabel_.c_str());

  // Go over jet constituents 
  for (const edm::Ptr<reco::Candidate> constit : jet.getJetConstituents()) {
    if (chargedOnly_ && constit->charge() == 0) continue;
    if (constit->pt() < ptCut_) continue;
    
    // Look for particle in ipTracks
    auto itIPTrack = std::find(ipTracks.begin(), ipTracks.end(), constit);
    if (itIPTrack == ipTracks.end()) continue;

    // [DEBUG] 
    // std::cout << "constit pt " << constit->pt() << " itIPTrack pt " << (*itIPTrack)->pt() << std::endl;

    // Get status of particle (1 = not from HF, 100+ = from HF) 
    int status = 1;

    if (withTruthInfo_) {
      if (candToGenParticleMap.find(constit) != candToGenParticleMap.end()) { 
        edm::Ptr<pat::PackedGenParticle> matchGenParticle = candToGenParticleMap.at(constit);
        status = matchGenParticle->status();
      }
    } else {
      // Initialize values 
      const double missing_value = -1000000.;

      float ip3dSig = missing_value;
      bool inSV = false;
      float svtxdls = missing_value;
      float svtxm = missing_value;

      // Get IP info 
      int trkIPIndex = itIPTrack - ipTracks.begin();
      const reco::btag::TrackIPData trkIPdata = ipData[trkIPIndex];
      ip3dSig = trkIPdata.ip3d.significance();
      // if nan go back to missing_value
      if (ip3dSig != ip3dSig) ip3dSig = missing_value;
         
      // Get SV info
      for (uint ivtx = 0; ivtx < svTagInfo->nVertices(); ivtx++) {
        std::vector<edm::Ptr<reco::Candidate>> isvTracks = svTagInfo->vertexTracks(ivtx);
        auto itSVTrack = std::find(isvTracks.begin(), isvTracks.end(), constit);
        if (itSVTrack == isvTracks.end()) continue;

        inSV = true;
        
        Measurement1D m1D = svTagInfo->flightDistance(ivtx);
        svtxdls = m1D.significance();

        Measurement1D m2D = svTagInfo->flightDistance(ivtx, true);
        
        const reco::VertexCompositePtrCandidate svtx = svTagInfo->secondaryVertex(ivtx);
        svtxm = svtx.p4().mass();

        break;
      } // end vtx loop

      if (withCuts_) {
        if (inSV || (ip3dSig > 3.)) {
          status = 100;
        } 
      } else if (withBDT_) {
        float inSVBDT = (inSV) ? 1. : 0.;

        // Initialize BDT related variables
        bst_ulong out_len = 0; // bst_ulong is a typedef of unsigned long
        const float *f; // array to store predictions
        float threshold = 0.44;

        DMatrixHandle data_;
        const int nFeatures = 4;
        const int nEntries = 1;
        float trackDataBDT[nEntries][nFeatures] = {{inSVBDT, svtxdls, svtxm, ip3dSig}};
        
        XGDMatrixCreateFromMat((float *)trackDataBDT, nEntries, nFeatures, missing_value, &data_);
        XGBoosterPredict(booster_, data_, 0, 0, &out_len, &f);
        float prediction = f[0];
        // std::cout << "trk with pt " << constit->pt()
        //           << ", ip3dsig " << ip3dSig
        //           << ", in SV " << inSVBDT
        //           << ", svtxm " << svtxm
        //           << ", svtxdls " << svtxdls
        //           << ", model prediction " << prediction 
        //           << std::endl;

        if (prediction > threshold) {
          status = 100;
        }
      }
    }

    // Add particle to output collection or from HF map
    if (status == 1) {
      fastjet::PseudoJet outConstit(constit->px(), constit->py(), constit->pz(), constit->energy());
      outputJetConstituents.push_back(outConstit);
    } else if (status >= 100) {
      hfConstituentsMap[status].push_back(constit);
    }
  } // end jet constituents loop
  // Aggregate particles coming from HF decays into pseudo-B/C's 
  // and add them to the collection
  
  for (auto itTrackFromHF = hfConstituentsMap.begin(); itTrackFromHF != hfConstituentsMap.end(); itTrackFromHF++) {
    reco::Candidate::PolarLorentzVector pseudoHF(0., 0., 0., 0.);
    for (edm::Ptr<reco::Candidate> hfConstituent : itTrackFromHF->second) {
      reco::Candidate::PolarLorentzVector productLorentzVector(0., 0., 0., 0.);
      productLorentzVector.SetPt(hfConstituent->pt());
      productLorentzVector.SetEta(hfConstituent->eta());
      productLorentzVector.SetPhi(hfConstituent->phi());
      productLorentzVector.SetM(hfConstituent->mass());
      pseudoHF += productLorentzVector;
      totalPseudoHF += productLorentzVector;
    }
    fastjet::PseudoJet outPseudoHFConstit(pseudoHF.Px(), pseudoHF.Py(), pseudoHF.Pz(), pseudoHF.E());
    outputJetConstituents.push_back(outPseudoHFConstit);
  } // end tracks from B loop
  
  outputPseudoHF.setP4(totalPseudoHF);
  return dynGroomedJets<T>::jetConstituentsPseudoHFPair {outputJetConstituents, outputPseudoHF};
} // end aggregateHFReco()

// Function to aggregate gen particles coming from B decays into pseudo-B's 
template <class T>
typename dynGroomedJets<T>::jetConstituentsPseudoHFPair dynGroomedJets<T>::aggregateHFGen(const T& genJet,
                                                                                          reco::TrackToGenParticleMap candToGenParticleMap) const
{
  // std::cout << "Aggregating HF at gen level" << std::endl;

  // Input and output particle collections
  std::vector<edm::Ptr<reco::Candidate>> inputJetConstituents = genJet.getJetConstituents();
  std::vector<fastjet::PseudoJet> outputJetConstituents = {};
  pat::PackedCandidate outputPseudoHF;

  // Particle collection to aggregate into pseudo-Bs
  std::map<int, std::vector<edm::Ptr<reco::Candidate>>> hfConstituentsMap;
  reco::Candidate::PolarLorentzVector totalPseudoHF(0., 0., 0., 0.);

  // Go over jet constituents
  for (edm::Ptr<reco::Candidate> constit : inputJetConstituents) {
    if(chargedOnly_ && constit->charge() == 0) continue;
    if(constit->pt() < ptCut_) continue;
  
    // Get status of matched gen particle
    edm::Ptr<pat::PackedGenParticle> matchGenParticle = candToGenParticleMap.at(constit);
    int status = matchGenParticle->status();

    // Add particle to output collection or from HF map
    if (status == 1) {
      fastjet::PseudoJet outConstit(constit->px(), constit->py(), constit->pz(), constit->energy());
      outputJetConstituents.push_back(outConstit);
    } else if (status >= 100) {
      hfConstituentsMap[status].push_back(constit);
    }
  } // end constit loop 

  // Aggregate particles coming from HF decays into pseudo-B/D's 
  // and add them to the collection
  for (auto it = hfConstituentsMap.begin(); it != hfConstituentsMap.end(); ++it) {
    // std::cout << "B with code " << it->first << " has " << (it->second).size() << " constituents" << std::endl;
    reco::Candidate::PolarLorentzVector pseudoHF(0., 0., 0., 0.);
    for (edm::Ptr<reco::Candidate> hfConstituent : it->second) {
      reco::Candidate::PolarLorentzVector productLorentzVector(0., 0., 0., 0.);
      productLorentzVector.SetPt(hfConstituent->pt());
      productLorentzVector.SetEta(hfConstituent->eta());
      productLorentzVector.SetPhi(hfConstituent->phi());
      productLorentzVector.SetM(hfConstituent->mass());
      pseudoHF += productLorentzVector;
      totalPseudoHF += productLorentzVector;
      //std::cout << "hf product lorentz vector: " << productLorentzVector << std::endl;
    }
    // std::cout << "added HF to collection" << std::endl;
    fastjet::PseudoJet outPseudoHFConstit(pseudoHF.Px(), pseudoHF.Py(), pseudoHF.Pz(), pseudoHF.E());
    outputJetConstituents.push_back(outPseudoHFConstit);
  } // end map loop
  
  // std::cout << "genjet has " << hfConstituentsMap.size() << " B's " << std::endl;

  outputPseudoHF.setP4(totalPseudoHF);
  // std::cout << "vector mass " << totalPseudoHF.M() << " candidate mass " << outputPseudoHF.mass() << std::endl;
  return dynGroomedJets<T>::jetConstituentsPseudoHFPair {outputJetConstituents, outputPseudoHF};
}

template <class T>
void dynGroomedJets<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setComment("Dynamically groomed jets");

  // Input collections
  desc.add<edm::InputTag>("jetSrc", edm::InputTag("slimmedJets"));
  desc.add<edm::InputTag>("constitSrc", edm::InputTag("packedPFCandidates"));
  desc.add<edm::InputTag>("candToGenParticleMap", edm::InputTag("TrackToGenParticleMapProducer"));

  // Configuration parameters
  desc.add<bool>("writeConstits", false);
  desc.add<bool>("doLateSD", false);
  desc.add<bool>("chargedOnly", false);
  
  desc.add<double>("zcut", 0.1);
  desc.add<double>("beta", 0.0);
  desc.add<double>("dynktcut", 1.0);
  desc.add<double>("rParam", 0.4);
  desc.add<double>("ptCut", 1.);

  desc.add<bool>("doGenJets", false);

  desc.add<bool>("aggregateHF", false);
  desc.add<bool>("aggregateWithTruthInfo", true);
  desc.add<bool>("aggregateWithCuts", false);
  desc.add<bool>("aggregateWithBDT", false);
  desc.add<edm::FileInPath>("model_path", edm::FileInPath("RecoHI/HiJetAlgos/data/dummy.model"));

  // Tag info labels
  desc.add<std::string>("ipTagInfoLabel", "pfImpactParameter");
  desc.add<std::string>("svTagInfoLabel", "pfInclusiveSecondaryVertexFinder");

  descriptions.add("dynGroomedPATJets", desc);
}

using dynGroomedPATJets = dynGroomedJets<pat::Jet>;

// define this as a plug-in
DEFINE_FWK_MODULE(dynGroomedPATJets);
