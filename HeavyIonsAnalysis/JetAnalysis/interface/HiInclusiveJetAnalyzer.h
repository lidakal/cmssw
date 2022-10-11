#ifndef MNguyen_HiInclusiveJetAnalyzer_inclusiveJetAnalyzer_
#define MNguyen_HiInclusiveJetAnalyzer_inclusiveJetAnalyzer_

// system include files
#include <memory>
#include <string>
#include <iostream>

// ROOT headers
#include "TTree.h"

// user include files
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfo.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

/**\class HiInclusiveJetAnalyzer

   \author Matt Nguyen
   \date   November 2010
*/

class HiInclusiveJetAnalyzer : public edm::EDAnalyzer {
public:
  explicit HiInclusiveJetAnalyzer(const edm::ParameterSet&);
  ~HiInclusiveJetAnalyzer();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginRun(const edm::Run & r, const edm::EventSetup & c);
  virtual void beginJob();

private:
  int  getGroomedJetIndex(const reco::Jet jet) const;
  int  getGroomedGenJetIndex(const reco::GenJet jet) const;
  int  getGroomedPartonJetIndex(const reco::GenJet jet) const;
  bool isHardProcess(const int);
  bool isFromGSP(const reco::Candidate* c);
  int trkGenPartMatch(const reco::CandidatePtr track, reco::GenParticleCollection genParticles, float ptCut);
  int trkInVector(reco::CandidatePtr trk, std::vector<reco::CandidatePtr> tracks) const;


  edm::InputTag   jetTagLabel_;
  edm::EDGetTokenT<reco::JetView>                jetTag_;
  edm::EDGetTokenT<pat::JetCollection>           jetTagPat_;
  edm::EDGetTokenT<pat::JetCollection>           subJetTagPat_;
  edm::EDGetTokenT<reco::JetView>                matchTag_;
  edm::EDGetTokenT<reco::PFCandidateCollection>  pfCandidateLabel_;
  edm::EDGetTokenT<std::vector<reco::GenJet>>    genjetTag_;
  edm::EDGetTokenT<std::vector<reco::GenJet>>    partonjetTag_;
  edm::EDGetTokenT<reco::JetView>                subjetGenTag_;

  std::string jetName_; //used as prefix for jet structures
  edm::Handle<reco::JetView> gensubjets_;
  edm::Handle<edm::View<reco::Jet> > groomedJets; 
  edm::Handle<edm::View<reco::Jet> > groomedGenJets; 
  edm::Handle<edm::View<reco::Jet> > groomedPartonJets; 

  bool doMatch_;
  bool useJEC_;
  bool isMC_;
  bool fillGenJets_;
  bool useQuality_;
  std::string trackQuality_;

  double genPtMin_;
  bool doLifeTimeTagging_;
  bool doLifeTimeTaggingExtra_;
  bool skipCorrections_;
  bool doHiJetID_;
  bool doStandardJetID_;

  double rParam;
  double hardPtMin_;
  double jetPtMin_;
  bool doSubJets_;

  bool doTrackMatching_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticleSrc_;

  TTree *t;
  edm::Service<TFileService> fs1;

  std::string bTagJetName_;
  std::string svTagInfos_;
  std::string ipTagInfos_;
  std::string jetPBJetTags_;
  std::string combinedSVV2BJetTags_;
  std::string deepCSVBJetTags_;
  std::string deepCSVBBJetTags_;
  std::string deepCSVCJetTags_;
  std::string deepFlavourBJetTags_;
  std::string deepFlavourBBJetTags_;
  std::string deepFlavourLEPBJetTags_;
  std::string deepFlavourCJetTags_;
  std::string deepCSVBSubJetTags_;
  std::string deepCSVBBSubJetTags_;
  std::string deepCSVCSubJetTags_;

  bool doExtendedFlavorTagging_;
  edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> jetFlavourInfosToken_;
  edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> subjetFlavourInfosToken_;
  edm::EDGetTokenT<edm::View<reco::Jet> >                  groomedJetsToken_;
  edm::EDGetTokenT<edm::View<reco::Jet> >                  groomedGenJetsToken_;
  edm::EDGetTokenT<edm::View<reco::Jet> >                  groomedPartonJetsToken_;
  bool                                                     useSubjets_;

  static const int MAXJETS = 100;
  static const int MAXTRACKS = 1000;
  //static const int MAXSVS = MAXJETS*3;
  //static const int MAXTRACKSINSVS = MAXTRACKS*3;

  struct JRA{

    int nref;

    float rawpt[MAXJETS];
    float jtpt[MAXJETS];
    float jteta[MAXJETS];
    float jtphi[MAXJETS];


    float jty[MAXJETS];
    float jtpu[MAXJETS];
    float jtm[MAXJETS];
    float jtarea[MAXJETS];

    float jtPfCHF[MAXJETS];
    float jtPfNHF[MAXJETS];
    float jtPfCEF[MAXJETS];
    float jtPfNEF[MAXJETS];
    float jtPfMUF[MAXJETS];

    int jtPfCHM[MAXJETS];
    int jtPfNHM[MAXJETS];
    int jtPfCEM[MAXJETS];
    int jtPfNEM[MAXJETS];
    int jtPfMUM[MAXJETS];
    
    float jttau1[MAXJETS];
    float jttau2[MAXJETS];
    float jttau3[MAXJETS];

    float jtHadFlav[MAXJETS];
    float jtParFlav[MAXJETS];
    int jtNbHad[MAXJETS];
    int jtNcHad[MAXJETS];
    int jtNbPar[MAXJETS];
    int jtNcPar[MAXJETS];
    bool jtHasGSPB[MAXJETS];
    bool jtHasGSPC[MAXJETS];

    bool jtIsHardest[MAXJETS];
    float jtptG[MAXJETS];
    float jtetaG[MAXJETS];
    float jtphiG[MAXJETS];
    float jtmassG[MAXJETS];

    float sjt1E[MAXJETS];
    float sjt1Pt[MAXJETS];
    float sjt1Eta[MAXJETS];
    float sjt1Phi[MAXJETS];
    float sjt1Mass[MAXJETS];
    float sjt1HadFlav[MAXJETS];
    float sjt1ParFlav[MAXJETS];
    float sjt1DiscDeepCSVB[MAXJETS];
    float sjt1DiscDeepCSVBB[MAXJETS];
    float sjt1DiscDeepCSVC[MAXJETS];

    float sjt2E[MAXJETS];
    float sjt2Pt[MAXJETS];
    float sjt2Eta[MAXJETS];
    float sjt2Phi[MAXJETS];
    float sjt2Mass[MAXJETS];
    float sjt2HadFlav[MAXJETS];
    float sjt2ParFlav[MAXJETS];
    float sjt2DiscDeepCSVB[MAXJETS];
    float sjt2DiscDeepCSVBB[MAXJETS];
    float sjt2DiscDeepCSVC[MAXJETS];


    float gsjt1E[MAXJETS];
    float gsjt1Pt[MAXJETS];
    float gsjt1Eta[MAXJETS];
    float gsjt1Phi[MAXJETS];
    float gsjt1Mass[MAXJETS];
    float gsjt2E[MAXJETS];
    float gsjt2Pt[MAXJETS];
    float gsjt2Eta[MAXJETS];
    float gsjt2Phi[MAXJETS];
    float gsjt2Mass[MAXJETS];
    float rsjt1E[MAXJETS];
    float rsjt1Pt[MAXJETS];
    float rsjt1Eta[MAXJETS];
    float rsjt1Phi[MAXJETS];
    float rsjt1Mass[MAXJETS];
    float rsjt2E[MAXJETS];
    float rsjt2Pt[MAXJETS];
    float rsjt2Eta[MAXJETS];
    float rsjt2Phi[MAXJETS];
    float rsjt2Mass[MAXJETS];

    float chargedMax[MAXJETS];
    float chargedSum[MAXJETS];
    int chargedN[MAXJETS];

    float photonMax[MAXJETS];
    float photonSum[MAXJETS];
    int photonN[MAXJETS];

    float chargedHardSum[MAXJETS];
    float photonHardSum[MAXJETS];

    int trackHardN[MAXJETS];
    int chargedHardN[MAXJETS];
    int photonHardN[MAXJETS];

    float neutralMax[MAXJETS];
    float neutralSum[MAXJETS];
    int neutralN[MAXJETS];

    float eMax[MAXJETS];
    float eSum[MAXJETS];
    int eN[MAXJETS];

    float muMax[MAXJETS];
    float muSum[MAXJETS];
    int muN[MAXJETS];

    float hcalSum[MAXJETS];
    float ecalSum[MAXJETS];

    float fHPD[MAXJETS];
    float fRBX[MAXJETS];
    int n90[MAXJETS];
    float fSubDet1[MAXJETS];
    float fSubDet2[MAXJETS];
    float fSubDet3[MAXJETS];
    float fSubDet4[MAXJETS];
    float restrictedEMF[MAXJETS];
    int nHCAL[MAXJETS];
    int nECAL[MAXJETS];
    float apprHPD[MAXJETS];
    float apprRBX[MAXJETS];

    int n2RPC[MAXJETS];
    int n3RPC[MAXJETS];
    int nRPC[MAXJETS];

    float fEB[MAXJETS];
    float fEE[MAXJETS];
    float fHB[MAXJETS];
    float fHE[MAXJETS];
    float fHO[MAXJETS];
    float fLong[MAXJETS];
    float fShort[MAXJETS];
    float fLS[MAXJETS];
    float fHFOOT[MAXJETS];

    float matchedPt[MAXJETS];
    float matchedRawPt[MAXJETS];
    float matchedR[MAXJETS];
    float matchedPu[MAXJETS];

    float jtDiscCSVV2[MAXJETS];
    float jtDiscDeepCSVB[MAXJETS];
    float jtDiscDeepCSVBB[MAXJETS];
    float jtDiscDeepCSVC[MAXJETS];
    float jtDiscDeepFlavourB[MAXJETS];
    float jtDiscDeepFlavourBB[MAXJETS];
    float jtDiscDeepFlavourC[MAXJETS];
    float jtDiscDeepFlavourLEPB[MAXJETS];
    float jtDiscProb[MAXJETS];
 

    int nsvtx[MAXJETS];
    
    std::vector<std::vector<int> >svtxntrk;
    std::vector<std::vector<float> >svtxdl;
    std::vector<std::vector<float> >svtxdls;
    std::vector<std::vector<float> >svtxdl2d;
    std::vector<std::vector<float> >svtxdls2d;
    std::vector<std::vector<float> >svtxm;
    std::vector<std::vector<float> >svtxmcorr;
    std::vector<std::vector<float> >svtxpt;
    std::vector<std::vector<float> >svJetDeltaR;

    std::vector<std::vector<float> >svtxTrPt;
    std::vector<std::vector<float> >svtxTrEta;
    std::vector<std::vector<float> >svtxTrPhi;
    
    /*
    float svtxntrk[MAXSVS];
    float svtxdl[MAXSVS];
    float svtxdls[MAXSVS];
    float svtxdl2d[MAXSVS];
    float svtxdls2d[MAXSVS];
    float svtxm[MAXSVS];
    float svtxmcorr[MAXSVS];
    float svtxpt[MAXSVS];
    float svJetDeltaR[MAXSVS];

    float svtxTrPt[MAXTRACKSINSVS];
    float svtxTrEta[MAXTRACKSINSVS];
    float svtxTrPhi[MAXTRACKSINSVS];
    */

    float mue[MAXJETS];
    float mupt[MAXJETS];
    float mueta[MAXJETS];
    float muphi[MAXJETS];
    float mudr[MAXJETS];
    float muptrel[MAXJETS];
    int muchg[MAXJETS];


    float refpt[MAXJETS];
    float refeta[MAXJETS];
    float refphi[MAXJETS];
    float refm[MAXJETS];
    float refarea[MAXJETS];
    float refy[MAXJETS];
    float refdphijt[MAXJETS];
    float refdrjt[MAXJETS];
    float refparton_pt[MAXJETS];
    int refparton_flavor[MAXJETS];

    bool refIsHardest[MAXJETS];
    float refptG[MAXJETS];
    float refetaG[MAXJETS];
    float refphiG[MAXJETS];
    float refmG[MAXJETS];

    float pthat;
    int ngen;
    int genmatchindex[MAXJETS];
    float genpt[MAXJETS];
    float geneta[MAXJETS];
    float genphi[MAXJETS];
    float genm[MAXJETS];
    float gendphijt[MAXJETS];
    float gendrjt[MAXJETS];

    bool genIsHardest[MAXJETS];
    float genptG[MAXJETS];
    float genetaG[MAXJETS];
    float genphiG[MAXJETS];
    float genmG[MAXJETS];

    int npar;
    float parpt[MAXJETS];
    float pareta[MAXJETS];
    float parphi[MAXJETS];
    float parm[MAXJETS];

    int parNb[MAXJETS];
    int parNc[MAXJETS];
    bool parHasGSPB[MAXJETS];
    bool parHasGSPC[MAXJETS];

    float parptG[MAXJETS];
    float paretaG[MAXJETS];
    float parphiG[MAXJETS];
    float parmG[MAXJETS];

    float psjt1E[MAXJETS];
    float psjt1Pt[MAXJETS];
    float psjt1Eta[MAXJETS];
    float psjt1Phi[MAXJETS];
    float psjt1Mass[MAXJETS];
    float psjt2E[MAXJETS];
    float psjt2Pt[MAXJETS];
    float psjt2Eta[MAXJETS];
    float psjt2Phi[MAXJETS];
    float psjt2Mass[MAXJETS];

    int nselIPtrk[MAXJETS];
    
    int nIP;
    float ipPt[MAXTRACKS];
    float ipEta[MAXTRACKS];
    float ipPhi[MAXTRACKS];
    float ipProb[MAXTRACKS];
    float ip3dSig[MAXTRACKS];

    float ipPtMatch[MAXTRACKS];
    float ipEtaMatch[MAXTRACKS];
    float ipPhiMatch[MAXTRACKS];
    int ipMatchStatus[MAXTRACKS];

    int ipInSV[MAXTRACKS];
    float ipSvtxdls[MAXTRACKS];
    float ipSvtxdls2d[MAXTRACKS];
    float ipSvtxm[MAXTRACKS];
    float ipSvtxmcorr[MAXTRACKS];

  };

  JRA jets_;

};

#endif
