/*
  Based on the jet response analyzer
  Modified by Matt Nguyen, November 2010

*/

#include "HeavyIonsAnalysis/JetAnalysis/interface/HiInclusiveJetAnalyzer.h"

#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "TMath.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

using namespace std;
using namespace edm;
using namespace reco;
using namespace pat;

HiInclusiveJetAnalyzer::HiInclusiveJetAnalyzer(const edm::ParameterSet& iConfig){

  doMatch_ = iConfig.getUntrackedParameter<bool>("matchJets",false);
  jetTagLabel_ = iConfig.getParameter<InputTag>("jetTag");
  jetTagPat_ = consumes<JetCollection> (jetTagLabel_);
  subJetTagPat_ = consumes<JetCollection> (iConfig.getParameter<InputTag>("subJetTag"));
  matchTag_ = consumes<JetView> (iConfig.getUntrackedParameter<InputTag>("matchTag"));
  jetName_ = iConfig.getUntrackedParameter<string>("jetName");
  doSubJets_ = iConfig.getUntrackedParameter<bool>("doSubJets",0);
  subjetGenTag_ = consumes<JetView> (iConfig.getUntrackedParameter<InputTag>("subjetGenTag"));
  isMC_ = iConfig.getUntrackedParameter<bool>("isMC",false);
  fillGenJets_ = iConfig.getUntrackedParameter<bool>("fillGenJets",false);
  doExtendedFlavorTagging_ = iConfig.getUntrackedParameter<bool>("doExtendedFlavorTagging",true);
  doHiJetID_ = iConfig.getUntrackedParameter<bool>("doHiJetID",false);
  doStandardJetID_ = iConfig.getUntrackedParameter<bool>("doStandardJetID",false);
  rParam = iConfig.getParameter<double>("rParam");
  hardPtMin_ = iConfig.getUntrackedParameter<double>("hardPtMin",4);
  jetPtMin_ = iConfig.getParameter<double>("jetPtMin");
  if(isMC_) genjetTag_ = consumes<vector<GenJet>>(iConfig.getParameter<InputTag>("genjetTag"));
  if(isMC_) partonjetTag_ = consumes<vector<GenJet>>(iConfig.getParameter<InputTag>("partonjetTag"));
  useJEC_ = iConfig.getUntrackedParameter<bool>("useJEC",true);
  doLifeTimeTagging_ = iConfig.getUntrackedParameter<bool>("doLifeTimeTagging",false);
  doLifeTimeTaggingExtra_ = iConfig.getUntrackedParameter<bool>("doLifeTimeTaggingExtra",false);
  doTrackMatching_ = iConfig.getUntrackedParameter<bool>("doTrackMatching", false);
  genParticleSrc_ = mayConsume<reco::GenParticleCollection>(iConfig.getParameter<InputTag>("genParticlesTag"));
  skipCorrections_  = iConfig.getUntrackedParameter<bool>("skipCorrections",false);
  pfCandidateLabel_ = consumes<PFCandidateCollection> (iConfig.getUntrackedParameter<edm::InputTag>("pfCandidateLabel",edm::InputTag("particleFlow")));  
  if (doSubJets_){ 
    groomedJetsToken_ = consumes<edm::View<reco::Jet> >(iConfig.getParameter<edm::InputTag>("groomedJets"));
    if(isMC_)groomedGenJetsToken_ = consumes<edm::View<reco::Jet> >(iConfig.getParameter<edm::InputTag>("groomedGenJets"));
    if(isMC_)groomedPartonJetsToken_ = consumes<edm::View<reco::Jet> >(iConfig.getParameter<edm::InputTag>("groomedPartonJets"));
  }
  useSubjets_ = iConfig.exists("subjetFlavourInfos") && iConfig.exists("groomedJets");
  if (isMC_&& doExtendedFlavorTagging_) {
    jetFlavourInfosToken_ = consumes<JetFlavourInfoMatchingCollection>( iConfig.getParameter<edm::InputTag>("jetFlavourInfos") );
    //if (useSubjets_) subjetFlavourInfosToken_ = consumes<JetFlavourInfoMatchingCollection>(iConfig.getParameter<edm::InputTag>("subjetFlavourInfos"));
  }

  if(doLifeTimeTagging_){
    bTagJetName_ = iConfig.getUntrackedParameter<string>("bTagJetName");
    if(doLifeTimeTaggingExtra_){
      svTagInfos_ = bTagJetName_+"PfInclusiveSecondaryVertexFinder";
      ipTagInfos_ = bTagJetName_+"PfImpactParameter";
    }
    jetPBJetTags_ = bTagJetName_+"JetProbabilityBJetTags";
    combinedSVV2BJetTags_ = bTagJetName_+"CombinedSecondaryVertexV2BJetTags";
    deepCSVBJetTags_ = bTagJetName_+"PfDeepCSVJetTags:probb";
    deepCSVBBJetTags_ = bTagJetName_+"PfDeepCSVJetTags:probbb";
    deepCSVCJetTags_ = bTagJetName_+"PfDeepCSVJetTags:probc";
    deepFlavourBJetTags_ = bTagJetName_+"PfDeepFlavourJetTags:probb";
    deepFlavourBBJetTags_ = bTagJetName_+"PfDeepFlavourJetTags:probbb";
    deepFlavourLEPBJetTags_ = bTagJetName_+"PfDeepFlavourJetTags:problepb";
    deepFlavourCJetTags_ = bTagJetName_+"PfDeepFlavourJetTags:probc";
    
    deepCSVBSubJetTags_ = bTagJetName_+"PfDeepCSVSubJetTags:probb";
    deepCSVBBSubJetTags_ = bTagJetName_+"PfDeepCSVSubJetTags:probbb";
    deepCSVCSubJetTags_ = bTagJetName_+"PfDeepCSVSubJetTags:probc";
  }


  if(isMC_)genPtMin_ = iConfig.getUntrackedParameter<double>("genPtMin",10);

}



HiInclusiveJetAnalyzer::~HiInclusiveJetAnalyzer() { }



void
HiInclusiveJetAnalyzer::beginRun(const edm::Run& run,
				 const edm::EventSetup & es) {}

void
HiInclusiveJetAnalyzer::beginJob() {

  //string jetTagName = jetTag_.label()+"_tree";
  string jetTagTitle = jetTagLabel_.label()+" Jet Analysis Tree";
  t = fs1->make<TTree>("t",jetTagTitle.c_str());

  t->Branch("nref",&jets_.nref,"nref/I");
  t->Branch("rawpt",jets_.rawpt,"rawpt[nref]/F");
  if(!skipCorrections_) t->Branch("jtpt",jets_.jtpt,"jtpt[nref]/F");
  t->Branch("jteta",jets_.jteta,"jteta[nref]/F");
  t->Branch("jtphi",jets_.jtphi,"jtphi[nref]/F");

  if(!doSubJets_){
    t->Branch("jtPfCHF",jets_.jtPfCHF,"jtPfCHF[nref]/F");
    t->Branch("jtPfNHF",jets_.jtPfNHF,"jtPfNHF[nref]/F");
    t->Branch("jtPfCEF",jets_.jtPfCEF,"jtPfCEF[nref]/F");
    t->Branch("jtPfNEF",jets_.jtPfNEF,"jtPfNEF[nref]/F");
    t->Branch("jtPfMUF",jets_.jtPfMUF,"jtPfMUF[nref]/F");
    
    t->Branch("jtPfCHM",jets_.jtPfCHM,"jtPfCHM[nref]/I");
    t->Branch("jtPfNHM",jets_.jtPfNHM,"jtPfNHM[nref]/I");
    t->Branch("jtPfCEM",jets_.jtPfCEM,"jtPfCEM[nref]/I");
    t->Branch("jtPfNEM",jets_.jtPfNEM,"jtPfNEM[nref]/I");
    t->Branch("jtPfMUM",jets_.jtPfMUM,"jtPfMUM[nref]/I");
  }
  if(doExtendedFlavorTagging_){
    t->Branch("jtHadFlav",jets_.jtHadFlav,"jtHadFlav[nref]/F");
    t->Branch("jtParFlav",jets_.jtParFlav,"jtParFlav[nref]/F");
    t->Branch("jtNbHad",jets_.jtNbHad,"jtNbHad[nref]/I");
    t->Branch("jtNcHad",jets_.jtNcHad,"jtNcHad[nref]/I");
    t->Branch("jtNbPar",jets_.jtNbPar,"jtNbPar[nref]/I");
    t->Branch("jtNcPar",jets_.jtNcPar,"jtNcPar[nref]/I");
    t->Branch("jtHasGSPB",jets_.jtHasGSPB,"jtHasGSPB[nref]/O");
    t->Branch("jtHasGSPC",jets_.jtHasGSPC,"jtHasGSPC[nref]/O");
  } 
 
  if(doSubJets_) {
    /*
    t->Branch("jtptG",jets_.jtptG,"jtptG[nref]/F");
    t->Branch("jtetaG",jets_.jtetaG,"jtetaG[nref]/F");
    t->Branch("jtphiG",jets_.jtphiG,"jtphiG[nref]/F");
    */
    t->Branch("jtIsHardest",jets_.jtIsHardest,"jtIsHardest[nref]/O");
    t->Branch("sjt1E",jets_.sjt1E,"sjt1E[nref]/F");
    t->Branch("sjt1Pt",jets_.sjt1Pt,"sjt1Pt[nref]/F");
    t->Branch("sjt1Eta",jets_.sjt1Eta,"sjt1Eta[nref]/F");
    t->Branch("sjt1Phi",jets_.sjt1Phi,"sjt1Phi[nref]/F");
    t->Branch("sjt1HadFlav",jets_.sjt1HadFlav,"sjt1HadFlav[nref]/F");
    t->Branch("sjt1ParFlav",jets_.sjt1ParFlav,"sjt1HadFlav[nref]/F");
    t->Branch("sjt2E",jets_.sjt2E,"sjt2E[nref]/F");
    t->Branch("sjt2Pt",jets_.sjt2Pt,"sjt2Pt[nref]/F");
    t->Branch("sjt2Eta",jets_.sjt2Eta,"sjt2Eta[nref]/F");
    t->Branch("sjt2Phi",jets_.sjt2Phi,"sjt2Phi[nref]/F");
    t->Branch("sjt2HadFlav",jets_.sjt2HadFlav,"sjt2HadFlav[nref]/F");
    t->Branch("sjt2ParFlav",jets_.sjt2ParFlav,"sjt2HadFlav[nref]/F");
    /*
    t->Branch("sjt1DiscDeepCSVB",jets_.sjt1DiscDeepCSVB,"sjt1DiscDeepCSVB[nref]/F");
    t->Branch("sjt1DiscDeepCSVBB",jets_.sjt1DiscDeepCSVBB,"sjt1DiscDeepCSVBB[nref]/F");
    t->Branch("sjt1DiscDeepCSVC",jets_.sjt1DiscDeepCSVC,"sjt1DiscDeepCSVC[nref]/F");
    t->Branch("sjt2DiscDeepCSVB",jets_.sjt2DiscDeepCSVB,"sjt2DiscDeepCSVB[nref]/F");
    t->Branch("sjt2DiscDeepCSVBB",jets_.sjt2DiscDeepCSVBB,"sjt2DiscDeepCSVBB[nref]/F");
    t->Branch("sjt2DiscDeepCSVC",jets_.sjt2DiscDeepCSVC,"sjt2DiscDeepCSVC[nref]/F");
    */
    t->Branch("rsjt1E",jets_.rsjt1E,"rsjt1E[nref]/F");
    t->Branch("rsjt1Pt",jets_.rsjt1Pt,"rsjt1Pt[nref]/F");
    t->Branch("rsjt1Eta",jets_.rsjt1Eta,"rsjt1Eta[nref]/F");
    t->Branch("rsjt1Phi",jets_.rsjt1Phi,"rsjt1Phi[nref]/F");
    t->Branch("rsjt2E",jets_.rsjt2E,"rsjt2E[nref]/F");
    t->Branch("rsjt2Pt",jets_.rsjt2Pt,"rsjt2Pt[nref]/F");
    t->Branch("rsjt2Eta",jets_.rsjt2Eta,"rsjt2Eta[nref]/F");
    t->Branch("rsjt2Phi",jets_.rsjt2Phi,"rsjt2Phi[nref]/F");

  }


  // jet ID information, jet composition
  if(doHiJetID_){

    t->Branch("chargedMax", jets_.chargedMax,"chargedMax[nref]/F");
    t->Branch("chargedSum", jets_.chargedSum,"chargedSum[nref]/F");
    t->Branch("chargedN", jets_.chargedN,"chargedN[nref]/I");
    t->Branch("chargedHardSum", jets_.chargedHardSum,"chargedHardSum[nref]/F");
    t->Branch("chargedHardN", jets_.chargedHardN,"chargedHardN[nref]/I");

    t->Branch("photonMax", jets_.photonMax,"photonMax[nref]/F");
    t->Branch("photonSum", jets_.photonSum,"photonSum[nref]/F");
    t->Branch("photonN", jets_.photonN,"photonN[nref]/I");
    t->Branch("photonHardSum", jets_.photonHardSum,"photonHardSum[nref]/F");
    t->Branch("photonHardN", jets_.photonHardN,"photonHardN[nref]/I");

    t->Branch("neutralMax", jets_.neutralMax,"neutralMax[nref]/F");
    t->Branch("neutralSum", jets_.neutralSum,"neutralSum[nref]/F");
    t->Branch("neutralN", jets_.neutralN,"neutralN[nref]/I");

    t->Branch("hcalSum", jets_.hcalSum,"hcalSum[nref]/F");
    t->Branch("ecalSum", jets_.ecalSum,"ecalSum[nref]/F");

    t->Branch("eMax", jets_.eMax,"eMax[nref]/F");
    t->Branch("eSum", jets_.eSum,"eSum[nref]/F");
    t->Branch("eN", jets_.eN,"eN[nref]/I");

    t->Branch("muMax", jets_.muMax,"muMax[nref]/F");
    t->Branch("muSum", jets_.muSum,"muSum[nref]/F");
    t->Branch("muN", jets_.muN,"muN[nref]/I");
  }

  if(doStandardJetID_){
    t->Branch("fHPD",jets_.fHPD,"fHPD[nref]/F");
    t->Branch("fRBX",jets_.fRBX,"fRBX[nref]/F");
    t->Branch("n90",jets_.n90,"n90[nref]/I");
    t->Branch("fSubDet1",jets_.fSubDet1,"fSubDet1[nref]/F");
    t->Branch("fSubDet2",jets_.fSubDet2,"fSubDet2[nref]/F");
    t->Branch("fSubDet3",jets_.fSubDet3,"fSubDet3[nref]/F");
    t->Branch("fSubDet4",jets_.fSubDet4,"fSubDet4[nref]/F");
    t->Branch("restrictedEMF",jets_.restrictedEMF,"restrictedEMF[nref]/F");
    t->Branch("nHCAL",jets_.nHCAL,"nHCAL[nref]/I");
    t->Branch("nECAL",jets_.nECAL,"nECAL[nref]/I");
    t->Branch("apprHPD",jets_.apprHPD,"apprHPD[nref]/F");
    t->Branch("apprRBX",jets_.apprRBX,"apprRBX[nref]/F");

    //  t->Branch("hitsInN90",jets_.n90,"hitsInN90[nref]");
    t->Branch("n2RPC",jets_.n2RPC,"n2RPC[nref]/I");
    t->Branch("n3RPC",jets_.n3RPC,"n3RPC[nref]/I");
    t->Branch("nRPC",jets_.nRPC,"nRPC[nref]/I");

    t->Branch("fEB",jets_.fEB,"fEB[nref]/F");
    t->Branch("fEE",jets_.fEE,"fEE[nref]/F");
    t->Branch("fHB",jets_.fHB,"fHB[nref]/F");
    t->Branch("fHE",jets_.fHE,"fHE[nref]/F");
    t->Branch("fHO",jets_.fHO,"fHO[nref]/F");
    t->Branch("fLong",jets_.fLong,"fLong[nref]/F");
    t->Branch("fShort",jets_.fShort,"fShort[nref]/F");
    t->Branch("fLS",jets_.fLS,"fLS[nref]/F");
    t->Branch("fHFOOT",jets_.fHFOOT,"fHFOOT[nref]/F");
  }

  // Jet ID
  if(doMatch_){
    if(!skipCorrections_) t->Branch("matchedPt", jets_.matchedPt,"matchedPt[nref]/F");
    t->Branch("matchedR", jets_.matchedR,"matchedR[nref]/F");
  }

  // b-jet discriminators
  if (doLifeTimeTagging_) {
    t->Branch("jtDiscCSVV2",jets_.jtDiscCSVV2,"jtDiscCSVV2[nref]/F");
    t->Branch("jtDiscDeepCSVB",jets_.jtDiscDeepCSVB,"jtDiscDeepCSVB[nref]/F");
    t->Branch("jtDiscDeepCSVBB",jets_.jtDiscDeepCSVBB,"jtDiscDeepCSVBB[nref]/F");
    t->Branch("jtDiscDeepCSVC",jets_.jtDiscDeepCSVC,"jtDiscDeepCSVC[nref]/F");
    t->Branch("jtDiscDeepFlavourB",jets_.jtDiscDeepFlavourB,"jtDiscDeepFlavourB[nref]/F");
    t->Branch("jtDiscDeepFlavourBB",jets_.jtDiscDeepFlavourBB,"jtDiscDeepFlavourBB[nref]/F");
    t->Branch("jtDiscDeepFlavourLEPB",jets_.jtDiscDeepFlavourLEPB,"jtDiscDeepFlavourLEPB[nref]/F");
    t->Branch("jtDiscDeepFlavourC",jets_.jtDiscDeepFlavourC,"jtDiscDeepFlavourC[nref]/F");
    t->Branch("jtDiscProb",jets_.jtDiscProb,"jtDiscProb[nref]/F");


    if (doLifeTimeTaggingExtra_) {
      t->Branch("nsvtx",    jets_.nsvtx,    "nsvtx[nref]/I");
      t->Branch("svtxntrk", &jets_.svtxntrk);
      //t->Branch("svtxdl",   &jets_.svtxdl);
      t->Branch("svtxdls",  &jets_.svtxdls);
      //t->Branch("svtxdl2d", &jets_.svtxdl2d);
      t->Branch("svtxdls2d", &jets_.svtxdls2d);
      t->Branch("svtxm",    &jets_.svtxm);
      t->Branch("svtxpt",   &jets_.svtxpt);
      t->Branch("svtxmcorr", &jets_.svtxmcorr);
      
      t->Branch("svtxTrPt",   &jets_.svtxTrPt);
      t->Branch("svtxTrEta",   &jets_.svtxTrEta);
      t->Branch("svtxTrPhi",   &jets_.svtxTrPhi);
   
      t->Branch("nselIPtrk", jets_.nselIPtrk,"nselIPtrk[nref]/I");
      t->Branch("nIP", &jets_.nIP,"nIP/I");
      t->Branch("ipPt",jets_.ipPt,"ipPt[nIP]/F");
      t->Branch("ipEta",jets_.ipEta,"ipEta[nIP]/F");
      t->Branch("ipPhi",jets_.ipPhi,"ipPhi[nIP]/F");
      t->Branch("ipProb",jets_.ipProb,"ipProb[nIP]/F");
      t->Branch("ip3dSig",jets_.ip3dSig,"ip3dSig[nIP]/F");

      if (doTrackMatching_) {
	      t->Branch("ipPtMatch", jets_.ipPtMatch, "ipPtMatch[nIP]/F");
        t->Branch("ipEtaMatch",jets_.ipEtaMatch,"ipEtaMatch[nIP]/F");
        t->Branch("ipPhiMatch",jets_.ipPhiMatch,"ipPhiMatch[nIP]/F");
	      t->Branch("ipMatchStatus", jets_.ipMatchStatus, "ipMatchStatus[nIP]/I");

        t->Branch("ipInSV", jets_.ipInSV, "ipInSV[nIP]/I");
        t->Branch("ipSvtxdls", jets_.ipSvtxdls, "ipSvtxdls[nIP]/F");
        t->Branch("ipSvtxdls2d", jets_.ipSvtxdls2d, "ipSvtxdls2d[nIP]/F");
        t->Branch("ipSvtxm", jets_.ipSvtxm, "ipSvtxm[nIP]/F");
        t->Branch("ipSvtxmcorr", jets_.ipSvtxmcorr, "ipSvtxmcorr[nIP]/F");
      }
    }
  }

  if(isMC_){

    // Only matched gen jets
    t->Branch("refpt",jets_.refpt,"refpt[nref]/F");
    t->Branch("refeta",jets_.refeta,"refeta[nref]/F");
    t->Branch("refphi",jets_.refphi,"refphi[nref]/F");

    t->Branch("refdphijt",jets_.refdphijt,"refdphijt[nref]/F");
    t->Branch("refdrjt",jets_.refdrjt,"refdrjt[nref]/F");
    // matched parton
    t->Branch("refparton_pt",jets_.refparton_pt,"refparton_pt[nref]/F");
    t->Branch("refparton_flavor",jets_.refparton_flavor,"refparton_flavor[nref]/I");

    if(doSubJets_) {
      t->Branch("refIsHardest",jets_.refIsHardest,"refIsHardest[nref]/O");
      //t->Branch("refptG",jets_.refptG,"refptG[nref]/F");
      //t->Branch("refetaG",jets_.refetaG,"refetaG[nref]/F");
      //t->Branch("refphiG",jets_.refphiG,"refphiG[nref]/F");
    }

    if(fillGenJets_){
      // For all gen jets, matched or unmatched
      t->Branch("ngen",&jets_.ngen,"ngen/I");
      t->Branch("genmatchindex",jets_.genmatchindex,"genmatchindex[ngen]/I");
      t->Branch("genpt",jets_.genpt,"genpt[ngen]/F");
      t->Branch("geneta",jets_.geneta,"geneta[ngen]/F");
      t->Branch("genphi",jets_.genphi,"genphi[ngen]/F");

      t->Branch("npar",&jets_.npar,"npar/I");
      t->Branch("parpt",jets_.parpt,"parpt[npar]/F");
      t->Branch("pareta",jets_.pareta,"pareta[npar]/F");
      t->Branch("parphi",jets_.parphi,"parphi[npar]/F");

      t->Branch("parNb",jets_.parNb,"parNb[npar]/I");
      t->Branch("parNc",jets_.parNc,"parNc[npar]/I");
      t->Branch("parHasGSPB",jets_.parHasGSPB,"parHasGSPB[npar]/O");
      t->Branch("parHasGSPC",jets_.parHasGSPC,"parHasGSPC[npar]/O");
    
      if(doSubJets_) {
	t->Branch("genIsHardest",jets_.genIsHardest,"genIsHardest[ngen]/O");
	/*
        t->Branch("genptG",jets_.genptG,"genptG[ngen]/F");
        t->Branch("genetaG",jets_.genetaG,"genetaG[ngen]/F");
        t->Branch("genphiG",jets_.genphiG,"genphiG[ngen]/F");
	*/
	t->Branch("gsjt1E",jets_.gsjt1E,"gsjt1E[ngen]/F");
	t->Branch("gsjt1Pt",jets_.gsjt1Pt,"gsjt1Pt[ngen]/F");
	t->Branch("gsjt1Eta",jets_.gsjt1Eta,"gsjt1Eta[ngen]/F");
	t->Branch("gsjt1Phi",jets_.gsjt1Phi,"gsjt1Phi[ngen]/F");
	t->Branch("gsjt2Pt",jets_.gsjt2Pt,"gsjt2Pt[ngen]/F");
	t->Branch("gsjt2Eta",jets_.gsjt2Eta,"gsjt2Eta[ngen]/F");
	t->Branch("gsjt2Phi",jets_.gsjt2Phi,"gsjt2Phi[ngen]/F");
	/*
        t->Branch("parptG",jets_.parptG,"parptG[npar]/F");
        t->Branch("paretaG",jets_.paretaG,"paretaG[npar]/F");
        t->Branch("parphiG",jets_.parphiG,"parphiG[npar]/F");
	*/
	t->Branch("psjt1E",jets_.psjt1E,"psjt1E[npar]/F");
	t->Branch("psjt1Pt",jets_.psjt1Pt,"psjt1Pt[npar]/F");
	t->Branch("psjt1Eta",jets_.psjt1Eta,"psjt1Eta[npar]/F");
	t->Branch("psjt1Phi",jets_.psjt1Phi,"psjt1Phi[npar]/F");
	t->Branch("psjt2Pt",jets_.psjt2Pt,"psjt2Pt[npar]/F");
	t->Branch("psjt2Eta",jets_.psjt2Eta,"psjt2Eta[npar]/F");
	t->Branch("psjt2Phi",jets_.psjt2Phi,"psjt2Phi[npar]/F");
      }
    }
  }

}


void
HiInclusiveJetAnalyzer::analyze(const Event& iEvent,
				const EventSetup& iSetup) {

  edm::Handle<pat::JetCollection> patjets;
  iEvent.getByToken(jetTagPat_, patjets);

  edm::Handle<pat::JetCollection> patsubjets;
  iEvent.getByToken(subJetTagPat_, patsubjets);

  edm::Handle<JetView> matchedjets;
  iEvent.getByToken(matchTag_, matchedjets);

  edm::Handle<PFCandidateCollection> pfCandidates;
  iEvent.getByToken(pfCandidateLabel_,pfCandidates);

  edm::Handle<JetFlavourInfoMatchingCollection> jetFlavourInfos;
  //edm::Handle<JetFlavourInfoMatchingCollection> subjetFlavourInfos;

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleSrc_, genParticles);
  
  if(doSubJets_){ 
    iEvent.getByToken(groomedJetsToken_, groomedJets);
    iEvent.getByToken(groomedGenJetsToken_, groomedGenJets);
    iEvent.getByToken(groomedPartonJetsToken_, groomedPartonJets);
  }
  if(doExtendedFlavorTagging_){
    iEvent.getByToken(jetFlavourInfosToken_, jetFlavourInfos );
    //if(useSubjets_)iEvent.getByToken(subjetFlavourInfosToken_, subjetFlavourInfos);     
  }


  jets_.nref = 0;
  jets_.nIP = 0;
  int npatjets = (int) patjets->size();
  for(int j = 0; j < npatjets; j++){

    const pat::Jet& jet = (*patjets)[j];

    if(jet.pt() < jetPtMin_) continue;
    if (useJEC_ ){
      jets_.rawpt[jets_.nref]=(*patjets)[j].correctedJet("Uncorrected").pt();
    }

    math::XYZVector jetDir = jet.momentum().Unit();

    if(doLifeTimeTagging_){
      jets_.jtDiscCSVV2[jets_.nref]=(*patjets)[j].bDiscriminator(combinedSVV2BJetTags_);
      jets_.jtDiscDeepCSVB[jets_.nref]=(*patjets)[j].bDiscriminator(deepCSVBJetTags_);
      jets_.jtDiscDeepCSVBB[jets_.nref]=(*patjets)[j].bDiscriminator(deepCSVBBJetTags_);
      jets_.jtDiscDeepCSVC[jets_.nref]=(*patjets)[j].bDiscriminator(deepCSVCJetTags_);

      jets_.jtDiscDeepFlavourB[jets_.nref]=(*patjets)[j].bDiscriminator(deepFlavourBJetTags_);
      jets_.jtDiscDeepFlavourBB[jets_.nref]=(*patjets)[j].bDiscriminator(deepFlavourBBJetTags_);
      jets_.jtDiscDeepFlavourLEPB[jets_.nref]=(*patjets)[j].bDiscriminator(deepFlavourLEPBJetTags_);
      jets_.jtDiscDeepFlavourC[jets_.nref]=(*patjets)[j].bDiscriminator(deepFlavourCJetTags_);

      jets_.jtDiscProb[jets_.nref]=(*patjets)[j].bDiscriminator(jetPBJetTags_);

      if(doLifeTimeTaggingExtra_){
	if((*patjets)[j].hasTagInfo(svTagInfos_.c_str()) ){
	  
	  const CandSecondaryVertexTagInfo *svTagInfo = (*patjets)[j].tagInfoCandSecondaryVertex(svTagInfos_.c_str());
	  
	  jets_.nsvtx[jets_.nref]     = svTagInfo->nVertices();
	  
	  std::vector<int> svtxntrks;
	  std::vector<float> svjetDr, svtxdl, svtxdls, svtxdl2d, svtxdls2d, svtxm, svtxmcorr, svtxpt, svtxTrPt, svtxTrEta, svtxTrPhi;
	  std::vector<int> svType;
	  if(jets_.nsvtx[jets_.nref]==0){
	    svtxntrks.push_back(-1);
	    svjetDr.push_back(-1);
	    svtxdl.push_back(-999);
	    svtxdls.push_back(-999);
	    svtxdl2d.push_back(-999);
	    svtxm.push_back(-10);
	    svtxmcorr.push_back(-10);
	    svtxpt.push_back(-1);	
	    svType.push_back(-9);
	    svtxTrPt.push_back(-1);
	    svtxTrEta.push_back(-1);
	    svtxTrPhi.push_back(-1);
	  }
	  
	  for (int ivtx=0; ivtx<jets_.nsvtx[jets_.nref]; ivtx++) {
	    svtxntrks.push_back(svTagInfo->nVertexTracks(ivtx));
	    svjetDr.push_back(deltaR(svTagInfo->flightDirection(ivtx),jetDir));
	    
	    // this is the 3d flight distance, for 2-D use (0,true)
	    Measurement1D m1D = svTagInfo->flightDistance(ivtx);
	    svtxdl.push_back(m1D.value());
	    svtxdls.push_back(m1D.significance());
	    Measurement1D m2D = svTagInfo->flightDistance(ivtx,true);
	    svtxdl2d.push_back(m2D.value());
	    svtxdls2d.push_back(m2D.significance());
	    //const Vertex& svtx = svTagInfo->secondaryVertex(ivtx);
	    //'const Vertex&' from expression of type 'const VertexCompositePtrCandidate'
	    const VertexCompositePtrCandidate svtx = svTagInfo->secondaryVertex(ivtx);
	    double svtxM = svtx.p4().mass();
	    double svtxPt = svtx.p4().pt();
	    //if(jets_.nsvtx[jets_.nref]>1) cout << "svtxm: "<< svtxM << " svtxpt: "<< svtxPt << endl;
	    svtxm.push_back(svtxM); 
	    svtxpt.push_back(svtxPt);
	    //try out the corrected mass (http://arxiv.org/pdf/1504.07670v1.pdf)
	    //mCorr=srqt(m^2+p^2sin^2(th)) + p*sin(th)
	    double sinth = svtx.p4().Vect().Unit().Cross(svTagInfo->flightDirection(0).unit()).Mag2();
	    sinth = sqrt(sinth);
	    svtxmcorr.push_back(sqrt(pow(svtxM,2)+(pow(svtxPt,2)*pow(sinth,2)))+svtxPt*sinth);
	    
	    std::vector<edm::Ptr<reco::Candidate> > svtxTracks = svTagInfo->vertexTracks(ivtx);
	    for(unsigned int itrk=0; itrk<svtxTracks.size(); itrk++){
	      svtxTrPt.push_back(svtxTracks.at(itrk)->pt());
	      svtxTrEta.push_back(svtxTracks.at(itrk)->eta());
	      svtxTrPhi.push_back(svtxTracks.at(itrk)->phi());
	    }
	    
	    
	  } // SV loop
	  jets_.svtxntrk.push_back(svtxntrks);
	  jets_.svJetDeltaR.push_back(svjetDr);
	  jets_.svtxdl.push_back(svtxdl);
	  jets_.svtxdls.push_back(svtxdls);
	  jets_.svtxdl2d.push_back(svtxdl2d);
	  jets_.svtxdls2d.push_back(svtxdls2d);
	  jets_.svtxm.push_back(svtxm);
	  jets_.svtxmcorr.push_back(svtxmcorr);
	  jets_.svtxpt.push_back(svtxpt);
	  jets_.svtxTrPt.push_back(svtxTrPt);
	  jets_.svtxTrEta.push_back(svtxTrEta);
	  jets_.svtxTrPhi.push_back(svtxTrPhi);
	}

	int counts = 0;
	for (auto genPart : *genParticles) { 
	  double partPt = genPart.pt();
	  if (partPt >= 1.) counts++;
	}
	//std::cout << "nb of gen particles total, ntuplizer: " << (*genParticles).size() << std::endl;
	//std::cout << "nb of gen particles with pt >= 1 GeV, ntuplizer: " << counts << std::endl;
	
	if((*patjets)[j].hasTagInfo(ipTagInfos_.c_str()) ){
	  const CandIPTagInfo *ipTagInfo = (*patjets)[j].tagInfoCandIP(ipTagInfos_.c_str());
	  int it = 0;
	  for (const auto & trk : ipTagInfo->selectedTracks()){
	    jets_.ipPt[jets_.nIP + it] = trk->pt();
	    jets_.ipEta[jets_.nIP + it] = trk->eta();
	    jets_.ipPhi[jets_.nIP + it] = trk->phi();
	    jets_.ipProb[jets_.nIP + it] = ipTagInfo->probabilities(0)[it];
	    reco::btag::TrackIPData data = ipTagInfo->impactParameterData()[it];
	    jets_.ip3dSig[jets_.nIP + it] = data.ip3d.significance();
	    if (doTrackMatching_) {
	      int partID = trkGenPartMatch(trk, *genParticles, 1.);
	      jets_.ipPtMatch[jets_.nIP + it] = (*genParticles)[partID].pt();
	      jets_.ipEtaMatch[jets_.nIP + it] = (*genParticles)[partID].eta();
	      jets_.ipPhiMatch[jets_.nIP + it] = (*genParticles)[partID].phi();
	      jets_.ipMatchStatus[jets_.nIP + it] = (*genParticles)[partID].status();

        jets_.ipInSV[jets_.nIP + it] = 0;
        jets_.ipSvtxdls[jets_.nIP + it] = - 1000000.;
        jets_.ipSvtxdls2d[jets_.nIP + it] = - 1000000.;
        jets_.ipSvtxm[jets_.nIP + it] = - 1000000.;
        jets_.ipSvtxmcorr[jets_.nIP + it] = - 1000000.;

        // Add if in SV
        if((*patjets)[j].hasTagInfo(svTagInfos_.c_str()) ){
	        const CandSecondaryVertexTagInfo *svTagInfo = (*patjets)[j].tagInfoCandSecondaryVertex(svTagInfos_.c_str());
          //std::vector<reco::CandidatePtr> svTracks = svTagInfo->vertexTracks();
          //int whichTrackInSV = trkInVector(trk, svTracks);

          // Need to find in which SV the track is and save the properties of the SV
          for (int ivtx = 0; ivtx < jets_.nsvtx[jets_.nref]; ivtx++) { 
            // look for track in SV
            std::vector<reco::CandidatePtr> isvTracks = svTagInfo->vertexTracks(ivtx);
            int whichTrackInSV = trkInVector(trk, isvTracks);
            if (whichTrackInSV < 0) continue;
            jets_.ipInSV[jets_.nIP + it] = 1;

            // if track is found in this SV, save the SV properties 
            // this is the 3d flight distance, for 2-D use (ivtx, true)
            Measurement1D m1D = svTagInfo->flightDistance(ivtx);
            jets_.ipSvtxdls[jets_.nIP + it] = m1D.significance();
            Measurement1D m2D = svTagInfo->flightDistance(ivtx, true);
            jets_.ipSvtxdls2d[jets_.nIP + it] = m2D.significance();
            const VertexCompositePtrCandidate svtx = svTagInfo->secondaryVertex(ivtx);
            double svtxM = svtx.p4().mass();
            double svtxPt = svtx.p4().pt();
            jets_.ipSvtxm[jets_.nIP + it] = svtxM; 
            //try out the corrected mass (http://arxiv.org/pdf/1504.07670v1.pdf)
            //mCorr=srqt(m^2+p^2sin^2(th)) + p*sin(th)
            double sinth = svtx.p4().Vect().Unit().Cross(svTagInfo->flightDirection(0).unit()).Mag2();
            sinth = sqrt(sinth);
            double svtxMcorr = sqrt(pow(svtxM,2)+(pow(svtxPt,2)*pow(sinth,2)))+svtxPt*sinth;
            jets_.ipSvtxmcorr[jets_.nIP + it] = svtxMcorr;
            break;
          } // end SV loop
        } // end if has SV tag info
	    } // end if track matching
	    it++;
	  } // end track loop
	  jets_.nIP +=  it;
	  jets_.nselIPtrk[jets_.nref] = it;
	}
      }
    }
    if(doHiJetID_){
      // Jet ID variables

      jets_.muMax[jets_.nref] = 0;
      jets_.muSum[jets_.nref] = 0;
      jets_.muN[jets_.nref] = 0;

      jets_.eMax[jets_.nref] = 0;
      jets_.eSum[jets_.nref] = 0;
      jets_.eN[jets_.nref] = 0;

      jets_.neutralMax[jets_.nref] = 0;
      jets_.neutralSum[jets_.nref] = 0;
      jets_.neutralN[jets_.nref] = 0;

      jets_.photonMax[jets_.nref] = 0;
      jets_.photonSum[jets_.nref] = 0;
      jets_.photonN[jets_.nref] = 0;
      jets_.photonHardSum[jets_.nref] = 0;
      jets_.photonHardN[jets_.nref] = 0;

      jets_.chargedMax[jets_.nref] = 0;
      jets_.chargedSum[jets_.nref] = 0;
      jets_.chargedN[jets_.nref] = 0;
      jets_.chargedHardSum[jets_.nref] = 0;
      jets_.chargedHardN[jets_.nref] = 0;

      jets_.hcalSum[jets_.nref] = 0;
      jets_.ecalSum[jets_.nref] = 0;


      for(unsigned int icand = 0; icand < pfCandidates->size(); ++icand){
        const PFCandidate& track = (*pfCandidates)[icand];
        double dr = deltaR(jet,track);
        if(dr < rParam){
	  double ptcand = track.pt();
	  int pfid = track.particleId();

	  switch(pfid){
	  case 1:
	    jets_.chargedSum[jets_.nref] += ptcand;
	    jets_.chargedN[jets_.nref] += 1;
	    if(ptcand > hardPtMin_){
	      jets_.chargedHardSum[jets_.nref] += ptcand;
	      jets_.chargedHardN[jets_.nref] += 1;
	    }
	    if(ptcand > jets_.chargedMax[jets_.nref]) jets_.chargedMax[jets_.nref] = ptcand;
	    break;

	  case 2:
	    jets_.eSum[jets_.nref] += ptcand;
	    jets_.eN[jets_.nref] += 1;
	    if(ptcand > jets_.eMax[jets_.nref]) jets_.eMax[jets_.nref] = ptcand;
	    break;

	  case 3:
	    jets_.muSum[jets_.nref] += ptcand;
	    jets_.muN[jets_.nref] += 1;
	    if(ptcand > jets_.muMax[jets_.nref]) jets_.muMax[jets_.nref] = ptcand;
	    break;

	  case 4:
	    jets_.photonSum[jets_.nref] += ptcand;
	    jets_.photonN[jets_.nref] += 1;
	    if(ptcand > hardPtMin_){
	      jets_.photonHardSum[jets_.nref] += ptcand;
	      jets_.photonHardN[jets_.nref] += 1;
	    }
	    if(ptcand > jets_.photonMax[jets_.nref]) jets_.photonMax[jets_.nref] = ptcand;
	    break;

	  case 5:
	    jets_.neutralSum[jets_.nref] += ptcand;
	    jets_.neutralN[jets_.nref] += 1;
	    if(ptcand > jets_.neutralMax[jets_.nref]) jets_.neutralMax[jets_.nref] = ptcand;
	    break;

	  default:
	    break;

	  }
	}
      }

    }

    if(doMatch_){

      // Alternative reconstruction matching (PF for calo, calo for PF)
      double drMin = 100;
      for(unsigned int imatch = 0 ; imatch < matchedjets->size(); ++imatch){
	const reco::Jet& mjet = (*matchedjets)[imatch];

	double dr = deltaR(jet,mjet);
	if(dr < drMin){
	  jets_.matchedPt[jets_.nref] = mjet.pt();
	  jets_.matchedR[jets_.nref] = dr;
	  drMin = dr;
	}
      }

    }

    jets_.jtpt[jets_.nref] = jet.pt();
    jets_.jteta[jets_.nref] = jet.eta();
    jets_.jtphi[jets_.nref] = jet.phi();

   if(doSubJets_){
     int igj =  getGroomedJetIndex(jet);
     //if(igj != (int)j ) cout<<" mismatched (reco) index !"<<endl;     
     if(igj > -1){

       const reco::Jet& gjet = (*groomedJets)[igj];
       //jets_.jtptG[jets_.nref] = gjet.pt();
       //jets_.jtetaG[jets_.nref] = gjet.eta();
       //jets_.jtphiG[jets_.nref] = gjet.phi();
       if(gjet.jetArea()>0.5) jets_.jtIsHardest[jets_.nref] = true;
       else jets_.jtIsHardest[jets_.nref] = false;

       if(gjet.numberOfDaughters()>0) {
	 const Candidate & sjt1 = *gjet.daughter(0);       
	 jets_.sjt1E[jets_.nref] = sjt1.energy();
	 jets_.sjt1Pt[jets_.nref] = sjt1.pt();
	 jets_.sjt1Eta[jets_.nref] = sjt1.eta();
	 jets_.sjt1Phi[jets_.nref] = sjt1.phi();

	 for(unsigned int j1 = 0; j1 < patsubjets->size(); j1++){
	   const pat::Jet& subjet = (*patsubjets)[j1];
	   if(acos(cos(subjet.phi() - sjt1.phi()))<0.001 && fabs(subjet.eta() - sjt1.eta())<0.001)
	     {
	       jets_.sjt1HadFlav[jets_.nref] = subjet.hadronFlavour();
	       jets_.sjt1ParFlav[jets_.nref] = subjet.partonFlavour();	       
	       //jets_.sjt1DiscDeepCSVB[jets_.nref]=subjet.bDiscriminator(deepCSVBSubJetTags_);
	       //jets_.sjt1DiscDeepCSVBB[jets_.nref]=subjet.bDiscriminator(deepCSVBBSubJetTags_);
	       //jets_.sjt1DiscDeepCSVC[jets_.nref]=subjet.bDiscriminator(deepCSVCSubJetTags_);
	     }
	 }  	     
	 	 
	 if(gjet.numberOfDaughters()>1) {
	   const Candidate & sjt2 = *gjet.daughter(1);
	   jets_.sjt2E[jets_.nref] = sjt2.energy();
	   jets_.sjt2Pt[jets_.nref] = sjt2.pt();
	   jets_.sjt2Eta[jets_.nref] = sjt2.eta();
	   jets_.sjt2Phi[jets_.nref] = sjt2.phi();

	   for(unsigned int j2 = 0; j2 < patsubjets->size(); j2++){
	     const pat::Jet& subjet = (*patsubjets)[j2];
	     if(acos(cos(subjet.phi() - sjt2.phi()))<0.001 && fabs(subjet.eta() - sjt2.eta())<0.001)
	       {
		 jets_.sjt2HadFlav[jets_.nref] = subjet.hadronFlavour();
		 jets_.sjt2ParFlav[jets_.nref] = subjet.partonFlavour();	       
		 //jets_.sjt2DiscDeepCSVB[jets_.nref]=subjet.bDiscriminator(deepCSVBSubJetTags_);
		 //jets_.sjt2DiscDeepCSVBB[jets_.nref]=subjet.bDiscriminator(deepCSVBBSubJetTags_);
		 //jets_.sjt2DiscDeepCSVC[jets_.nref]=subjet.bDiscriminator(deepCSVCSubJetTags_);
	       }
	   }  	     
	   
	 }
	 else{
	   jets_.sjt2E[jets_.nref] = -1;
	   jets_.sjt2Pt[jets_.nref] = -1;
	   jets_.sjt2Eta[jets_.nref] = -999;
	   jets_.sjt2Phi[jets_.nref] = -999;
	 }
       }
       else{
	 jets_.sjt1E[jets_.nref] = -1;
	 jets_.sjt1Pt[jets_.nref] = -1;
	 jets_.sjt1Eta[jets_.nref] = -999;
	 jets_.sjt1Phi[jets_.nref] = -999;
       }     
     }
     else{
       jets_.sjt2E[jets_.nref] = -1;
       jets_.sjt2Pt[jets_.nref] = -1;
       jets_.sjt2Eta[jets_.nref] = -999;
       jets_.sjt2Phi[jets_.nref] = -999;
       jets_.sjt1E[jets_.nref] = -1;
       jets_.sjt1Pt[jets_.nref] = -1;
       jets_.sjt1Eta[jets_.nref] = -999;
       jets_.sjt1Phi[jets_.nref] = -999;       
     }
   }
   
   if( (*patjets)[j].isPFJet() && !doSubJets_){
     jets_.jtPfCHF[jets_.nref] = (*patjets)[j].chargedHadronEnergyFraction();
     jets_.jtPfNHF[jets_.nref] = (*patjets)[j].neutralHadronEnergyFraction();
     jets_.jtPfCEF[jets_.nref] = (*patjets)[j].chargedEmEnergyFraction();
     jets_.jtPfNEF[jets_.nref] = (*patjets)[j].neutralEmEnergyFraction();
     jets_.jtPfMUF[jets_.nref] = (*patjets)[j].muonEnergyFraction();
     
     jets_.jtPfCHM[jets_.nref] = (*patjets)[j].chargedHadronMultiplicity();
     jets_.jtPfNHM[jets_.nref] = (*patjets)[j].neutralHadronMultiplicity();
     jets_.jtPfCEM[jets_.nref] = (*patjets)[j].electronMultiplicity();
     jets_.jtPfNEM[jets_.nref] = (*patjets)[j].photonMultiplicity();
     jets_.jtPfMUM[jets_.nref] = (*patjets)[j].muonMultiplicity();
   }
   if(isMC_ ){
     const GenJet * genjet = (*patjets)[j].genJet();
     if(genjet){

       jets_.refpt[jets_.nref] = genjet->pt();
       jets_.refeta[jets_.nref] = genjet->eta();
       jets_.refphi[jets_.nref] = genjet->phi();
       jets_.refdphijt[jets_.nref] = reco::deltaPhi(jet.phi(), genjet->phi());
       jets_.refdrjt[jets_.nref] = deltaR(jet.eta(),jet.phi(),genjet->eta(),genjet->phi());
     
       int igj = getGroomedGenJetIndex(*genjet);
       //if(igj != (int)j ) cout<<" mismatched index (ref) !"<<endl;

       if(igj > -1){
	 const reco::Jet& gjet = (*groomedGenJets)[igj];
	 if(gjet.jetArea()>0.5) jets_.refIsHardest[jets_.nref] = true;
	 else jets_.refIsHardest[jets_.nref] = false;
	 //jets_.refptG[jets_.nref]  = gjet.pt();
	 //jets_.refetaG[jets_.nref] = gjet.eta();
	 //jets_.refphiG[jets_.nref] = gjet.phi();

	 if(gjet.numberOfDaughters()>0) {
	   const Candidate & gsjt1 = *gjet.daughter(0);       
	   jets_.rsjt1E[jets_.nref] = gsjt1.energy();
	   jets_.rsjt1Pt[jets_.nref] = gsjt1.pt();
	   jets_.rsjt1Eta[jets_.nref] = gsjt1.eta();
	   jets_.rsjt1Phi[jets_.nref] = gsjt1.phi();
	   
	   if(gjet.numberOfDaughters()>1) {
	     const Candidate & gsjt2 = *gjet.daughter(1);
	     jets_.rsjt2E[jets_.nref] = gsjt2.energy();
	     jets_.rsjt2Pt[jets_.nref] = gsjt2.pt();
	     jets_.rsjt2Eta[jets_.nref] = gsjt2.eta();
	     jets_.rsjt2Phi[jets_.nref] = gsjt2.phi();
	   }
	   else{
	     jets_.rsjt2Pt[jets_.nref] = -1;
	     jets_.rsjt2Eta[jets_.nref] = -999;
	     jets_.rsjt2Phi[jets_.nref] = -999;
	   }
	 }
	 else{
	   jets_.rsjt1E[jets_.nref] = -1;
	   jets_.rsjt1Pt[jets_.nref] = -1;
	   jets_.rsjt1Eta[jets_.nref] = -999;
	   jets_.rsjt1Phi[jets_.nref] = -999;
	 }     
       }else{
	 jets_.refpt[jets_.nref] = -999.;
	 jets_.refeta[jets_.nref] = -999.;
	 jets_.refphi[jets_.nref] = -999.;
	 jets_.refm[jets_.nref] = -999.;
	 jets_.refdphijt[jets_.nref] = -999.;
	 jets_.refdrjt[jets_.nref] = -999.;
	 
	 /*
	   if(doSubJets_) {
	   jets_.refptG[jets_.nref]  = -999.;
	   jets_.refetaG[jets_.nref] = -999.;
	   jets_.refphiG[jets_.nref] = -999.;	 
	   }
	 */
       }
     }
     
     jets_.jtParFlav[jets_.nref] = (*patjets)[j].partonFlavour();
     jets_.jtHadFlav[jets_.nref] = (*patjets)[j].hadronFlavour();

     bool foundMatch = false;
     for (const JetFlavourInfoMatching& jetFlavourInfoMatching : *jetFlavourInfos) {
       if (deltaR(jet.p4(), jetFlavourInfoMatching.first->p4()) < 1e-6) {
	 JetFlavourInfo jetInfo = jetFlavourInfoMatching.second;
	 const GenParticleRefVector &bHadronsInJet = jetInfo.getbHadrons();
	 const GenParticleRefVector &cHadronsInJet = jetInfo.getcHadrons();

	 jets_.jtNbHad[jets_.nref] = bHadronsInJet.size();
	 jets_.jtNcHad[jets_.nref] = cHadronsInJet.size();

	 const GenParticleRefVector &partonsInJet = jetInfo.getPartons();

	 int nb=0;
	 int nc=0;
	 bool hasBfromGSP = false;
	 bool hasCfromGSP = false;

	 for (GenParticleRefVector::const_iterator it = partonsInJet.begin(); it != partonsInJet.end(); ++it) {
	   int parFlav = (*it)->pdgId();
	   const Candidate* c = (*it).get();

	   if(abs(parFlav)==5){
	     nb++;
	     if(isFromGSP(c)) hasBfromGSP = true;
	   }
	   else if(abs(parFlav)==4){
	     nc++;
	     if(isFromGSP(c)) hasCfromGSP = true;
	   }
	 }

	 jets_.jtNbPar[jets_.nref] = nb;
	 jets_.jtNcPar[jets_.nref] = nc;
	 jets_.jtHasGSPB[jets_.nref] = hasBfromGSP;
	 jets_.jtHasGSPC[jets_.nref] = hasCfromGSP;

	 foundMatch = true;
	 break;
       }
     }
   
     if(!foundMatch) cout<<" didn't find match! "<<endl;
      // matched partons
      const GenParticle & parton = *(*patjets)[j].genParton();

      if((*patjets)[j].genParton()){
	jets_.refparton_pt[jets_.nref] = parton.pt();
	jets_.refparton_flavor[jets_.nref] = parton.pdgId();
      } else {
	jets_.refparton_pt[jets_.nref] = -999;
	jets_.refparton_flavor[jets_.nref] = -999;
      }
    }
    jets_.nref++;
  }

  if(isMC_){

    edm::Handle<vector<GenJet> >genjets;
    iEvent.getByToken(genjetTag_, genjets);
    
    edm::Handle<vector<GenJet> >partonJets;
    iEvent.getByToken(partonjetTag_, partonJets);
  
    jets_.ngen = 0;
    for(unsigned int igen = 0 ; igen < genjets->size(); ++igen){
      const GenJet & genjet = (*genjets)[igen];
      float genjet_pt = genjet.pt();
      jets_.genmatchindex[jets_.ngen] = -1;
      
      if(genjet_pt<genPtMin_) continue;
	
      for(int ijet = 0 ; ijet < jets_.nref; ++ijet){
	double deltaPt = fabs(genjet.pt()-jets_.refpt[ijet]); 
	if(deltaPt > 0.01) continue;
	double deltaEta = fabs(genjet.eta()-jets_.refeta[ijet]);
	if(deltaEta > 0.0001) continue;
	double deltaPhi = fabs(reco::deltaPhi(genjet.phi(), jets_.refphi[ijet])); 
	if(deltaPhi > 0.0001) continue;
	jets_.genmatchindex[jets_.ngen] = (int)ijet;
	break;
      }
      
      jets_.genpt [jets_.ngen] = genjet_pt;
      jets_.geneta[jets_.ngen] = genjet.eta();
      jets_.genphi[jets_.ngen] = genjet.phi();
      
      if(doSubJets_){
	int igj = getGroomedGenJetIndex(genjet);
	if(igj != (int)igen ) cout<<" mismatched index (gen) !"<<endl;
	if(igj > -1){
	  const reco::Jet& gjet = (*groomedGenJets)[igj];
	  if(gjet.jetArea()>0.5) jets_.genIsHardest[jets_.ngen] = true;
	  else jets_.genIsHardest[jets_.ngen] = false;
	  //jets_.genptG[jets_.ngen]  = gjet.pt();
	  //jets_.genetaG[jets_.ngen] = gjet.eta();
	  //jets_.genphiG[jets_.ngen] = gjet.phi();

	  if(gjet.numberOfDaughters()>0) {
	    const Candidate & gsjt1 = *gjet.daughter(0);       
	    jets_.gsjt1E[jets_.ngen] = gsjt1.energy();
	    jets_.gsjt1Pt[jets_.ngen] = gsjt1.pt();
	    jets_.gsjt1Eta[jets_.ngen] = gsjt1.eta();
	    jets_.gsjt1Phi[jets_.ngen] = gsjt1.phi();

	    if(gjet.numberOfDaughters()>1) {
	      const Candidate & gsjt2 = *gjet.daughter(1);
	      jets_.gsjt2E[jets_.ngen] = gsjt2.energy();
	      jets_.gsjt2Pt[jets_.ngen] = gsjt2.pt();
	      jets_.gsjt2Eta[jets_.ngen] = gsjt2.eta();
	      jets_.gsjt2Phi[jets_.ngen] = gsjt2.phi();
	    }
	    else{
	      jets_.gsjt2E[jets_.ngen] = -1;
	      jets_.gsjt2Pt[jets_.ngen] = -1;
	      jets_.gsjt2Eta[jets_.ngen] = -999;
	      jets_.gsjt2Phi[jets_.ngen] = -999;
	    }
	  }
	  else{
	    jets_.gsjt1E[jets_.ngen] = -1;
	    jets_.gsjt1Pt[jets_.ngen] = -1;
	    jets_.gsjt1Eta[jets_.ngen] = -999;
	    jets_.gsjt1Phi[jets_.ngen] = -999;
	  }     
	}
	else{
	    jets_.gsjt1E[jets_.ngen] = -1;
	    jets_.gsjt1Pt[jets_.ngen] = -1;
	    jets_.gsjt1Eta[jets_.ngen] = -999;
	    jets_.gsjt1Phi[jets_.ngen] = -999;
	    jets_.gsjt2E[jets_.ngen] = -1;
	    jets_.gsjt2Pt[jets_.ngen] = -1;
	    jets_.gsjt2Eta[jets_.ngen] = -999;
	    jets_.gsjt2Phi[jets_.ngen] = -999;	    
	}
      }	    
      jets_.ngen++;
    }

    jets_.npar = 0;
    for(unsigned int ipar = 0 ; ipar < partonJets->size(); ++ipar){
      const reco::GenJet & pjet = (*partonJets)[ipar];
      
      if(pjet.pt() < genPtMin_) continue;
      jets_.parpt[jets_.npar]  = pjet.pt();
      jets_.pareta[jets_.npar] = pjet.eta();
      jets_.parphi[jets_.npar] = pjet.phi();
    
      int nb=0;
      int nc=0;
      bool hasBfromGSP = false;
      bool hasCfromGSP = false;

      for (auto constit : pjet.getJetConstituents()) {
	if(abs(constit->pdgId()) ==5){
	  if(isFromGSP(constit.get()) ) hasBfromGSP = true;
	  nb++;
	}
	if(abs(constit->pdgId()) ==4){
	  if(isFromGSP(constit.get()) ) hasCfromGSP = true;
	  nc++;	
	}
      }
      
      jets_.parNb[jets_.npar] = nb;
      jets_.parNc[jets_.npar] = nc;
      jets_.parHasGSPB[jets_.npar] = hasBfromGSP;
      jets_.parHasGSPC[jets_.npar] = hasCfromGSP;
      
      int ipj = getGroomedPartonJetIndex(pjet);
      if(ipj != (int)ipar ) cout<<" mismatched index (gen) !"<<endl;
      if(ipj > -1){
	const reco::Jet& gpjet = (*groomedPartonJets)[ipj];
	//jets_.parptG[jets_.npar]  = gpjet.pt();
	//jets_.paretaG[jets_.npar] = gpjet.eta();
	//jets_.parphiG[jets_.npar] = gpjet.phi();
      		
	if(gpjet.numberOfDaughters()>0) {
	  const Candidate & psjt1 = *gpjet.daughter(0);       
	  jets_.psjt1E[jets_.npar] = psjt1.energy();
	  jets_.psjt1Pt[jets_.npar] = psjt1.pt();
	  jets_.psjt1Eta[jets_.npar] = psjt1.eta();
	  jets_.psjt1Phi[jets_.npar] = psjt1.phi();

	  if(gpjet.numberOfDaughters()>1) {
	    const Candidate & psjt2 = *gpjet.daughter(1);
	    jets_.psjt2E[jets_.npar] = psjt2.energy();
	    jets_.psjt2Pt[jets_.npar] = psjt2.pt();
	    jets_.psjt2Eta[jets_.npar] = psjt2.eta();
	    jets_.psjt2Phi[jets_.npar] = psjt2.phi();
	  }
	  else{
	    jets_.psjt2E[jets_.npar] = -1;
	    jets_.psjt2Pt[jets_.npar] = -1;
	    jets_.psjt2Eta[jets_.npar] = -999;
	    jets_.psjt2Phi[jets_.npar] = -999;
	  }
	}
	else{
	  jets_.psjt1E[jets_.npar] = -1;
	  jets_.psjt1Pt[jets_.npar] = -1;
	  jets_.psjt1Eta[jets_.npar] = -999;
	  jets_.psjt1Phi[jets_.npar] = -999;
	}           
      }
      else{
	jets_.psjt2E[jets_.npar] = -1;
	jets_.psjt2Pt[jets_.npar] = -1;
	jets_.psjt2Eta[jets_.npar] = -999;
	jets_.psjt2Phi[jets_.npar] = -999;
	jets_.psjt1E[jets_.npar] = -1;
	jets_.psjt1Pt[jets_.npar] = -1;
	jets_.psjt1Eta[jets_.npar] = -999;
	jets_.psjt1Phi[jets_.npar] = -999;
      }
      jets_.npar++;
    }
  }
  
  t->Fill();

  memset(&jets_,0,sizeof jets_);
}

//--------------------------------------------------------------------------------------------------
int HiInclusiveJetAnalyzer::getGroomedGenJetIndex(const GenJet jet) const {

  //Find closest soft-dropped gen jet
  double drMin = 0.2;
  int imatch = -1;

  for(unsigned int i = 0 ; i < groomedGenJets->size(); ++i) {
    const reco::Jet& mjet = (*groomedGenJets)[i];

    double dr = deltaR(jet,mjet);
    if(dr < drMin){
      imatch = i;
      drMin = dr;
    }
    if(drMin < 1.0e-6) break;
  }
  //if(drMin>0.2)cout<<" drMin gen : "<<drMin<<endl;  
  return imatch;
}

//--------------------------------------------------------------------------------------------------
int HiInclusiveJetAnalyzer::getGroomedPartonJetIndex(const GenJet jet) const {

  //Find closest soft-dropped gen jet
  double drMin = 0.2;
  int imatch = -1;

  for(unsigned int i = 0 ; i < groomedPartonJets->size(); ++i) {
    const reco::Jet& mjet = (*groomedPartonJets)[i];

    double dr = deltaR(jet,mjet);
    if(dr < drMin){
      imatch = i;
      drMin = dr;
    }
    if(drMin < 1.0e-6) break;
  }
  //if(drMin>0.2)cout<<" drMin parton : "<<drMin<<endl;  
  return imatch;
}

//--------------------------------------------------------------------------------------------------
int HiInclusiveJetAnalyzer::getGroomedJetIndex(const reco::Jet jet) const {
  //Find closest soft-dropped jet
  double drMin = 0.2;
  int imatch = -1;
  for(unsigned int i = 0 ; i < groomedJets->size(); ++i) {
    const reco::Jet& mjet = (*groomedJets)[i];
    
    double dr = deltaR(jet,mjet);
    if(dr < drMin){
      imatch = i;
      drMin = dr;
    }
    if(drMin < 1.0e-6) break;
  }
  //if(drMin>0.2)cout<<" drMin reco : "<<drMin<<endl;  
  return imatch;
}

// from BTagAnalyzer
bool HiInclusiveJetAnalyzer::isFromGSP(const Candidate* c)
{  
  bool isFromGSP = false;
  
  if( c->numberOfMothers() == 1 ) {
    const Candidate* dau = c;
    const Candidate* mom = c->mother();
    while( dau->numberOfMothers() == 1 && !( isHardProcess(mom->status()) && (abs(mom->pdgId())==4 || abs(mom->pdgId())==5) ) ) {
      if( mom->pdgId()==21 )
	{
	  isFromGSP = 1;
	  break;
	}
      dau = mom;
      mom = dau->mother();
    }
  }
  
  return isFromGSP;
}
//---------------------------------------------------------------------
int HiInclusiveJetAnalyzer::trkGenPartMatch(const reco::CandidatePtr track, reco::GenParticleCollection genParticles, float ptCut)
{
  /* 
     Match reconstructed track to gen particle 
  */

  // Closest gen particle ID
  int partID = -1;

  double dRmin = 100.;
  double dR = dRmin;
  double trkEta = track->eta();
  double trkPhi = track->phi();

  for (size_t i = 0; i < genParticles.size(); i++) {
    reco::GenParticle genParticle = genParticles[i];

    double partPt = genParticle.pt();
    if (partPt < ptCut) continue;
	if (genParticle.charge() == 0) continue;

    double partEta = genParticle.eta();
    double partPhi = genParticle.phi();

    //std::cout << "i: " << i << ", pt: " << partPt << ", eta: " << partEta << ", phi: " << partPhi << std::endl;

    double dEta = trkEta - partEta;
    double dPhi = std::acos(std::cos(trkPhi - partPhi));
    dR = std::sqrt((dEta * dEta) + (dPhi * dPhi));

    if (dR < dRmin) {
      dRmin = dR;
      partID = i;
    }
  }
  return partID;
}
//---------------------------------------------------------------------
int HiInclusiveJetAnalyzer::trkInVector(reco::CandidatePtr trk, std::vector<reco::CandidatePtr> tracks) const {
  int indexInVector = -1;
  auto it = find(tracks.begin(), tracks.end(), trk);
  if (it != tracks.end()) {
	indexInVector = it - tracks.begin();
  }
  return indexInVector;
}

//-------------------------------------------------------------------
bool HiInclusiveJetAnalyzer::isHardProcess(const int status)
{
  // assumes Pythia8
  if( status>=21 && status<=29 ) return true;
  return false;
}



DEFINE_FWK_MODULE(HiInclusiveJetAnalyzer);
