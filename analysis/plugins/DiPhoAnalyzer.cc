#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "flashgg/MicroAODFormats/interface/Photon.h"
#include "flashgg/MicroAODFormats/interface/DiPhotonCandidate.h"

#include "TMath.h"
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#define MAX_PU_REWEIGHT 60

using namespace std;
using namespace edm;
using namespace flashgg;

using pat::PackedGenParticle;   

// diphoton tree
struct diphoTree_struc_ {

  int run;
  int event;
  int lumi;
  int nvtx;
  float rho;
  int sampleID;
  float totXsec;
  float pu_weight;
  float pu_n;
  float ptgg;
  float mgg;
  float pt1; 
  float ptOverM1; 
  float eta1; 
  float phi1;
  float r91; 
  float sieie1; 
  float hoe1; 
  float scRawEne1;
  float chiso1; 
  float phoiso1; 
  // float neuiso1;
  float pt2; 
  float ptOverM2; 
  float eta2; 
  float phi2;
  float r92; 
  float sieie2; 
  float hoe2; 
  float scRawEne2;
  float chiso2; 
  float phoiso2; 
  // float neuiso2;
  float vtxX; 
  float vtxY; 
  float vtxZ;
  int genmatch1;   
  int genmatch2;
  float genVtxX; 
  float genVtxY; 
  float genVtxZ;
};

class DiPhoAnalyzer : public edm::EDAnalyzer {
  
public:
  
  explicit DiPhoAnalyzer(const edm::ParameterSet&);
  ~DiPhoAnalyzer();
  
private:
  
  edm::Service<TFileService> fs_;
  
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  void initTreeStructure();

  void SetPuWeights(std::string puWeightFile);
  float GetPUWeight(float pun);

  // collections
  EDGetTokenT<View<reco::Vertex> > vertexToken_;
  edm::InputTag packedGenGamma_;  
  EDGetTokenT<edm::View<flashgg::DiPhotonCandidate> > diPhotonToken_; 
  EDGetTokenT<edm::View<PileupSummaryInfo> > PileUpToken_; 
  edm::InputTag rhoFixedGrid_;
  
  // sample-dependent parameters needed for the analysis
  int dopureweight_;
  int sampleIndex_;
  string puWFileName_;
  float xsec_;    // pb
  float kfac_;

  // to compute weights for pileup
  std::vector<Double_t> puweights_;

  // output tree with several diphoton infos
  TTree *DiPhotonTree;
  diphoTree_struc_ treeDipho_;

  // to keep track of the original number of events
  TH1F *h_entries;
};
   

DiPhoAnalyzer::DiPhoAnalyzer(const edm::ParameterSet& iConfig):

  // collections
  vertexToken_(consumes<View<reco::Vertex> >(iConfig.getUntrackedParameter<InputTag> ("VertexTag", InputTag("offlineSlimmedPrimaryVertices")))),
  packedGenGamma_(iConfig.getUntrackedParameter<InputTag>("packedGenParticles", InputTag("flashggGenPhotons"))),
  diPhotonToken_(consumes<View<flashgg::DiPhotonCandidate> >(iConfig.getUntrackedParameter<InputTag> ("DiPhotonTag", InputTag("flashggDiPhotons")))),
  PileUpToken_(consumes<View<PileupSummaryInfo> >(iConfig.getUntrackedParameter<InputTag> ("PileUpTag", InputTag("addPileupInfo"))))
{ 
  dopureweight_ = iConfig.getUntrackedParameter<int>("dopureweight", 0);
  sampleIndex_  = iConfig.getUntrackedParameter<int>("sampleIndex",0);
  puWFileName_  = iConfig.getParameter<std::string>("puWFileName");   
  xsec_         = iConfig.getUntrackedParameter<double>("xsec",1.); 
  kfac_         = iConfig.getUntrackedParameter<double>("kfac",1.); 
};

DiPhoAnalyzer::~DiPhoAnalyzer() { };

void DiPhoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  // access edm objects                                                                                    
  Handle<View<reco::Vertex> > primaryVertices;
  iEvent.getByToken(vertexToken_,primaryVertices);
  const PtrVector<reco::Vertex>& vtxs = primaryVertices->ptrVector();

  Handle<vector<PackedGenParticle> > genGammas;     
  iEvent.getByLabel(packedGenGamma_, genGammas);       

  Handle<View<flashgg::DiPhotonCandidate> > diPhotons;
  iEvent.getByToken(diPhotonToken_,diPhotons);
  const PtrVector<flashgg::DiPhotonCandidate>& dipho = diPhotons->ptrVector();

  Handle<View< PileupSummaryInfo> > PileupInfos;
  iEvent.getByToken(PileUpToken_,PileupInfos);
  const PtrVector<PileupSummaryInfo>& PileupInfoPointers = PileupInfos->ptrVector();

  Handle<double> objs_rho;
  iEvent.getByLabel("fixedGridRhoAll",objs_rho);


  // trees re-initialization                                                                          
  initTreeStructure();        

  // To keep track of the total number of events
  h_entries->Fill(5);

  // Sample index
  int sampleID = sampleIndex_;

  // Event info
  int run   = iEvent.eventAuxiliary().run();
  int event = iEvent.eventAuxiliary().event();
  int lumi  = iEvent.eventAuxiliary().luminosityBlock(); 

  // Vertices
  int nvtx = vtxs.size(); 

  // Energy density
  float rho = *(objs_rho.product());

  // PU weight for MC only and if requested
  float pu_weight = 1.;
  float pu_n      = -1.;
  if (genGammas->size()>0) {   
    pu_n = 0.;
    for( unsigned int PVI = 0; PVI < PileupInfoPointers.size(); ++PVI) {
      Int_t pu_bunchcrossing = PileupInfoPointers[PVI]->getBunchCrossing();
      if (pu_bunchcrossing ==0) {
	pu_n = PileupInfoPointers[PVI]->getPU_NumInteractions();
      }
    }
    if (dopureweight_) 
      pu_weight = GetPUWeight(pu_n);         
  }

  // x-sec * kFact for MC only 
  float totXsec = -1.;
  if (genGammas->size()>0) totXsec = xsec_ * kfac_;

  // Diphoton candidates
  // chiara: per il momento quello a pT del sistema piu' alto
  // chiara: questi sono tutti i candidati, che se capisco non sono selezionati (da controllare). 
  // La selezione va applicata prima della scelta
  // La scelta e' chiaramente bacata per un po' di eventi, va capito perche'
  float maxDiphoPt = -999.;
  int candIndex = 9999;    // This int will store the index of the best diphoton candidate
  for (unsigned int diphotonlooper =0; diphotonlooper < dipho.size() ; diphotonlooper++){
    float thisPt = dipho[diphotonlooper]->pt();
    if (thisPt>maxDiphoPt) {
      maxDiphoPt = thisPt;
      candIndex  = diphotonlooper;
    }
  }

  // chiara: qui applicare tagli analisi

  // to be kept in the tree
  float ptgg, mgg;
  float pt1, ptOverM1, eta1, phi1;
  float r91, sieie1, hoe1, scRawEne1;
  float chiso1, phoiso1; //, neuiso1;
  float pt2, ptOverM2, eta2, phi2;
  float r92, sieie2, hoe2, scRawEne2;
  float chiso2, phoiso2;  //, neuiso2;
  float vtxX, vtxY, vtxZ;
  int genmatch1, genmatch2;
  float genVtxX, genVtxY, genVtxZ;   // chiara: vuoti

  // run the analysis only if the resonance candidate is found
  if(candIndex<9999) {

    //-------> diphoton system properties 
    ptgg = dipho[candIndex]->pt();
    mgg  = dipho[candIndex]->mass();

    //-------> individual photon properties, cominciamo con questi
    pt1       = dipho[candIndex]->leadingPhoton()->et();
    ptOverM1  = dipho[candIndex]->leadingPhoton()->pt()/mgg;
    eta1      = dipho[candIndex]->leadingPhoton()->eta();
    phi1      = dipho[candIndex]->leadingPhoton()->phi();
    r91       = dipho[candIndex]->leadingPhoton()->r9();
    sieie1    = dipho[candIndex]->leadingPhoton()->sigmaIetaIeta();
    hoe1      = dipho[candIndex]->leadingPhoton()->hadronicOverEm();
    scRawEne1 = dipho[candIndex]->leadingPhoton()->superCluster()->rawEnergy();
    chiso1    = dipho[candIndex]->leadingPhoton()->getpfChgIso03WrtVtx(dipho[candIndex]->getVertex());
    phoiso1   = dipho[candIndex]->leadingPhoton()->getpfPhoIso03();
    // neuiso1   = dipho[candIndex]->leadingPhoton()->getpfNeutIso03();

    pt2       = dipho[candIndex]->subLeadingPhoton()->et();
    ptOverM2  = dipho[candIndex]->subLeadingPhoton()->pt()/mgg;
    eta2      = dipho[candIndex]->subLeadingPhoton()->eta();
    phi2      = dipho[candIndex]->subLeadingPhoton()->phi();
    r92       = dipho[candIndex]->subLeadingPhoton()->r9();
    sieie2    = dipho[candIndex]->subLeadingPhoton()->sigmaIetaIeta();
    hoe2      = dipho[candIndex]->subLeadingPhoton()->hadronicOverEm();
    scRawEne2 = dipho[candIndex]->subLeadingPhoton()->superCluster()->rawEnergy();
    chiso2    = dipho[candIndex]->subLeadingPhoton()->getpfChgIso03WrtVtx(dipho[candIndex]->getVertex());
    phoiso2   = dipho[candIndex]->subLeadingPhoton()->getpfPhoIso03();
    // neuiso2   = dipho[candIndex]->subLeadingPhoton()->getpfNeutIso03();

    //-------> vtx info
    vtxX= dipho[candIndex]->getVertex()->x();
    vtxY= dipho[candIndex]->getVertex()->y();
    vtxZ= dipho[candIndex]->getVertex()->z();

    //-------> generated vtx info
    // chiara: rigirare per accedere alla mc truth completa (check con Pasquale)
    /*
    genVtxX = 0.;
    genVtxY = 0.;
    genVtxZ = 0.;
    if (genParticles->size()>0) {   
      for( unsigned int genLoop =0 ; genLoop < genParticles->size(); genLoop++){
	if( genParticles[genLoop]->pdgId() == 5100039) { 
	  genVtxX = genParticles[genLoop]->vx();
	  genVtxY = genParticles[genLoop]->vy();
	  genVtxZ = genParticles[genLoop]->vz();
	}
	break;
      }
    }
    */

    //-------> photons, MC truth match
    // chiara: se nel tree fossero ok basterebbe fare come le righe commentate sotto.
    // a me paiono bacate nel tree. Faccio il loop a mano, ma devo togliere il check sulla madre che non e' salvata
    // check con Pasquale
    genmatch1 = 0;
    genmatch2 = 0;
    // genmatch1 = (dipho[candIndex]->leadingPhoton()->genMatchType() == Photon::kPrompt); 
    // genmatch2 = (dipho[candIndex]->subLeadingPhoton()->genMatchType() == Photon::kPrompt); 
    if (genGammas->size()>0) {   
      for(vector<PackedGenParticle>::const_iterator igen=genGammas->begin(); igen!=genGammas->end(); ++igen) {
	if( igen->status() != 1 || igen->pdgId() != 22 ) continue; 
	// if( igen->motherRef()->pdgId() <= 25 ) {
	float deta = dipho[candIndex]->leadingPhoton()->eta() - igen->eta();
	float dphi = deltaPhi(dipho[candIndex]->leadingPhoton()->phi(),igen->phi());
	float dr = sqrt(deta*deta + dphi*dphi);
	float pt_change = (dipho[candIndex]->leadingPhoton()->et() - igen->et())/igen->et();
	if (dr<0.3 && fabs(pt_change) < 0.5) {
	  genmatch1 = 1;
	  break;
	}
	//}
      }
      for(vector<PackedGenParticle>::const_iterator igen=genGammas->begin(); igen!=genGammas->end(); ++igen) {
	if( igen->status() != 1 || igen->pdgId() != 22 ) continue; 
	// if( igen->motherRef()->pdgId() <= 25 ) {
	float deta = dipho[candIndex]->subLeadingPhoton()->eta() - igen->eta();
	float dphi = deltaPhi(dipho[candIndex]->subLeadingPhoton()->phi(),igen->phi());
	float dr = sqrt(deta*deta + dphi*dphi);
	float pt_change = (dipho[candIndex]->subLeadingPhoton()->et() - igen->et())/igen->et();
	if (dr<0.3 && fabs(pt_change) < 0.5) {
	  genmatch2 = 1;
	  break;
	}
	//}
      }

    } // MC truth

  } else {    // if we could not find any diphoton candidate

    ptgg     = -999.;
    mgg      = -999.;
    pt1      = -999.;
    ptOverM1 = -999.;
    eta1     = -999.;
    phi1     = -999.;
    r91      = -999.;
    sieie1   = -999.;
    hoe1     = -999.;
    scRawEne1  = -999.;
    chiso1   = -999.;
    phoiso1  = -999.;
    // neuiso1  = -999.;

    pt2      = -999.;
    ptOverM2 = -999.;
    eta2     = -999.;
    phi2     = -999.;
    r92      = -999.;
    sieie2   = -999.;
    hoe2     = -999.;
    scRawEne2  = -999.;
    chiso2   = -999.;
    phoiso2  = -999.;
    // neuiso2  = -999.;

    vtxX = -999.;
    vtxY = -999.;
    vtxZ = -999.;

    genmatch1 = -999;
    genmatch2 = -999;

    genVtxX = -999.;
    genVtxY = -999.;
    genVtxZ = -999.;
  }

  // Variables for the tree
  treeDipho_.run = run;
  treeDipho_.event = event;
  treeDipho_.lumi = lumi;
  treeDipho_.nvtx = nvtx;
  treeDipho_.rho = rho;
  treeDipho_.sampleID = sampleID;  
  treeDipho_.totXsec = totXsec;  
  treeDipho_.pu_weight = pu_weight;
  treeDipho_.pu_n = pu_n;
  treeDipho_.ptgg = ptgg;
  treeDipho_.mgg = mgg;
  treeDipho_.pt1 = pt1;
  treeDipho_.ptOverM1 = ptOverM1;
  treeDipho_.eta1 = eta1;
  treeDipho_.phi1 = phi1;
  treeDipho_.r91 = r91;
  treeDipho_.sieie1 = sieie1;
  treeDipho_.hoe1 = hoe1; 
  treeDipho_.scRawEne1 = scRawEne1;
  treeDipho_.chiso1 = chiso1; 
  treeDipho_.phoiso1 = phoiso1; 
  // treeDipho_.neuiso1 = neuiso1;
  treeDipho_.pt2 = pt2;
  treeDipho_.ptOverM2 = ptOverM2;
  treeDipho_.eta2 = eta2;
  treeDipho_.phi2 = phi2;
  treeDipho_.r92 = r92;
  treeDipho_.sieie2 = sieie2;
  treeDipho_.hoe2 = hoe2; 
  treeDipho_.scRawEne2 = scRawEne2;
  treeDipho_.chiso2 = chiso2; 
  treeDipho_.phoiso2 = phoiso2; 
  // treeDipho_.neuiso2 = neuiso2;
  treeDipho_.vtxX = vtxX;
  treeDipho_.vtxY = vtxY;
  treeDipho_.vtxZ = vtxZ;
  treeDipho_.genmatch1 = genmatch1; 
  treeDipho_.genmatch2 = genmatch2; 
  treeDipho_.genVtxX = genVtxX;
  treeDipho_.genVtxY = genVtxY;
  treeDipho_.genVtxZ = genVtxZ;

  // Filling the trees
  DiPhotonTree->Fill();
}

void DiPhoAnalyzer::beginJob() {

  // checking parameters:
  cout << "chiara, test => doPUrew = " << dopureweight_ 
       << ", sampleIndex = " << sampleIndex_ 
       << ", puWeights = " << puWFileName_ 
       << ", xSec = " << xsec_ 
       << ", kFactor = " << kfac_ << endl; 
  cout << endl;

  // loading weights for pileup if needed
  if (dopureweight_) 
    SetPuWeights(puWFileName_);
  
  // to keep track of the original number of events
  h_entries = fs_->make<TH1F>("h_entries", "h_entries", 10,  0., 10.);

  // Trees
  DiPhotonTree = fs_->make<TTree>("DiPhotonTree","di-photon tree");

  // with all infos
  TString treeEvent  = "run/I:event/I:lumi/I:nvtx/I:rho/F";
  TString treeSample = "sampleID/I:totXsec/F:pu_weight/F:pu_n/F";
  TString treeDiPho  = "ptgg/F:mgg/F";
  // TString treePho1    = "pt1:F/ptOverM1:F/eta1:F/phi1:F/r91:F/sieie1:F/hoe1:F/scRawEne1:F/chiso1:F/phoiso1:F/neuiso1:F/genmatch1:I";
  // TString treePho2    = "pt2:F/ptOverM2:F/eta2:F/phi2:F/r92:F/sieie2:F/hoe2:F/scRawEne2:F/chiso2:F/phoiso2:F/neuiso2:F/genmatch2:I";
  TString treePho1   = "pt1/F:ptOverM1/F:eta1/F:phi1/F:r91/F:sieie1/F:hoe1/F:scRawEne1/F:chiso1/F:phoiso1/F";
  TString treePho2   = "pt2/F:ptOverM2/F:eta2/F:phi2/F:r92/F:sieie2/F:hoe2/F:scRawEne2/F:chiso2/F:phoiso2/F";
  TString treeMcPho1 = "genmatch1/I";
  TString treeMcPho2 = "genmatch2/I";
  TString treeVtx    = "vtxX/F:vtxY/F:vtxZ/F";
  TString treeGenVtx = "genVtxX/F:genVtxY/F:genVtxZ/F";
  DiPhotonTree->Branch("eventinfo",&(treeDipho_.run),treeEvent);
  DiPhotonTree->Branch("sample",&(treeDipho_.sampleID),treeSample);
  DiPhotonTree->Branch("diphotons",&(treeDipho_.ptgg),treeDiPho);
  DiPhotonTree->Branch("gamma1",&(treeDipho_.pt1),treePho1);
  DiPhotonTree->Branch("gamma2",&(treeDipho_.pt2),treePho2);
  DiPhotonTree->Branch("mcgamma1",&(treeDipho_.genmatch1),treeMcPho1);
  DiPhotonTree->Branch("mcgamma2",&(treeDipho_.genmatch2),treeMcPho2);
  DiPhotonTree->Branch("vertex",&(treeDipho_.vtxX),treeVtx);
  DiPhotonTree->Branch("genVertex",&(treeDipho_.genVtxX),treeGenVtx);
}

void DiPhoAnalyzer::endJob() {

}

void DiPhoAnalyzer::initTreeStructure() {

  treeDipho_.run   = -500;
  treeDipho_.event = -500;
  treeDipho_.lumi  = -500;
  treeDipho_.nvtx  = -500;
  treeDipho_.rho   = -500.;
  treeDipho_.sampleID  = -500;
  treeDipho_.totXsec   = -500.;
  treeDipho_.pu_weight = -500.; 
  treeDipho_.pu_n = -500.;
  treeDipho_.ptgg = -500.;
  treeDipho_.mgg  = -500.;
  treeDipho_.pt1  = -500.;
  treeDipho_.ptOverM1 = -500.;
  treeDipho_.eta1 = -500.;
  treeDipho_.phi1 = -500.;
  treeDipho_.r91  = -500.;
  treeDipho_.sieie1 = -500.;
  treeDipho_.hoe1   = -500.;
  treeDipho_.scRawEne1 = -500.;
  treeDipho_.chiso1  = -500.;
  treeDipho_.phoiso1 = -500.;
  // treeDipho_.neuiso1 = -500.;
  treeDipho_.pt2  = -500.;
  treeDipho_.ptOverM2 = -500.;
  treeDipho_.eta2 = -500.;
  treeDipho_.phi2 = -500.;
  treeDipho_.r92  = -500.;
  treeDipho_.sieie2 = -500.;
  treeDipho_.hoe2   = -500.;
  treeDipho_.scRawEne2 = -500.;
  treeDipho_.chiso2  = -500.;
  treeDipho_.phoiso2 = -500.;
  // treeDipho_.neuiso2 = -500.;
  treeDipho_.vtxX = -500.;
  treeDipho_.vtxY = -500.;
  treeDipho_.vtxZ = -500.;
  treeDipho_.genmatch1 = -500;
  treeDipho_.genmatch2 = -500;
  treeDipho_.genVtxX = -500.;
  treeDipho_.genVtxY = -500.;
  treeDipho_.genVtxZ = -500.;
}

void DiPhoAnalyzer::SetPuWeights(std::string puWeightFile) {

  if (puWeightFile == "") {
    std::cout << "you need a weights file to use this function" << std::endl;
    return;
  }
  std::cout << "PU REWEIGHTING:: Using file " << puWeightFile << std::endl;

  TFile *f_pu  = new TFile(puWeightFile.c_str(),"READ");
  f_pu->cd();

  TH1D *puweights = 0;
  TH1D *gen_pu = 0;
  gen_pu    = (TH1D*) f_pu->Get("generated_pu");
  puweights = (TH1D*) f_pu->Get("weights");

  if (!puweights || !gen_pu) {
    std::cout << "weights histograms  not found in file " << puWeightFile << std::endl;
    return;
  }
  TH1D* weightedPU= (TH1D*)gen_pu->Clone("weightedPU");
  weightedPU->Multiply(puweights);

  // Rescaling weights in order to preserve same integral of events                               
  TH1D* weights = (TH1D*)puweights->Clone("rescaledWeights");
  weights->Scale( gen_pu->Integral(1,MAX_PU_REWEIGHT) / weightedPU->Integral(1,MAX_PU_REWEIGHT) );

  float sumPuWeights=0.;
  for (int i = 0; i<MAX_PU_REWEIGHT; i++) {
    float weight=1.;
    weight=weights->GetBinContent(i+1);
    sumPuWeights+=weight;
    puweights_.push_back(weight);
  }
}

float DiPhoAnalyzer::GetPUWeight(float pun) {
  
  float weight=1;
  if (sampleIndex_!=0 && pun<MAX_PU_REWEIGHT && puweights_.size()>0 && dopureweight_) 
    weight = puweights_[pun];
  return weight;
}

DEFINE_FWK_MODULE(DiPhoAnalyzer);
