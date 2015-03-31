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
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

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
  float neuiso1;
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
  float neuiso2;
  int vtxIndex;
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
  bool isGammaPresel( float pt, float r9, float chiso);
  bool isGammaSelected( float rho, float pt, float sceta, float r9, float chiso, float nhiso, float phoiso, float hoe, float sieie);
  int effectiveAreaRegion(float sceta);

  // collections
  edm::EDGetTokenT<EcalRecHitCollection> ecalHitEBToken_;
  edm::EDGetTokenT<EcalRecHitCollection> ecalHitEEToken_;
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

  // events breakdown
  TH1F *h_selection;
};
   

DiPhoAnalyzer::DiPhoAnalyzer(const edm::ParameterSet& iConfig):

  // collections
  ecalHitEBToken_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedBarrelRecHitCollection"))),
  ecalHitEEToken_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedEndcapRecHitCollection"))),
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


  // to recompute not-zero-suppressed cluster shapes 
  noZS::EcalClusterLazyTools *lazyToolnoZS = new noZS::EcalClusterLazyTools(iEvent, iSetup, ecalHitEBToken_, ecalHitEEToken_);

  // To keep track of the total number of events
  h_entries->Fill(5);

  // Events breakdown
  h_selection->Fill(0);

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
  if (sampleID>0) {   
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
  if (sampleID>0) totXsec = xsec_ * kfac_;

  
  // analysis cuts: trigger (chiara: qualcosa per emulare?)
  
  // Diphoton candidates before any selection
  if (dipho.size()>0) {

    // Diphoton candidates: preselection
    vector<int> preselDipho;
    for (unsigned int diphotonlooper = 0; diphotonlooper < dipho.size() ; diphotonlooper++){
      
      float leadPt     = dipho[diphotonlooper]->leadingPhoton()->et();
      float leadChIso  = dipho[diphotonlooper]->leadingPhoton()->egChargedHadronIso();
      float leadR9noZS = lazyToolnoZS->e3x3(*(dipho[diphotonlooper]->leadingPhoton()->superCluster()->seed())) / dipho[diphotonlooper]->leadingPhoton()->superCluster()->rawEnergy();
      bool leadPresel  = isGammaPresel( leadPt, leadR9noZS, leadChIso); 
      
      float subleadPt     = dipho[diphotonlooper]->subLeadingPhoton()->et();
      float subleadChIso  = dipho[diphotonlooper]->subLeadingPhoton()->egChargedHadronIso();
      float subleadR9noZS = lazyToolnoZS->e3x3(*(dipho[diphotonlooper]->subLeadingPhoton()->superCluster()->seed())) / dipho[diphotonlooper]->subLeadingPhoton()->superCluster()->rawEnergy();
      bool subleadPresel  = isGammaPresel( subleadPt, subleadR9noZS, subleadChIso); 
      if (!leadPresel || !subleadPresel) continue;
      
      preselDipho.push_back(diphotonlooper);
    }

    if (preselDipho.size()>0) {
      h_selection->Fill(1);  
  
      
      // Diphoton candidates: Id/isolation selection
      vector<int> selectedDipho;
      for (unsigned int diphotonlooper =0; diphotonlooper < preselDipho.size() ; diphotonlooper++){
	
	int theDiphoton = preselDipho[diphotonlooper];
	
	std::vector<float> leadCovnoZS    = lazyToolnoZS->localCovariances(*(dipho[theDiphoton]->leadingPhoton()->superCluster()->seed())) ;
	std::vector<float> subleadCovnoZS = lazyToolnoZS->localCovariances(*(dipho[theDiphoton]->subLeadingPhoton()->superCluster()->seed())) ;
	
	float leadSieienoZS;
	if (!isnan(leadCovnoZS[0]))
	  leadSieienoZS = sqrt (leadCovnoZS[0]); 
	else 
	  continue;
	
	float subleadSieienoZS;
	if (!isnan(subleadCovnoZS[0]))
	  subleadSieienoZS = sqrt (subleadCovnoZS[0]); 
	else 
	  continue;
	float leadPt     = dipho[theDiphoton]->leadingPhoton()->et();
	float leadScEta  = (dipho[theDiphoton]->leadingPhoton()->superCluster())->eta();   
	float leadR9noZS = lazyToolnoZS->e3x3(*(dipho[theDiphoton]->leadingPhoton()->superCluster()->seed())) / dipho[theDiphoton]->leadingPhoton()->superCluster()->rawEnergy();
	float leadHoE    = dipho[theDiphoton]->leadingPhoton()->hadronicOverEm();
	float leadChIso  = dipho[theDiphoton]->leadingPhoton()->egChargedHadronIso();
	float leadNeuIso = dipho[theDiphoton]->leadingPhoton()->egNeutralHadronIso();
	float leadPhoIso = dipho[theDiphoton]->leadingPhoton()->egPhotonIso();
	bool  leadSelel  = isGammaSelected( rho, leadPt, leadScEta, leadR9noZS, leadChIso, leadNeuIso, leadPhoIso, leadHoE, leadSieienoZS); 
	
	float subleadPt     = dipho[theDiphoton]->subLeadingPhoton()->et();
	float subleadScEta  = (dipho[theDiphoton]->subLeadingPhoton()->superCluster())->eta();   
	float subleadR9noZS = lazyToolnoZS->e3x3(*(dipho[theDiphoton]->subLeadingPhoton()->superCluster()->seed())) / dipho[theDiphoton]->subLeadingPhoton()->superCluster()->rawEnergy();
	float subleadHoE    = dipho[theDiphoton]->subLeadingPhoton()->hadronicOverEm();
	float subleadChIso  = dipho[theDiphoton]->subLeadingPhoton()->egChargedHadronIso();
	float subleadNeuIso = dipho[theDiphoton]->subLeadingPhoton()->egNeutralHadronIso();
	float subleadPhoIso = dipho[theDiphoton]->subLeadingPhoton()->egPhotonIso();
	bool  subleadSelel  = isGammaSelected( rho, subleadPt, subleadScEta, subleadR9noZS, subleadChIso, subleadNeuIso, subleadPhoIso, subleadHoE, subleadSieienoZS);  
	if (!leadSelel || !subleadSelel) continue;
	
	selectedDipho.push_back(theDiphoton);    
      }

      if (selectedDipho.size()>0) {
	h_selection->Fill(2);  


	// Diphoton candidates: pT
	vector<int> kineDipho;
	for (unsigned int diphotonlooper =0; diphotonlooper < selectedDipho.size() ; diphotonlooper++){
	  
	  int theDiphoton = selectedDipho[diphotonlooper];
	  
	  float leadPt     = dipho[theDiphoton]->leadingPhoton()->et();
	  float subleadPt  = dipho[theDiphoton]->subLeadingPhoton()->et();
	  if (leadPt<80 || subleadPt<80) continue;           // chiara: hardcoded
	  
	  kineDipho.push_back(theDiphoton);
	}
	
	if (kineDipho.size()>0) {
	  h_selection->Fill(3);
  
	  
	  // Diphoton candidates choice
	  // chiara: per il momento quello a pT del sistema piu' alto
	  float maxDiphoPt = -999.;
	  int candIndex = 9999;    // This int will store the index of the best diphoton candidate
	  for (unsigned int diphotonlooper =0; diphotonlooper < kineDipho.size() ; diphotonlooper++){
	    
	    int theDiphoton = kineDipho[diphotonlooper];
	    
	    float thisSystemPt = dipho[theDiphoton]->pt();
	    if (thisSystemPt>maxDiphoPt) {
	      maxDiphoPt = thisSystemPt;
	      candIndex  = theDiphoton;
	    }
	  }
	  
	  if (candIndex<999) {
	    h_selection->Fill(4);

	    // analysis cuts: good vertex 
	    // Since diphoton candidates have already an associated vtx, I check it only and discard the event if bad 
	    // this is why I put this selection AFTER the diphoton choice
	    bool goodVtx = true;
	    int theVertex = dipho[candIndex]->vertex_index();
	    if (vtxs[theVertex]->ndof()<=4)  goodVtx = false;
	    float d0vtx = sqrt( vtxs[theVertex]->position().x()*vtxs[theVertex]->position().x() + vtxs[theVertex]->position().y()*vtxs[theVertex]->position().y() );
	    if (fabs(d0vtx)>2) goodVtx = false;
	    if (fabs(vtxs[theVertex]->position().z())>=24) goodVtx = false;
	    bool isVtxFake = (vtxs[theVertex]->ndof()==0) && (vtxs[theVertex]->chi2()==0);   // chiara: also && tracks.empty, but can not be used here
	    if (isVtxFake) goodVtx = false;

	    if (goodVtx) {
	      h_selection->Fill(5);	    

	      // to be kept in the tree
	      float ptgg, mgg;
	      float pt1, ptOverM1, eta1, phi1;
	      float r91, sieie1, hoe1, scRawEne1;
	      float chiso1, phoiso1, neuiso1;
	      float pt2, ptOverM2, eta2, phi2;
	      float r92, sieie2, hoe2, scRawEne2;
	      float chiso2, phoiso2, neuiso2;
	      int vtxIndex;
	      float vtxX, vtxY, vtxZ;
	      int genmatch1, genmatch2;
	      float genVtxX, genVtxY, genVtxZ;   // chiara: vuoti
	      
	      
	      // analysis cuts: mgg > 300  
	      if (dipho[candIndex]->mass()>=300) {
		h_selection->Fill(6);
		
		// fully selected event: tree re-initialization                                                                          
		initTreeStructure();        
		
		//-------> diphoton system properties 
		ptgg = dipho[candIndex]->pt();
		mgg  = dipho[candIndex]->mass();
		
		
		//-------> individual photon properties, cominciamo con questi
		std::vector<float> leadCovnoZS = lazyToolnoZS->localCovariances(*(dipho[candIndex]->leadingPhoton()->superCluster()->seed())) ;
		std::vector<float> subleadCovnoZS = lazyToolnoZS->localCovariances(*(dipho[candIndex]->subLeadingPhoton()->superCluster()->seed())) ;
		
		pt1       = dipho[candIndex]->leadingPhoton()->et();
		ptOverM1  = pt1/mgg;
		eta1      = dipho[candIndex]->leadingPhoton()->eta();
		phi1      = dipho[candIndex]->leadingPhoton()->phi();
		r91       = lazyToolnoZS->e3x3(*(dipho[candIndex]->leadingPhoton()->superCluster()->seed())) / dipho[candIndex]->leadingPhoton()->superCluster()->rawEnergy();
		sieie1    = sqrt(leadCovnoZS[0]);
		hoe1      = dipho[candIndex]->leadingPhoton()->hadronicOverEm();
		scRawEne1 = dipho[candIndex]->leadingPhoton()->superCluster()->rawEnergy();
		chiso1    = dipho[candIndex]->leadingPhoton()->egChargedHadronIso();
		phoiso1   = dipho[candIndex]->leadingPhoton()->egPhotonIso();
		neuiso1   = dipho[candIndex]->leadingPhoton()->egNeutralHadronIso();
		
		pt2       = dipho[candIndex]->subLeadingPhoton()->et();
		ptOverM2  = pt2/mgg;
		eta2      = dipho[candIndex]->subLeadingPhoton()->eta();
		phi2      = dipho[candIndex]->subLeadingPhoton()->phi();
		r92       = lazyToolnoZS->e3x3(*(dipho[candIndex]->subLeadingPhoton()->superCluster()->seed())) / dipho[candIndex]->subLeadingPhoton()->superCluster()->rawEnergy();;
		sieie2    = sqrt(subleadCovnoZS[0]);
		hoe2      = dipho[candIndex]->subLeadingPhoton()->hadronicOverEm();
		scRawEne2 = dipho[candIndex]->subLeadingPhoton()->superCluster()->rawEnergy();
		chiso2    = dipho[candIndex]->subLeadingPhoton()->egChargedHadronIso();
		phoiso2   = dipho[candIndex]->subLeadingPhoton()->egPhotonIso();
		neuiso2   = dipho[candIndex]->subLeadingPhoton()->egNeutralHadronIso();
		
		//-------> vtx info
		vtxIndex = dipho[candIndex]->vertex_index();
		vtxX= dipho[candIndex]->getVertex()->x();
		vtxY= dipho[candIndex]->getVertex()->y();
		vtxZ= dipho[candIndex]->getVertex()->z();
		
		//-------> generated vtx info
		// chiara: rigirare per accedere alla mc truth completa (check con Pasquale)
		// controllare, ma dovrebbe esserci associata al diphoton (ma e' vuota)
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
		genmatch1 = 0;
		genmatch2 = 0;
		// genmatch1 = (dipho[candIndex]->leadingPhoton()->genMatchType() == Photon::kPrompt); 
		// genmatch2 = (dipho[candIndex]->subLeadingPhoton()->genMatchType() == Photon::kPrompt); 
		if (sampleID>0) {   
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
		neuiso1  = -999.;
		
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
		neuiso2  = -999.;
		
		vtxIndex = -999;
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
	      treeDipho_.neuiso1 = neuiso1;
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
	      treeDipho_.neuiso2 = neuiso2;
	      treeDipho_.vtxIndex = vtxIndex;
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

	    } // final diphoton candidate with good vtx
	  }   // final diphoton candidate found
	}     // diphoton candidate passing pt cuts
      }       // diphoton candidate passing ID+iso
    }         // diphoton candidate passing preselection
  }           // at least 1 reco diphoton candidate  

  // delete
  delete lazyToolnoZS;
}

void DiPhoAnalyzer::beginJob() {

  // loading weights for pileup if needed
  if (dopureweight_) 
    SetPuWeights(puWFileName_);
  
  // to keep track of the original number of events
  h_entries = fs_->make<TH1F>("h_entries", "h_entries", 10,  0., 10.);

  // for the event breakdown
  h_selection = fs_->make<TH1F>("h_selection", "h_selection", 10, 0., 10.);
  
  // Trees
  DiPhotonTree = fs_->make<TTree>("DiPhotonTree","di-photon tree");

  // with all infos
  DiPhotonTree->Branch("run",&(treeDipho_.run),"run/I");
  DiPhotonTree->Branch("event",&(treeDipho_.event),"event/I");
  DiPhotonTree->Branch("lumi",&(treeDipho_.lumi),"lumi/I");
  DiPhotonTree->Branch("nvtx",&(treeDipho_.nvtx),"nvtx/I");
  DiPhotonTree->Branch("rho",&(treeDipho_.rho),"rho/F");
  DiPhotonTree->Branch("sampleID",&(treeDipho_.sampleID),"sampleID/I");
  DiPhotonTree->Branch("totXsec",&(treeDipho_.totXsec),"totXsec/F");
  DiPhotonTree->Branch("pu_weight",&(treeDipho_.pu_weight),"pu_weight/F");
  DiPhotonTree->Branch("pu_n",&(treeDipho_.pu_n),"pu_n/F");
  DiPhotonTree->Branch("ptgg",&(treeDipho_.ptgg),"ptgg/F");
  DiPhotonTree->Branch("mgg",&(treeDipho_.mgg),"mgg/F");
  DiPhotonTree->Branch("pt1",&(treeDipho_.pt1),"pt1/F");
  DiPhotonTree->Branch("ptOverM1",&(treeDipho_.ptOverM1),"ptOverM1/F");
  DiPhotonTree->Branch("eta1",&(treeDipho_.eta1),"eta1/F");
  DiPhotonTree->Branch("phi1",&(treeDipho_.phi1),"phi1/F");
  DiPhotonTree->Branch("r91",&(treeDipho_.r91),"r91/F");
  DiPhotonTree->Branch("sieie1",&(treeDipho_.sieie1),"sieie1/F");
  DiPhotonTree->Branch("hoe1",&(treeDipho_.hoe1),"hoe1/F");
  DiPhotonTree->Branch("scRawEne1",&(treeDipho_.scRawEne1),"scRawEne1/F");
  DiPhotonTree->Branch("chiso1",&(treeDipho_.chiso1),"chiso1/F");
  DiPhotonTree->Branch("phoiso1",&(treeDipho_.phoiso1),"phoiso1/F");
  DiPhotonTree->Branch("neuiso1",&(treeDipho_.neuiso1),"neuiso1/F");
  DiPhotonTree->Branch("pt2",&(treeDipho_.pt2),"pt2/F");
  DiPhotonTree->Branch("ptOverM2",&(treeDipho_.ptOverM2),"ptOverM2/F");
  DiPhotonTree->Branch("eta2",&(treeDipho_.eta2),"eta2/F");
  DiPhotonTree->Branch("phi2",&(treeDipho_.phi2),"phi2/F");
  DiPhotonTree->Branch("r92",&(treeDipho_.r92),"r92/F");
  DiPhotonTree->Branch("sieie2",&(treeDipho_.sieie2),"sieie2/F");
  DiPhotonTree->Branch("hoe2",&(treeDipho_.hoe2),"hoe2/F");
  DiPhotonTree->Branch("scRawEne2",&(treeDipho_.scRawEne2),"scRawEne2/F");
  DiPhotonTree->Branch("chiso2",&(treeDipho_.chiso2),"chiso2/F");
  DiPhotonTree->Branch("phoiso2",&(treeDipho_.phoiso2),"phoiso2/F");
  DiPhotonTree->Branch("neuiso2",&(treeDipho_.neuiso2),"neuiso2/F");
  DiPhotonTree->Branch("genmatch1",&(treeDipho_.genmatch1),"genmatch1/I");
  DiPhotonTree->Branch("genmatch2",&(treeDipho_.genmatch2),"genmatch12/I");
  DiPhotonTree->Branch("vtxIndex",&(treeDipho_.vtxIndex),"vtxIndex/I");
  DiPhotonTree->Branch("vtxX",&(treeDipho_.vtxX),"vtxX/F");
  DiPhotonTree->Branch("vtxY",&(treeDipho_.vtxY),"vtxY/F");
  DiPhotonTree->Branch("vtxZ",&(treeDipho_.vtxZ),"vtxZ/F");
  DiPhotonTree->Branch("genVtxX",&(treeDipho_.genVtxX),"genVtxX/F");
  DiPhotonTree->Branch("genVtxY",&(treeDipho_.genVtxY),"genVtxY/F");
  DiPhotonTree->Branch("genVtxZ",&(treeDipho_.genVtxZ),"genVtxZ/F");
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
  treeDipho_.neuiso1 = -500.;
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
  treeDipho_.neuiso2 = -500.;
  treeDipho_.vtxIndex = -500;
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

// AOD preselection
bool DiPhoAnalyzer::isGammaPresel( float pt, float r9, float chiso) {

  bool isPresel = false;
  if (r9>0.8)         isPresel = true;
  if (chiso<20)       isPresel = true;
  if ((chiso/pt)<0.3) isPresel = true;

  return isPresel;
}

// chiara: still to be fully defined
bool DiPhoAnalyzer::isGammaSelected( float rho, float pt, float sceta, float r9, float chiso, float nhiso, float phoiso, float hoe, float sieie) {
  
  // classes: 0 = EB highR9, 1 = EB low R9, 2 = EE high R9, 3 = EE lowR9
  int etaclass = fabs(sceta)>1.5;
  int r9class  = r9<0.94;                   
  int theclass = 2.*etaclass + r9class;                  

  // cuts - hardcoded
  float chiso_cut[4]  = { 5.95, 7.08, 6.10, 5.07 };
  float phoiso_cut[4] = { 2.87, 5.47, 5.98, 3.44 };
  float nhiso_cut[4]  = { 27.4, 30.0, 30.0, 15.0 };   
  float sieie_cut[4]  = { 1.05e-02, 1.05e-02, 2.8e-02, 2.8e-02 };
  float hoe_cut[4]    = { 4.53e-01, 2.10e-01, 6.3e-02, 7.8e-02 };
  
  // effective areas - hardcoded 
  float chIsoAE[5] = { 0.00,0.000,0.00,0.00,0.00 };
  float phIsoAE[5] = { 0.21,0.200,0.14,0.22,0.31 };
  float nhIsoAE[5] = { 0.04,0.059,0.05,0.05,0.15 };

  // EA corrections 
  int theEAregion = effectiveAreaRegion(sceta);
  float corrChIso = chiso - rho*chIsoAE[theEAregion];
  float corrPhIso = phoiso - rho*phIsoAE[theEAregion];
  float corrNhIso = nhiso - rho*nhIsoAE[theEAregion];   

  if (corrChIso > chiso_cut[theclass])  return false;
  if (corrPhIso > phoiso_cut[theclass]) return false;
  if (corrNhIso > nhiso_cut[theclass])  return false;
  if (sieie > sieie_cut[theclass])      return false;
  if (hoe> hoe_cut[theclass])           return false;
  
  // chiara: manca il veto elettroni 

  return true;
} 

int DiPhoAnalyzer::effectiveAreaRegion(float sceta) {

  int theEAregion = 999;
  if (fabs(sceta)<=0.9) theEAregion = 0;
  if (fabs(sceta)<=1.5 && fabs(sceta)>0.9)  theEAregion = 1;
  if (fabs(sceta)<=2.0 && fabs(sceta)>1.5)  theEAregion = 2;
  if (fabs(sceta)<=2.2 && fabs(sceta)>2.0)  theEAregion = 3;
  if (fabs(sceta)<=2.5 && fabs(sceta)>2.2)  theEAregion = 4;
  return theEAregion;
}


DEFINE_FWK_MODULE(DiPhoAnalyzer);
