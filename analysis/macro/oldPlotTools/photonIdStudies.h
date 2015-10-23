//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Dec  9 12:22:26 2014 by ROOT version 5.34/18
// from TTree singlePhotons/single photon tree
// found on file: ../python/singlePhotonTree__kMpl01_M3000.root
//////////////////////////////////////////////////////////

#ifndef photonIdStudies_h
#define photonIdStudies_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class photonIdStudies {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
   // Declaration of leaf types
  Float_t         kinematics_pt;
  Float_t         kinematics_eta;
  Float_t         kinematics_phi;
  Int_t           kinematics_isEBEtaGap;
  Int_t           kinematics_isEBPhiGap;
  Int_t           kinematics_isEERingGap;
  Int_t           kinematics_isEEDeeGap;
  Int_t           kinematics_isEBEEGap;
  Float_t         supercluster_scEta;
  Float_t         supercluster_scPhi;
  Float_t         energy_eMax;
  Float_t         energy_e5x5;
  Float_t         energy_energy;
  Float_t         energy_energyInitial;
  Float_t         energy_energyRegression;
  Float_t         identification_e1x5;
  Float_t         identification_e2x5;
  Float_t         identification_sigmaIetaIeta;
  Float_t         identification_r9;
  Float_t         identification_hoe;
  Float_t         identification_h1oe;
  Float_t         identification_h2oe;
  Float_t         identification_htoe;
  Float_t         identification_ht1oe;
  Float_t         identification_ht2oe;
  Char_t          identification_passEleVeto;
  Char_t          identification_hasPixelSeed;
  Float_t         identificationNoZS_sieienoZS;
  Float_t         identificationNoZS_e5x5noZS;
  Float_t         identificationNoZS_e1x5noZS;
  Float_t         identificationNoZS_r9noZS;
  Float_t         isolation_trackIso;
  Float_t         isolation_ecalIso;
  Float_t         isolation_hcalIso;
  Float_t         isolation_chHadIso;
  Float_t         isolation_nHadIso;
  Float_t         isolation_photonIso;
  Float_t         isolation_rho;
  Float_t         mctruth_trueEnergy;
  Float_t         mctruth_truePt;
  Float_t         mctruth_trueEta;
  Float_t         mctruth_truePhi;
  Float_t         mctruth_minDR;
  Float_t         tree5x5_amplit[25];
  Int_t           tree5x5_ieta[25];
  Int_t           tree5x5_iphi[25];
  Int_t           tree5x5_ix[25];
  Int_t           tree5x5_iy[25];
  Int_t           tree5x5_iz[25];
  Int_t           tree5x5_kSaturated[25];
  Int_t           tree5x5_kLeRecovered[25];
  Int_t           tree5x5_kNeighRecovered[25];
  
  // List of branches
  TBranch        *b_kinematics;   //!
  TBranch        *b_supercluster;   //!
  TBranch        *b_energy;   //!
  TBranch        *b_identification;   //!
  TBranch        *b_identificationNoZS;   //!      
  TBranch        *b_isolation;   //!
  TBranch        *b_mctruth;   //!
  TBranch        *b_tree5x5;   //!
  
  photonIdStudies(TTree *tree=0);
  virtual ~photonIdStudies();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef photonIdStudies_cxx
photonIdStudies::photonIdStudies(TTree *tree) : fChain(0) 
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    // TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../python/singlePhotonTree__kMpl01_M3000.root");
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../python/singlePhotonTree_10000gj.root");
    if (!f || !f->IsOpen()) {
      // f = new TFile("../python/singlePhotonTree__kMpl01_M3000.root");
      f = new TFile("../python/singlePhotonTree_10000gj.root");
    }
    // TDirectory * dir = (TDirectory*)f->Get("../python/singlePhotonTree__kMpl01_M3000.root:/singlePhoAna");
    TDirectory * dir = (TDirectory*)f->Get("../python/singlePhotonTree_10000gj.root:/singlePhoAna");
    dir->GetObject("singlePhotons",tree);
    
  }
  Init(tree);
}

photonIdStudies::~photonIdStudies()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t photonIdStudies::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Long64_t photonIdStudies::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void photonIdStudies::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).
  
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
  
  fChain->SetBranchAddress("kinematics", &kinematics_pt, &b_kinematics);
  fChain->SetBranchAddress("supercluster", &supercluster_scEta, &b_supercluster);
  fChain->SetBranchAddress("energy", &energy_eMax, &b_energy);
  fChain->SetBranchAddress("identification", &identification_e1x5, &b_identification);
  fChain->SetBranchAddress("identificationNoZS", &identificationNoZS_sieienoZS, &b_identificationNoZS);
  fChain->SetBranchAddress("isolation", &isolation_trackIso, &b_isolation);
  fChain->SetBranchAddress("mctruth", &mctruth_trueEnergy, &b_mctruth);
  fChain->SetBranchAddress("tree5x5", tree5x5_amplit, &b_tree5x5);
  Notify();
}

Bool_t photonIdStudies::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.
  
  return kTRUE;
}

void photonIdStudies::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

Int_t photonIdStudies::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}

#endif // #ifdef photonIdStudies_cxx
