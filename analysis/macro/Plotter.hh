#ifndef _plottertools_
#define _plottertools_

#include "Style.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TROOT.h"
#include "TMath.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLatex.h"

#include <iostream>
#include <cmath>

#define nvar 30 //nvar has to be equal to NVARIABLES

struct PlotOptStruct{
public:
  Int_t		nbins;
  Double_t	xmin;
  Double_t	xmax; 
};

typedef std::vector<TString> 		TStrVec;
typedef std::map<TString,TH1D*>		TH1DMap;
typedef TH1DMap::iterator		TH1DMapIter;

class Plotter{
public:
  Plotter(const TString inName, const TString outName, const TString inSpecies, const Double_t lumi);
  ~Plotter();

  void DoPlots();  
  void getTree();
  void make1DHistos(); 
  void make2DHistos();
  void FindMinAndMax(TH1F *& h, int plotLog);
  void DrawWriteSave1DPlot(TH1F *& h, TString plotName, Bool_t DrawNorm); 
  void DrawWriteSave2DPlot(TH2F *& h, TString varX, TString varY); 

  void SetBranchAddresses();
  void SetUpPlots();
  TH1D * MakeTH1DPlot(const TString hname, const TString htitle, const Int_t nbins, const Double_t xlow, const Double_t xhigh, const TString xtitle, const TString ytitle);
  void DoAnalysis();
  void SavePlots(); 

  void DeleteBranches();  
  void DeleteHists();

  void InitTreeVar();
  void InitTreeEffVar();
  void InitPhotonIDSel();

private:
  TString 	name;
  TString 	fName;
  TString 	species;
  TFile * 	inFile;
  TFile * 	outFile;
  TCanvas * 	fTH1Canv;
  TCanvas * 	fTH2Canv;

  Double_t 	fLumi;

  TTree * 	tpho;
  Int_t		nphotons;

  TH1DMap	fTH1DMap;

//
  TStrVec 	varname;
  TStrVec	xaxisLabel;
  Int_t		NVARIABLES;
  Float_t	variable[nvar];

  TStrVec       effvar;	
  Int_t		N2DVARIABLES;
  Int_t		intvariable[nvar];

  TStrVec	selvar;
  Int_t		NSEL;
  Int_t		selvarPair[10];

  Int_t		nbins[nvar];
  Int_t		range[nvar][2];
//

  // variables for branches
  Int_t 	nvtx;
  Float_t	weight;
  Float_t	mgg;
  Float_t	ptgg;
  Float_t	t1pfmet; 
  Float_t	t1pfmetphi; 
  Float_t	pt1;
  Float_t	pt2;
  Float_t	phi1;
  Float_t	phi2;
  Float_t	eta1;
  Float_t	eta2;
  Float_t	r91;
  Float_t	r92;
  Float_t	phoiso1;
  Float_t	phoiso2;
  Float_t	chiso1;
  Float_t	chiso2;
  Float_t	neuiso1;
  Float_t	neuiso2;
  Float_t	sieie1;
  Float_t	sieie2;
  Float_t	hoe1;
  Float_t	hoe2;
  Int_t		passCHiso1;
  Int_t		passCHiso2;
  Int_t		passNHiso1;
  Int_t		passNHiso2;
  Int_t		passPHiso1;
  Int_t		passPHiso2;
  Int_t		passSieie1;
  Int_t		passSieie2;
  Int_t		passHoe1;
  Int_t		passHoe2;
 
  // branches
  TBranch 	*b_nvtx;
  TBranch	*b_weight;
  TBranch	*b_mgg;
  TBranch	*b_ptgg;
  TBranch	*b_t1pfmet;
  TBranch	*b_t1pfmetPhi;
  TBranch	*b_pt1;
  TBranch	*b_pt2;
  TBranch	*b_phi1;
  TBranch	*b_phi2;
  TBranch	*b_eta1;
  TBranch	*b_eta2;
  TBranch	*b_r91;
  TBranch	*b_r92;
  TBranch	*b_phoiso1;
  TBranch	*b_phoiso2;
  TBranch	*b_chiso1;
  TBranch	*b_chiso2;
  TBranch	*b_neuiso1;
  TBranch	*b_neuiso2;
  TBranch	*b_sieie1;
  TBranch	*b_sieie2;
  TBranch	*b_hoe1;
  TBranch	*b_hoe2;
  TBranch	*b_passCHiso1;
  TBranch	*b_passCHiso2;
  TBranch	*b_passNHiso1;
  TBranch	*b_passNHiso2;
  TBranch	*b_passPHiso1;
  TBranch	*b_passPHiso2;
  TBranch	*b_passSieie1;
  TBranch	*b_passSieie2;
  TBranch	*b_passHoe1;
  TBranch	*b_passHoe2;
};

#endif
