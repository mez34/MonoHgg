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


typedef std::vector<TString> 		TStrVec;
typedef std::pair<TString, Int_t>	SelPair;

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


};

#endif
