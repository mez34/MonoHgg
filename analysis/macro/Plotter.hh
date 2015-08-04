#ifndef _plottertools_
#define _plottertools_

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

class Plotter{
public:
  Plotter(const TString inName, const TString outName, const TString inSpecies, const Double_t lumi);

  void DoPlots();  
  void getTree();
  void make1DHistos(); 
  void make2DHistos();
  void FindMinAndMax(TH1F *& h, int plotLog);
  void DrawWriteSave1DPlot(TH1F *& h, TString plotName, Bool_t DrawNorm); 
  void DrawWriteSave2DPlot(TH2F *& h, TString varX, TString varY); 

  ~Plotter();

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

  Int_t		NVARIABLES;
  Int_t		N2DVARIABLES;
  Int_t		nphotons;

  Float_t	variable[nvar];
  Int_t		intvariable[nvar];
  TString 	varname[nvar];
  Int_t		nbins[nvar];
  Int_t		range[nvar][2];

};

#endif
