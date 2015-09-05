#ifndef _ABCDfunction_
#define _ABCDfunction_

#include "Style.hh"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
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

typedef std::vector<TFile*>   TFileVec;
typedef std::vector<TH1D*>    TH1DVec;
typedef std::vector<TH1DVec>  TH1DVecVec;
typedef std::vector<TH2D*>    TH2DVec;
typedef std::vector<TH2DVec>  TH2DVecVec;

class ABCDMethod{
public: 
  ABCDMethod(const SamplePairVec Samples, const Double_t inLumi, const TString outname);
  void DoAnalysis();
  void InitHists();
  void InitVariables();
  ~ABCDMethod();

private:
  Double_t	lumi;
  TString	fInDir;
  TString	fOutDir;
  TFile *	fOutFile;

  Int_t		fNBkg;
  Int_t 	fNSig;
  Int_t		fNData;
  Int_t		fNTH1D;
  Int_t		fNTH2D;

  TStrVec	fTH1DNames;
  TStrVec	fTH2DNames;

  TStrVec	fSigNames;
  TStrVec	fBkgNames;
  TStrVec	fDataNames;

  TFileVec	fDataFiles;
  TFileVec	fBkgFiles;
  TFileVec	fSigFiles;

  TH1DVecVec	fInDataTH1DHists;
  TH1DVecVec	fInBkgTH1DHists;
  TH1DVecVec	fInSigTH1DHists;

  TH2DVecVec	fInDataTH2DHists;
  TH2DVecVec	fInBkgTH2DHists;
  TH2DVecVec	fInSigTH2DHists;

  TH1DVec	fOutDataTH1DHists;
  TH1DVec	fOutBkgTH1DHists;
 

};
#endif
