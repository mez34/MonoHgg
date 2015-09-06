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

typedef std::vector<DblVec>   DblVecVec;
typedef std::vector<TFile*>   TFileVec;
typedef std::vector<TH1D*>    TH1DVec;
typedef std::vector<TH1DVec>  TH1DVecVec;
typedef std::vector<TH2D*>    TH2DVec;
typedef std::vector<TH2DVec>  TH2DVecVec;

class ABCDMethod{
public: 
  ABCDMethod(const SamplePairVec Samples, const Double_t inLumi, const TString outname);
  void DoAnalysis();
  Double_t ComputeIntAndErr(TH2D *& h, Double_t & error, const Double_t minX, const Double_t maxX, const Double_t minY, const Double_t maxY, const UInt_t isReg);
  void FillTable( const TString fSampleName, const UInt_t reg, const UInt_t Integral, const UInt_t Error);
  void InitHists();
  void InitVariables();
  ~ABCDMethod();

private:
  Double_t	mgg_minAB1;
  Double_t	mgg_minCD;
  Double_t	mgg_maxCD;
  Double_t	mgg_maxAB2; 
  Double_t	met_minB;
  Double_t	met_minD;
  Double_t	met_maxD;

  DblVecVec 	Sig_Int;
  DblVecVec 	Sig_IntErr;
  DblVecVec 	Bkg_Int;
  DblVecVec 	Bkg_IntErr;
  DblVecVec 	Data_Int;
  DblVecVec 	Data_IntErr; 

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

  TH2DVec	fOutDataTH2DHists; 
  TH2DVec	fOutBkgTH2DHists;
 

};
#endif
