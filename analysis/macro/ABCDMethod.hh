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

//#include "RooGlobalFunc.h"
#include "RooRealVar.h"

#include <iostream>
#include <fstream>
#include <cmath>

using namespace RooFit;

typedef std::vector<DblVec>   	 DblVecVec;
typedef std::vector<TFile*>   	 TFileVec;
typedef std::vector<TH1D*>    	 TH1DVec;
typedef std::vector<TH1DVec>  	 TH1DVecVec;
typedef std::vector<TH2D*>    	 TH2DVec;
typedef std::vector<TH2DVec>  	 TH2DVecVec;
typedef std::vector<RooRealVar*> RooVec;
typedef std::vector<RooVec>	 RooVecVec;

class ABCDMethod{
public: 
  ABCDMethod(const SamplePairVec Samples, const Double_t inLumi, const TString outname, Bool_t Blind);
  void DoAnalysis();
  Double_t ComputeIntAndErr(TH2D *& h, Double_t & error, const UInt_t minX, const UInt_t maxX, const UInt_t minY, const UInt_t maxY);
  void GetFinalValuesForABCDReg();
  void DoABCDCalculations();
  Double_t FindDiff(const Double_t NA, const Double_t NB, const Double_t NC, const Double_t ND);
  Double_t FindExpectedValuesInD(const Double_t NA, const Double_t NB, const Double_t NC, const Double_t NAerr, const Double_t NBerr, const Double_t NCerr, Double_t & NDerr);
  void SetRooVariables();
  void FillTable();
  void WriteDataCard(const TString fSampleName, const RooRealVar* sigrate, const Double_t expsig, const DblVecVec bkgrates, const RooVecVec bkgrate);
  void InitHists();
  void InitVariables();
  ~ABCDMethod();

private:
  //indices to keep track of bkg samples
  UInt_t	i_vh;
  UInt_t 	i_hgg;
  UInt_t	i_qcd;
  UInt_t	i_gg;
  UInt_t	i_gj;
  UInt_t	i_dy;

  Double_t	mgg_minAB1;
  Double_t	mgg_minCD;
  Double_t	mgg_maxCD;
  Double_t	mgg_maxAB2; 
  Double_t	met_minB;
  Double_t	met_minD;
  Double_t	met_maxD;

  TStrMap	fSampleTitleMap;

  RooVecVec	fRooData;
  RooVecVec	fRooBkg;
  RooVecVec	fRooSig;

  DblVecVec 	Sig_Int;
  DblVecVec 	Sig_IntErr;
  DblVecVec 	Bkg_Int;
  DblVecVec 	Bkg_IntErr;
  DblVecVec 	Data_Int;
  DblVecVec 	Data_IntErr; 

  DblVecVec 	fSig_Int;
  DblVecVec 	fSig_IntErr;
  DblVecVec 	fBkg_Int;
  DblVecVec 	fBkg_IntErr;
  DblVecVec 	fData_Int;
  DblVecVec 	fData_IntErr; 

  DblVec	fCorrData;
  DblVec	fCorrBkg;
  DblVec	fCorrSig;

  DblVec	fExpData;
  DblVec	fExpBkg;
  DblVec	fExpSig;
  DblVec	fExpErrData;
  DblVec	fExpErrBkg;
  DblVec	fExpErrSig;

  DblVec	fDiffData;
  DblVec	fDiffBkg;
  DblVec	fDiffSig;

  Double_t	lumi;
  Bool_t	doBlind;
  TString	fInDir;
  TString	fOutDir;
  TFile *	fOutFile;
  std::ofstream	fOutTxtFile;
  std::ofstream	fOutTableTxtFile;

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
  TH2DVec	fOutSelBkgTH2DHists;

};
#endif
