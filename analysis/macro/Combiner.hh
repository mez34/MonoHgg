#ifndef _combinertools_
#define _combinertools_

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
#include "THStack.h"

#include <iostream>
#include <cmath>

typedef std::pair<TString, Int_t> SamplePair;
typedef std::vector<SamplePair> SamplePairVec;
typedef SamplePairVec::iterator SamplePairVecIter;

typedef std::map<TString, Color_t> ColorMap;
typedef std::map<TString, TString> TStrMap;

typedef std::vector<TString>  TStrVec;
typedef std::vector<TFile*>   TFileVec;
typedef std::vector<TH1D*>    TH1DVec;
typedef std::vector<THStack*> THStackVec;
typedef std::vector<TLegend*> TLegVec;
typedef std::vector<TCanvas*> TCanvVec;

class Combiner{
public:
  Combiner(const SamplePairVec Samples, const Double_t inLumi, const ColorMap colorMap, const TString outname);
  void InitTH1DNames();
  void DoComb();
  void OverlayPlots();
  void StackPlots();
  ~Combiner();

private:
  Double_t 	lumi;

  UInt_t	fNData;
  UInt_t	fNBkg;
  UInt_t	fNSig;

  UInt_t	fNTH1D;
  TStrVec	fTH1DNames;

  TStrVec 	fDataNames; 
  TStrVec	fBkgNames;
  TStrVec	fSigNames;

  TFileVec	fDataFiles;
  TFileVec	fBkgFiles;
  TFileVec	fSigFiles;

  ColorMap	fColorMap;
  TStrMap	fSampleTitleMap;

  TH1DVec	fInDataTH1DHists;
  TH1DVec	fInBkgTH1DHists;
  TH1DVec	fInSigTH1DHists;

  TH1DVec	fOutBkgTH1DHists;
  TH1DVec	fOutSigTH1DHists;
  TH1DVec	fOutDataTH1DHists;

  TH1DVec	fOutTH1DComb;
  THStackVec    fOutTH1DStacks;
  TLegVec	fTH1DLegends;
  TCanvVec	fOutTH1DCanvases;

  TString	fOutDir;
  TString	fOutName;
  TFile *	fOutFile;


};
#endif
