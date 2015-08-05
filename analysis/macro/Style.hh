#ifndef _style_
#define _style_

#include "TROOT.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TString.h"
#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"

#include <math.h>
#include <iostream>

  void MakeOutDirectory(TString outdir);
  void CheckValidFile(TFile *& file, const TString fname);
  void CheckValidTree(TTree*& tree, const TString tname, const TString fname);
  void CheckValidTH1D(TH1D*& plot, const TString pname, const TString fname);
  void CMSLumi(TCanvas *& canvas, const Int_t iPosX, const Double_t inlumi);
  void setTDRStyle(); 

class Style{
public:
  Style(const Double_t inLumi);
  ~Style();

private:
  TCanvas *	canv; 
  Double_t	lumi;
};
#endif
