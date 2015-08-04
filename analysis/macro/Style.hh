#ifndef _style_
#define _style_

#include "TROOT.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TString.h"

#include <iostream>

class Style{
public:
  Style(const Double_t inLumi);
  void CMSLumi(TCanvas *& canvas, const Int_t iPosX);
  void MakeOutDirectory(TString outdir);
  ~Style();

private:
  TCanvas *	canv; 
  Double_t	lumi;
};
#endif
