#define singleGammaPlots_cxx
#include "singleGammaPlots.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TProfile.h>
#include <TGraphAsymmErrors.h>
#include <iostream.h> 

using namespace std;

int singleGammaPlots::effectiveAreaRegion(float theEta) {

  int theEAregion = 999;
  if (fabs(theEta)<1.) theEAregion = 0;
  if (fabs(theEta)<1.479 && fabs(theEta)>1.)    theEAregion = 1;
  if (fabs(theEta)<2.    && fabs(theEta)>1.479) theEAregion = 2;
  if (fabs(theEta)<2.2   && fabs(theEta)>2.0)   theEAregion = 3;
  if (fabs(theEta)<2.3   && fabs(theEta)>2.2)   theEAregion = 4;
  if (fabs(theEta)<2.4   && fabs(theEta)>2.3)   theEAregion = 5;
  if (fabs(theEta)>2.4) theEAregion = 6;
  return theEAregion;
}

float singleGammaPlots::corrChIso(float uncorrChIso, float etaSc) {

  float EA_charged[7] = { 0.012, 0.010, 0.014, 0.012, 0.016, 0.020, 0.012};
  int theEAregion_fG  = effectiveAreaRegion(etaSc);
  float corr = uncorrChIso - isolation_rho*EA_charged[theEAregion_fG];
  return corr;
}

float singleGammaPlots::corrNIso(float uncorrNIso, float etaSc) {

  float EA_neutral[7] = { 0.030, 0.057, 0.039, 0.015, 0.024, 0.039, 0.072};
  int theEAregion_fG  = effectiveAreaRegion(etaSc);
  float corr = uncorrNIso - isolation_rho*EA_neutral[theEAregion_fG];
  return corr;
}

float singleGammaPlots::corrGIso(float uncorrGIso, float etaSc) {

  float EA_photons[7] = { 0.148, 0.130, 0.112, 0.216, 0.262, 0.260, 0.266};
  int theEAregion_fG  = effectiveAreaRegion(etaSc);
  float corr = uncorrGIso - isolation_rho*EA_photons[theEAregion_fG];
  return corr;
}


void singleGammaPlots::Loop()
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;

  // counters for efficiencies
  float totalEB = 0;
  float totalEB_noGap = 0;
  float totalEB_noGap_highR9 = 0;
  float okHoE_totalEB = 0;
  float okHoE_totalEB_noGap = 0;
  float okHoE_totalEB_noGap_highR9 = 0;
  float okHoET_totalEB = 0;
  float okHoET_totalEB_noGap = 0;
  float okHoET_totalEB_noGap_highR9 = 0;
  //
  float totalEE = 0;
  float totalEE_noGap = 0;
  float totalEE_noGap_highR9 = 0;
  float okHoE_totalEE = 0;
  float okHoE_totalEE_noGap = 0;
  float okHoE_totalEE_noGap_highR9 = 0;
  float okHoET_totalEE = 0;
  float okHoET_totalEE_noGap = 0;
  float okHoET_totalEE_noGap_highR9 = 0;



  // booking histos - kinematics
  TH1F *H_pt_EB = new TH1F("H_pt_EB","p_{T}, barrel",90, 0.,4500.);
  TH1F *H_pt_EE = new TH1F("H_pt_EE","p_{T}, endcap",90, 0.,4500.);
  H_pt_EB->Sumw2();
  H_pt_EE->Sumw2();
  H_pt_EB->GetXaxis()->SetTitle("pT [GeV]");
  H_pt_EE->GetXaxis()->SetTitle("pT [GeV]");

  TH1F *H_eta = new TH1F("H_eta","#eta",50, -2.6,2.6);
  H_eta->Sumw2();
  H_eta->GetXaxis()->SetTitle("#eta");


  
  // booking histos - energies
  TH1F *H_E25OverEtrue_EB = new TH1F("H_E25OverEtrue_EB","H_E25OverEtrue_EB",50, 0.8,1.1);
  TH1F *H_E25OverEtrue_EE = new TH1F("H_E25OverEtrue_EE","H_E25OverEtrue_EE",50, 0.8,1.1);
  H_E25OverEtrue_EB->Sumw2();
  H_E25OverEtrue_EE->Sumw2();
  H_E25OverEtrue_EB->SetTitle("");
  H_E25OverEtrue_EE->SetTitle("");
  H_E25OverEtrue_EB->GetXaxis()->SetTitle("E5x5 / Etrue");
  H_E25OverEtrue_EE->GetXaxis()->SetTitle("E5x5 / Etrue");

  TH1F *H_EneOverEtrue_EB = new TH1F("H_EneOverEtrue_EB","H_EneOverEtrue_EB",50, 0.8,1.1);
  TH1F *H_EneOverEtrue_EE = new TH1F("H_EneOverEtrue_EE","H_EneOverEtrue_EE",50, 0.8,1.1);
  H_EneOverEtrue_EB->Sumw2();
  H_EneOverEtrue_EE->Sumw2();
  H_EneOverEtrue_EB->SetTitle("");
  H_EneOverEtrue_EE->SetTitle("");
  H_EneOverEtrue_EB->GetXaxis()->SetTitle("Ene / Etrue");
  H_EneOverEtrue_EE->GetXaxis()->SetTitle("Ene / Etrue");

  TProfile *P_E25OverEtrueVsETtrue_EB = new TProfile("P_E25OverEtrueVsETtrue_EB","P_E25OverEtrueVsETtrue_EB",90, 0.,4500.);
  TProfile *P_E25OverEtrueVsETtrue_EE = new TProfile("P_E25OverEtrueVsETtrue_EE","P_E25OverEtrueVsETtrue_EE",90, 0.,4500.);
  P_E25OverEtrueVsETtrue_EB->SetTitle("");
  P_E25OverEtrueVsETtrue_EB->SetMarkerStyle(20);
  P_E25OverEtrueVsETtrue_EB->SetMarkerSize(0.8);
  P_E25OverEtrueVsETtrue_EB->SetMarkerColor(2);
  P_E25OverEtrueVsETtrue_EB->GetXaxis()->SetTitle("ET true");
  P_E25OverEtrueVsETtrue_EB->GetYaxis()->SetTitle("E25 / Etrue");
  P_E25OverEtrueVsETtrue_EE->SetTitle("");
  P_E25OverEtrueVsETtrue_EE->SetMarkerStyle(20);
  P_E25OverEtrueVsETtrue_EE->SetMarkerSize(0.8);
  P_E25OverEtrueVsETtrue_EE->SetMarkerColor(2);
  P_E25OverEtrueVsETtrue_EE->GetXaxis()->SetTitle("ET true");
  P_E25OverEtrueVsETtrue_EE->GetYaxis()->SetTitle("E25 / Etrue");

  TProfile *P_EneOverEtrueVsETtrue_EB = new TProfile("P_EneOverEtrueVsETtrue_EB","P_EneOverEtrueVsETtrue_EB",90, 0.,4500.);
  TProfile *P_EneOverEtrueVsETtrue_EE = new TProfile("P_EneOverEtrueVsETtrue_EE","P_EneOverEtrueVsETtrue_EE",90, 0.,4500.);
  P_EneOverEtrueVsETtrue_EB->SetTitle("");
  P_EneOverEtrueVsETtrue_EB->SetMarkerStyle(20);
  P_EneOverEtrueVsETtrue_EB->SetMarkerSize(0.8);
  P_EneOverEtrueVsETtrue_EB->SetMarkerColor(4);
  P_EneOverEtrueVsETtrue_EB->GetXaxis()->SetTitle("ET true");
  P_EneOverEtrueVsETtrue_EB->GetYaxis()->SetTitle("Ene / Etrue");
  P_EneOverEtrueVsETtrue_EE->SetTitle("");
  P_EneOverEtrueVsETtrue_EE->SetMarkerStyle(20);
  P_EneOverEtrueVsETtrue_EE->SetMarkerSize(0.8);
  P_EneOverEtrueVsETtrue_EE->SetMarkerColor(4);
  P_EneOverEtrueVsETtrue_EE->GetXaxis()->SetTitle("ET true");
  P_EneOverEtrueVsETtrue_EE->GetYaxis()->SetTitle("Ene / Etrue");




  // booking histos - H/E studies
  TH1F *H_trueET_EB        = new TH1F("H_trueET_EB","H_trueET_EB",25, 0.,4500.);
  TH1F *H_trueET_EE        = new TH1F("H_trueET_EE","H_trueET_EE",25, 0.,4500.);
  TH1F *H_trueET_okHoE_EB  = new TH1F("H_trueET_okHoE_EB", "H_trueET_okHoE_EB", 25, 0.,4500.);
  TH1F *H_trueET_okHoE_EE  = new TH1F("H_trueET_okHoE_EE", "H_trueET_okHoE_EE", 25, 0.,4500.);
  TH1F *H_trueET_okHoET_EB = new TH1F("H_trueET_okHoET_EB","H_trueET_okHoET_EB",25, 0.,4500.);
  TH1F *H_trueET_okHoET_EE = new TH1F("H_trueET_okHoET_EE","H_trueET_okHoET_EE",25, 0.,4500.);
  H_trueET_EB->Sumw2();
  H_trueET_EE->Sumw2();
  H_trueET_okHoE_EB->Sumw2();
  H_trueET_okHoE_EE->Sumw2();
  H_trueET_okHoET_EB->Sumw2();
  H_trueET_okHoET_EE->Sumw2();

  TH1F *H_scEta            = new TH1F("H_scEta",           "H_scEta",           57,-2.5,2.5);
  TH1F *H_scEta_okHoE      = new TH1F("H_scEta_okHoE",     "H_scEta_okHoE",     57,-2.5,2.5);
  TH1F *H_scEta_okHoET     = new TH1F("H_scEta_okHoET",    "H_scEta_okHoET",    57,-2.5,2.5);
  TH1F *H_scEta_okHoEtrue  = new TH1F("H_scEta_okHoEtrue", "H_scEta_okHoEtrue", 57,-2.5,2.5);
  TH1F *H_scEta_okHoEtrueT = new TH1F("H_scEta_okHoEtrueT","H_scEta_okHoEtrueT",57,-2.5,2.5);
  H_scEta->Sumw2();
  H_scEta_okHoE->Sumw2();
  H_scEta_okHoET->Sumw2();
  H_scEta_okHoEtrue->Sumw2();
  H_scEta_okHoEtrueT->Sumw2();

  TH1F *H_scPhi            = new TH1F("H_scPhi",           "H_scPhi",            72,-3.14,3.14);
  TH1F *H_scPhi_okHoE      = new TH1F("H_scPhi_okHoE",     "H_scPhi_okHoE",      72,-3.14,3.14);
  TH1F *H_scPhi_okHoET     = new TH1F("H_scPhi_okHoET",    "H_scPhi_okHoET",     72,-3.14,3.14);
  TH1F *H_scPhi_okHoEtrue  = new TH1F("H_scPhi_okHoEtrue", "H_scPhi_okHoEtrue",  72,-3.14,3.14);
  TH1F *H_scPhi_okHoEtrueT = new TH1F("H_scPhi_okHoEtrueT","H_scPhi_okHoEtrueT", 72,-3.14,3.14);
  H_scPhi->Sumw2();
  H_scPhi_okHoE->Sumw2();
  H_scPhi_okHoET->Sumw2();
  H_scPhi_okHoEtrue->Sumw2();
  H_scPhi_okHoEtrueT->Sumw2();

  TProfile *P_HvsScEta = new TProfile("P_HvsScEta","P_HvsScEta",57, -2.5,2.5);
  P_HvsScEta->SetTitle("");
  P_HvsScEta->SetMarkerStyle(20);
  P_HvsScEta->SetMarkerSize(0.8);

  TProfile *P_HvsScPhi = new TProfile("P_HvsScPhi","P_HvsScPhi",72,-3.14,3.14);
  P_HvsScPhi->SetTitle("");
  P_HvsScPhi->SetMarkerStyle(20);
  P_HvsScPhi->SetMarkerSize(0.8);

  TH1F *H_EneOverEtrue_okHoE_EB = new TH1F("H_EneOverEtrue_okHoE_EB","H_EneOverEtrue_okHoE_EB",50, 0.,1.2);
  TH1F *H_EneOverEtrue_okHoE_EE = new TH1F("H_EneOverEtrue_okHoE_EE","H_EneOverEtrue_okHoE_EE",50, 0.,1.2);
  TH1F *H_EneOverEtrue_noHoE_EB = new TH1F("H_EneOverEtrue_noHoE_EB","H_EneOverEtrue_noHoE_EB",50, 0.,1.2);
  TH1F *H_EneOverEtrue_noHoE_EE = new TH1F("H_EneOverEtrue_noHoE_EE","H_EneOverEtrue_noHoE_EE",50, 0.,1.2);
  H_EneOverEtrue_okHoE_EB->Sumw2();
  H_EneOverEtrue_okHoE_EE->Sumw2();
  H_EneOverEtrue_okHoE_EB->SetTitle("");
  H_EneOverEtrue_okHoE_EE->SetTitle("");
  H_EneOverEtrue_okHoE_EB->GetXaxis()->SetTitle("Ene / Etrue");
  H_EneOverEtrue_okHoE_EE->GetXaxis()->SetTitle("Ene / Etrue");

  TH1F *H_EneOverEtrue_okHoET_EB = new TH1F("H_EneOverEtrue_okHoET_EB","H_EneOverEtrue_okHoET_EB",50, 0.,1.2);
  TH1F *H_EneOverEtrue_okHoET_EE = new TH1F("H_EneOverEtrue_okHoET_EE","H_EneOverEtrue_okHoET_EE",50, 0.,1.2);
  TH1F *H_EneOverEtrue_noHoET_EB = new TH1F("H_EneOverEtrue_noHoET_EB","H_EneOverEtrue_noHoET_EB",50, 0.,1.2);
  TH1F *H_EneOverEtrue_noHoET_EE = new TH1F("H_EneOverEtrue_noHoET_EE","H_EneOverEtrue_noHoET_EE",50, 0.,1.2);
  H_EneOverEtrue_okHoET_EB->Sumw2();
  H_EneOverEtrue_okHoET_EE->Sumw2();
  H_EneOverEtrue_okHoET_EB->SetTitle("");
  H_EneOverEtrue_okHoET_EE->SetTitle("");
  H_EneOverEtrue_okHoET_EB->GetXaxis()->SetTitle("Ene / Etrue");
  H_EneOverEtrue_okHoET_EE->GetXaxis()->SetTitle("Ene / Etrue");

  // same studies, but now removing gaps
  TH1F *HnogaptrueET_EB        = new TH1F("HnogaptrueET_EB","HnogaptrueET_EB",25, 0.,4500.);
  TH1F *HnogaptrueET_EE        = new TH1F("HnogaptrueET_EE","HnogaptrueET_EE",25, 0.,4500.);
  TH1F *HnogaptrueET_okHoE_EB  = new TH1F("HnogaptrueET_okHoE_EB", "HnogaptrueET_okHoE_EB", 25, 0.,4500.);
  TH1F *HnogaptrueET_okHoE_EE  = new TH1F("HnogaptrueET_okHoE_EE", "HnogaptrueET_okHoE_EE", 25, 0.,4500.);
  TH1F *HnogaptrueET_okHoET_EB = new TH1F("HnogaptrueET_okHoET_EB","HnogaptrueET_okHoET_EB",25, 0.,4500.);
  TH1F *HnogaptrueET_okHoET_EE = new TH1F("HnogaptrueET_okHoET_EE","HnogaptrueET_okHoET_EE",25, 0.,4500.);
  HnogaptrueET_EB->Sumw2();
  HnogaptrueET_EE->Sumw2();
  HnogaptrueET_okHoE_EB->Sumw2();
  HnogaptrueET_okHoE_EE->Sumw2();
  HnogaptrueET_okHoET_EB->Sumw2();
  HnogaptrueET_okHoET_EE->Sumw2();

  TH1F *Hnogap_scEta            = new TH1F("Hnogap_scEta",           "Hnogap_scEta",           57,-2.5,2.5);
  TH1F *Hnogap_scEta_okHoE      = new TH1F("Hnogap_scEta_okHoE",     "Hnogap_scEta_okHoE",     57,-2.5,2.5);
  TH1F *Hnogap_scEta_okHoET     = new TH1F("Hnogap_scEta_okHoET",    "Hnogap_scEta_okHoET",    57,-2.5,2.5);
  TH1F *Hnogap_scEta_okHoEtrue  = new TH1F("Hnogap_scEta_okHoEtrue", "Hnogap_scEta_okHoEtrue", 57,-2.5,2.5);
  TH1F *Hnogap_scEta_okHoEtrueT = new TH1F("Hnogap_scEta_okHoEtrueT","Hnogap_scEta_okHoEtrueT",57,-2.5,2.5);
  Hnogap_scEta->Sumw2();
  Hnogap_scEta_okHoE->Sumw2();
  Hnogap_scEta_okHoET->Sumw2();
  Hnogap_scEta_okHoEtrue->Sumw2();
  Hnogap_scEta_okHoEtrueT->Sumw2();

  TH1F *Hnogap_scPhi            = new TH1F("Hnogap_scPhi",           "Hnogap_scPhi",            72,-3.14,3.14);
  TH1F *Hnogap_scPhi_okHoE      = new TH1F("Hnogap_scPhi_okHoE",     "Hnogap_scPhi_okHoE",      72,-3.14,3.14);
  TH1F *Hnogap_scPhi_okHoET     = new TH1F("Hnogap_scPhi_okHoET",    "Hnogap_scPhi_okHoET",     72,-3.14,3.14);
  TH1F *Hnogap_scPhi_okHoEtrue  = new TH1F("Hnogap_scPhi_okHoEtrue", "Hnogap_scPhi_okHoEtrue",  72,-3.14,3.14);
  TH1F *Hnogap_scPhi_okHoEtrueT = new TH1F("Hnogap_scPhi_okHoEtrueT","Hnogap_scPhi_okHoEtrueT", 72,-3.14,3.14);
  Hnogap_scPhi->Sumw2();
  Hnogap_scPhi_okHoE->Sumw2();
  Hnogap_scPhi_okHoET->Sumw2();
  Hnogap_scPhi_okHoEtrue->Sumw2();
  Hnogap_scPhi_okHoEtrueT->Sumw2();

  TH1F *Hnogap_EneOverEtrue_okHoE_EB = new TH1F("Hnogap_EneOverEtrue_okHoE_EB","Hnogap_EneOverEtrue_okHoE_EB",50, 0.,1.2);
  TH1F *Hnogap_EneOverEtrue_okHoE_EE = new TH1F("Hnogap_EneOverEtrue_okHoE_EE","Hnogap_EneOverEtrue_okHoE_EE",50, 0.,1.2);
  TH1F *Hnogap_EneOverEtrue_noHoE_EB = new TH1F("Hnogap_EneOverEtrue_noHoE_EB","Hnogap_EneOverEtrue_noHoE_EB",50, 0.,1.2);
  TH1F *Hnogap_EneOverEtrue_noHoE_EE = new TH1F("Hnogap_EneOverEtrue_noHoE_EE","Hnogap_EneOverEtrue_noHoE_EE",50, 0.,1.2);
  Hnogap_EneOverEtrue_okHoE_EB->Sumw2();
  Hnogap_EneOverEtrue_okHoE_EE->Sumw2();
  Hnogap_EneOverEtrue_okHoE_EB->SetTitle("");
  Hnogap_EneOverEtrue_okHoE_EE->SetTitle("");
  Hnogap_EneOverEtrue_okHoE_EB->GetXaxis()->SetTitle("Ene / Etrue");
  Hnogap_EneOverEtrue_okHoE_EE->GetXaxis()->SetTitle("Ene / Etrue");

  TH1F *Hnogap_EneOverEtrue_okHoET_EB = new TH1F("Hnogap_EneOverEtrue_okHoET_EB","Hnogap_EneOverEtrue_okHoET_EB",50, 0.,1.2);
  TH1F *Hnogap_EneOverEtrue_okHoET_EE = new TH1F("Hnogap_EneOverEtrue_okHoET_EE","Hnogap_EneOverEtrue_okHoET_EE",50, 0.,1.2);
  TH1F *Hnogap_EneOverEtrue_noHoET_EB = new TH1F("Hnogap_EneOverEtrue_noHoET_EB","Hnogap_EneOverEtrue_noHoET_EB",50, 0.,1.2);
  TH1F *Hnogap_EneOverEtrue_noHoET_EE = new TH1F("Hnogap_EneOverEtrue_noHoET_EE","Hnogap_EneOverEtrue_noHoET_EE",50, 0.,1.2);
  Hnogap_EneOverEtrue_okHoET_EB->Sumw2();
  Hnogap_EneOverEtrue_okHoET_EE->Sumw2();
  Hnogap_EneOverEtrue_okHoET_EB->SetTitle("");
  Hnogap_EneOverEtrue_okHoET_EE->SetTitle("");
  Hnogap_EneOverEtrue_okHoET_EB->GetXaxis()->SetTitle("Ene / Etrue");
  Hnogap_EneOverEtrue_okHoET_EE->GetXaxis()->SetTitle("Ene / Etrue");


  // booking histos - 2012ID studies
  TH1F *H_trueET_okSee_EB  = new TH1F("H_trueET_okSee_EB", "H_trueET_okSee_EB", 25, 0.,4500.);
  TH1F *H_trueET_okSee_EE  = new TH1F("H_trueET_okSee_EE", "H_trueET_okSee_EE", 25, 0.,4500.);
  H_trueET_okSee_EB->Sumw2();
  H_trueET_okSee_EE->Sumw2();

  TH1F *H_trueET_okChIso_EB  = new TH1F("H_trueET_okChIso_EB", "H_trueET_okChIso_EB", 25, 0.,4500.);
  TH1F *H_trueET_okChIso_EE  = new TH1F("H_trueET_okChIso_EE", "H_trueET_okChIso_EE", 25, 0.,4500.);
  H_trueET_okChIso_EB->Sumw2();
  H_trueET_okChIso_EE->Sumw2();

  TH1F *H_trueET_okNhIso_EB  = new TH1F("H_trueET_okNhIso_EB", "H_trueET_okNhIso_EB", 25, 0.,4500.);
  TH1F *H_trueET_okNhIso_EE  = new TH1F("H_trueET_okNhIso_EE", "H_trueET_okNhIso_EE", 25, 0.,4500.);
  H_trueET_okNhIso_EB->Sumw2();
  H_trueET_okNhIso_EE->Sumw2();


  // loop over entries
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;


    // restricting the analysis to photons matching MC truth
    if (mctruth_trueEta<-800) continue;

    // true photon must be in the ECAL fiducial region
    if (fabs(mctruth_trueEta)>2.5) continue;
    if (fabs(mctruth_trueEta)>1.4442 && fabs(mctruth_trueEta)<1.566) continue;

    // also the supercluster of the reco photon must be in the ECAL fiducial region
    if (fabs(supercluster_scEta)>2.5) continue;
    if (fabs(supercluster_scEta)>1.4442 &&fabs(supercluster_scEta)<1.566) continue;

    // if (identification_r9<0.9) continue;


    // extended eta gaps for H/E studies
    // default is 1 xtal only, we remove 2.5 per side
    bool inExtEBEtaGap = false;
    int ieta = abs(kinematics_iEta);
    if ( ieta==0 ||  ieta==1 ||  ieta==2 ||  
	 ieta==23 || ieta==24 || ieta==25 || ieta==26 || ieta==27 || 
	 ieta==43 || ieta==44 || ieta==45 || ieta==46 || ieta==47 ||
	 ieta==63 || ieta==64 || ieta==65 || ieta==66 || ieta==67 ||
	 ieta==83 || ieta==84 || ieta==85) inExtEBEtaGap = true;

    // default is 1 xtal per side only, we remove 2 per side    
    bool inExtEBPhiGap = false;
    /*
    int iphi = abs(kinematics_iPhi);
    if (iphi==0 || iphi==1 ||
	iphi== 18 || iphi== 19 ||  iphi==20 ||  iphi==21 || iphi==22  ||  iphi==23  ||  
	iphi== 38 || iphi== 39 ||  iphi==40 ||  iphi==41 || iphi==42  ||  iphi==43  ||  
	iphi== 58 || iphi== 59 ||  iphi==60 ||  iphi==61 || iphi==62  ||  iphi==63  ||  
	iphi== 78 || iphi== 79 ||  iphi==80 ||  iphi==81 || iphi==82  ||  iphi==83  ||  
	iphi== 98 || iphi== 99 || iphi==100 || iphi==101 || iphi==102 ||  iphi==103 ||  
	iphi==118 || iphi==119 || iphi==120 || iphi==121 || iphi==122 ||  iphi==123 ||  
	iphi==138 || iphi==139 || iphi==140 || iphi==141 || iphi==142 ||  iphi==143 ||  
	iphi==158 || iphi==159 || iphi==160 || iphi==161 || iphi==162 ||  iphi==163 ||  
	iphi==178 || iphi==179 || iphi==180 || iphi==181 || iphi==182 ||  iphi==183 ||  
	iphi==198 || iphi==199 || iphi==200 || iphi==201 || iphi==202 ||  iphi==203 ||  
	iphi==218 || iphi==219 || iphi==220 || iphi==221 || iphi==222 ||  iphi==223 ||  
	iphi==238 || iphi==239 || iphi==240 || iphi==241 || iphi==242 ||  iphi==243 ||  
	iphi==258 || iphi==259 || iphi==260 || iphi==261 || iphi==262 ||  iphi==263 ||  
	iphi==278 || iphi==279 || iphi==280 || iphi==281 || iphi==282 ||  iphi==283 ||  
	iphi==298 || iphi==299 || iphi==300 || iphi==301 || iphi==302 ||  iphi==303 ||  
	iphi==318 || iphi==319 || iphi==320 || iphi==321 || iphi==322 ||  iphi==323 ||  
	iphi==338 || iphi==339 || iphi==340 || iphi==341 || iphi==342 ||  iphi==343 ||  
	iphi==358 || iphi==359 || iphi==360 ) inExtEBPhiGap = true;
    */

    // filling histos

    // full eta range
    H_eta -> Fill(kinematics_eta);

    if (mctruth_truePt>200) { 

      H_scEta -> Fill(supercluster_scEta);
      H_scPhi -> Fill(supercluster_scPhi);

      if (identification_hoe<0.05)  H_scEta_okHoE  -> Fill(supercluster_scEta);
      if (identification_htoe<0.05) H_scEta_okHoET -> Fill(supercluster_scEta);
      if (identification_hoe<0.05)  H_scPhi_okHoE  -> Fill(supercluster_scPhi);
      if (identification_htoe<0.05) H_scPhi_okHoET -> Fill(supercluster_scPhi);

      float hoetrue  = identification_hoe*(energy_energy/mctruth_trueEnergy);
      float htoetrue = identification_htoe*(energy_energy/mctruth_trueEnergy);
      if (hoetrue<0.05)  H_scEta_okHoEtrue  -> Fill(supercluster_scEta);
      if (htoetrue<0.05) H_scEta_okHoEtrueT -> Fill(supercluster_scEta);
      if (hoetrue<0.05)  H_scPhi_okHoEtrue  -> Fill(supercluster_scPhi);
      if (htoetrue<0.05) H_scPhi_okHoEtrueT -> Fill(supercluster_scPhi);

      float theH = identification_htoe*energy_energy;
      P_HvsScEta -> Fill(supercluster_scEta,theH);
      P_HvsScPhi -> Fill(supercluster_scPhi,theH);

      // same but removing (extended) gaps
      if (!kinematics_isEBEtaGap && !kinematics_isEBPhiGap && 
	  !kinematics_isEERingGap && !kinematics_isEEDeeGap && 
	  !kinematics_isEBEEGap && !inExtEBEtaGap && !inExtEBPhiGap ) {
	
	Hnogap_scEta -> Fill(supercluster_scEta);
	Hnogap_scPhi -> Fill(supercluster_scPhi);

	if (identification_hoe<0.05)  Hnogap_scEta_okHoE  -> Fill(supercluster_scEta);
	if (identification_htoe<0.05) Hnogap_scEta_okHoET -> Fill(supercluster_scEta);
	if (identification_hoe<0.05)  Hnogap_scPhi_okHoE  -> Fill(supercluster_scPhi);
	if (identification_htoe<0.05) Hnogap_scPhi_okHoET -> Fill(supercluster_scPhi);
	
	if (hoetrue<0.05)  Hnogap_scEta_okHoEtrue  -> Fill(supercluster_scEta);
	if (htoetrue<0.05) Hnogap_scEta_okHoEtrueT -> Fill(supercluster_scEta);
	if (hoetrue<0.05)  Hnogap_scPhi_okHoEtrue  -> Fill(supercluster_scPhi);
	if (htoetrue<0.05) Hnogap_scPhi_okHoEtrueT -> Fill(supercluster_scPhi);
      }
    }


    // isolation PU-corrections
    float correctedChIso = corrChIso(isolation_chHadIso, supercluster_scEta);
    float correctedNIso  = corrNIso(isolation_nHadIso, supercluster_scEta);
    if (correctedChIso<0) correctedChIso=0;
    if (correctedNIso<0)  correctedNIso=0;

    // EB / EE
    if (fabs(kinematics_eta)<1.5) {

      H_pt_EB  -> Fill(kinematics_pt);
      
      // energy studies
      H_E25OverEtrue_EB -> Fill(energy_e5x5/mctruth_trueEnergy);
      H_EneOverEtrue_EB -> Fill(energy_energy/mctruth_trueEnergy);
      P_E25OverEtrueVsETtrue_EB -> Fill(mctruth_truePt, (energy_e5x5/mctruth_trueEnergy));
      P_EneOverEtrueVsETtrue_EB -> Fill(mctruth_truePt, (energy_energy/mctruth_trueEnergy));

      // H/E cut efficiency studies
      H_trueET_EB -> Fill(mctruth_truePt);
      if (identification_hoe<0.05)  H_trueET_okHoE_EB  -> Fill(mctruth_truePt);
      if (identification_htoe<0.05) H_trueET_okHoET_EB -> Fill(mctruth_truePt);
      if (identification_hoe<0.05)  H_EneOverEtrue_okHoE_EB  -> Fill(energy_energy/mctruth_trueEnergy);
      if (identification_htoe<0.05) H_EneOverEtrue_okHoET_EB -> Fill(energy_energy/mctruth_trueEnergy);
      if (identification_hoe>0.05)  H_EneOverEtrue_noHoE_EB  -> Fill(energy_energy/mctruth_trueEnergy);
      if (identification_htoe>0.05) H_EneOverEtrue_noHoET_EB -> Fill(energy_energy/mctruth_trueEnergy);

      // H/E cut efficiency studies - no gap
      if (!kinematics_isEBEtaGap && !kinematics_isEBPhiGap && !kinematics_isEBEEGap && !inExtEBEtaGap && !inExtEBPhiGap) {	  
	HnogaptrueET_EB -> Fill(mctruth_truePt);
	if (identification_hoe<0.05)  HnogaptrueET_okHoE_EB  -> Fill(mctruth_truePt);
	if (identification_htoe<0.05) HnogaptrueET_okHoET_EB -> Fill(mctruth_truePt);
	if (identification_hoe<0.05)  Hnogap_EneOverEtrue_okHoE_EB  -> Fill(energy_energy/mctruth_trueEnergy);
	if (identification_htoe<0.05) Hnogap_EneOverEtrue_okHoET_EB -> Fill(energy_energy/mctruth_trueEnergy);
	if (identification_hoe>0.05)  Hnogap_EneOverEtrue_noHoE_EB  -> Fill(energy_energy/mctruth_trueEnergy);
	if (identification_htoe>0.05) Hnogap_EneOverEtrue_noHoET_EB -> Fill(energy_energy/mctruth_trueEnergy);
      }

      // numbers for H/E studies:
      totalEB++;
      if (identification_hoe<0.05)  okHoE_totalEB++;
      if (identification_htoe<0.05) okHoET_totalEB++;

      if (!kinematics_isEBEtaGap && !kinematics_isEBPhiGap && !kinematics_isEBEEGap && !inExtEBEtaGap && !inExtEBPhiGap) {
	totalEB_noGap++;
	if (identification_hoe<0.05)  okHoE_totalEB_noGap++;
	if (identification_htoe<0.05) okHoET_totalEB_noGap++;

	if (identification_r9>0.9){ 
	  totalEB_noGap_highR9++;
	  if (identification_hoe<0.05)  okHoE_totalEB_noGap_highR9++;
	  if (identification_htoe<0.05) okHoET_totalEB_noGap_highR9++;
	}
      }

      // 2012 loose photon ID. chiara, da capire che cono e' (ci vorrebbe 0.3)
      if (identification_sigmaIetaIeta<0.012)          H_trueET_okSee_EB   -> Fill(mctruth_truePt);
      if (correctedChIso<2.6)                          H_trueET_okChIso_EB -> Fill(mctruth_truePt);
      if (correctedNIso<(3.5 + 0.04*kinematics_pt))    H_trueET_okNhIso_EB -> Fill(mctruth_truePt);

    } else {      
      
      H_pt_EE  -> Fill(kinematics_pt);
      
      // energy studies
      H_E25OverEtrue_EE  -> Fill(energy_e5x5/mctruth_trueEnergy);
      H_EneOverEtrue_EE  -> Fill(energy_energy/mctruth_trueEnergy);
      P_E25OverEtrueVsETtrue_EE  -> Fill(mctruth_truePt, (energy_e5x5/mctruth_trueEnergy));
      P_EneOverEtrueVsETtrue_EE  -> Fill(mctruth_truePt, (energy_energy/mctruth_trueEnergy));

      // H/E cut efficiency studies    
      H_trueET_EE -> Fill(mctruth_truePt);
      if (identification_hoe<0.05)  H_trueET_okHoE_EE  -> Fill(mctruth_truePt);
      if (identification_htoe<0.05) H_trueET_okHoET_EE -> Fill(mctruth_truePt);
      if (identification_hoe<0.05)  H_EneOverEtrue_okHoE_EE  -> Fill(energy_energy/mctruth_trueEnergy);
      if (identification_htoe<0.05) H_EneOverEtrue_okHoET_EE -> Fill(energy_energy/mctruth_trueEnergy);
      if (identification_hoe>0.05)  H_EneOverEtrue_noHoE_EE  -> Fill(energy_energy/mctruth_trueEnergy);
      if (identification_htoe>0.05) H_EneOverEtrue_noHoET_EE -> Fill(energy_energy/mctruth_trueEnergy);

      // H/E cut efficiency studies - no gap   
      if (!kinematics_isEERingGap && !kinematics_isEEDeeGap && !kinematics_isEBEEGap) {
	HnogaptrueET_EE -> Fill(mctruth_truePt);
	if (identification_hoe<0.05)  HnogaptrueET_okHoE_EE  -> Fill(mctruth_truePt);
	if (identification_htoe<0.05) HnogaptrueET_okHoET_EE -> Fill(mctruth_truePt);
	if (identification_hoe<0.05)  Hnogap_EneOverEtrue_okHoE_EE  -> Fill(energy_energy/mctruth_trueEnergy);
	if (identification_htoe<0.05) Hnogap_EneOverEtrue_okHoET_EE -> Fill(energy_energy/mctruth_trueEnergy);
	if (identification_hoe>0.05)  Hnogap_EneOverEtrue_noHoE_EE  -> Fill(energy_energy/mctruth_trueEnergy);
	if (identification_htoe>0.05) Hnogap_EneOverEtrue_noHoET_EE -> Fill(energy_energy/mctruth_trueEnergy);
      }

      // numbers for H/E efficiencies
      totalEE++;
      if (identification_hoe<0.05)  okHoE_totalEE++;
      if (identification_htoe<0.05) okHoET_totalEE++;

      if (!kinematics_isEERingGap && !kinematics_isEEDeeGap && !kinematics_isEBEEGap) {
	totalEE_noGap++;
	if (identification_hoe<0.05)  okHoE_totalEE_noGap++;
	if (identification_htoe<0.05) okHoET_totalEE_noGap++;

	if (identification_r9>0.9){ 
	  totalEE_noGap_highR9++;
	  if (identification_hoe<0.05)  okHoE_totalEE_noGap_highR9++;
	  if (identification_htoe<0.05) okHoET_totalEE_noGap_highR9++;
	}
      }

      // 2012 loose photon ID. chiara, da capire che cono e' (ci vorrebbe 0.3)
      if (identification_sigmaIetaIeta<0.034)          H_trueET_okSee_EE   -> Fill(mctruth_truePt);
      if (correctedChIso<2.3)                          H_trueET_okChIso_EE -> Fill(mctruth_truePt);
      if (correctedNIso<(2.9 + 0.04*kinematics_pt))    H_trueET_okNhIso_EE -> Fill(mctruth_truePt);
    }
    
  }  // loop over entries


  // H/E efficiency numbers 
  cout << endl;
  cout << endl;
  cout << "EB: " << endl;
  cout << "no cut, H/E:" << okHoE_totalEB/totalEB             << ", H/E with tower: " << okHoET_totalEB/totalEB             << endl;
  cout << "no gap, H/E:" << okHoE_totalEB_noGap/totalEB_noGap << ", H/E with tower: " << okHoET_totalEB_noGap/totalEB_noGap << endl;
  cout << "no gap & high R9, H/E:" << okHoE_totalEB_noGap_highR9/totalEB_noGap_highR9 
       << ", H/E with tower: "     << okHoET_totalEB_noGap_highR9/totalEB_noGap_highR9 << endl;

  cout << "EE: " << endl;
  cout << "no cut, H/E:" << okHoE_totalEE/totalEE             << ", H/E with tower: " << okHoET_totalEE/totalEE             << endl;
  cout << "no gap, H/E:" << okHoE_totalEE_noGap/totalEE_noGap << ", H/E with tower: " << okHoET_totalEE_noGap/totalEE_noGap << endl;
  cout << "no gap & high R9, H/E:" << okHoE_totalEE_noGap_highR9/totalEE_noGap_highR9 
       << ", H/E with tower: "     << okHoET_totalEE_noGap_highR9/totalEE_noGap_highR9 << endl;
  cout << endl;
  cout << endl;


  // Efficiencies - H/E studies
  TGraphAsymmErrors *h1_effHoE_EB = new TGraphAsymmErrors(H_trueET_okHoE_EB,H_trueET_EB);
  TGraphAsymmErrors *h1_effHoE_EE = new TGraphAsymmErrors(H_trueET_okHoE_EE,H_trueET_EE);
  h1_effHoE_EB->SetTitle("H/E cut efficiency");
  h1_effHoE_EE->SetTitle("H/E cut efficiency");
  h1_effHoE_EB->GetXaxis()->SetTitle("ET true");
  h1_effHoE_EE->GetXaxis()->SetTitle("ET true");
  h1_effHoE_EB->SetMarkerStyle(20);
  h1_effHoE_EE->SetMarkerStyle(20);
  h1_effHoE_EB->SetMarkerColor(4);
  h1_effHoE_EE->SetMarkerColor(4);
  // 
  TGraphAsymmErrors *h1_effHoET_EB = new TGraphAsymmErrors(H_trueET_okHoET_EB,H_trueET_EB);
  TGraphAsymmErrors *h1_effHoET_EE = new TGraphAsymmErrors(H_trueET_okHoET_EE,H_trueET_EE);
  h1_effHoET_EB->SetTitle("H/E cut efficiency");
  h1_effHoET_EE->SetTitle("H/E cut efficiency");
  h1_effHoET_EB->GetXaxis()->SetTitle("ET true");
  h1_effHoET_EE->GetXaxis()->SetTitle("ET true");
  h1_effHoET_EB->SetMarkerStyle(20);
  h1_effHoET_EE->SetMarkerStyle(20);
  h1_effHoET_EB->SetMarkerColor(2);
  h1_effHoET_EE->SetMarkerColor(2);
  // 
  TGraphAsymmErrors *h1_effHoEvsScEta = new TGraphAsymmErrors(H_scEta_okHoE,H_scEta);
  h1_effHoEvsScEta->SetTitle("H/E cut efficiency");
  h1_effHoEvsScEta->SetTitle("SC #eta");
  h1_effHoEvsScEta->SetMarkerStyle(20);
  h1_effHoEvsScEta->SetMarkerColor(4);
  //
  TGraphAsymmErrors *h1_effHoETvsScEta = new TGraphAsymmErrors(H_scEta_okHoET,H_scEta);
  h1_effHoETvsScEta->SetTitle("H/E cut efficiency");
  h1_effHoETvsScEta->SetTitle("SC #eta");
  h1_effHoETvsScEta->SetMarkerStyle(20);
  h1_effHoETvsScEta->SetMarkerColor(2);
  //
  TGraphAsymmErrors *h1_effHoEvsScPhi = new TGraphAsymmErrors(H_scPhi_okHoE,H_scPhi);
  h1_effHoEvsScPhi->SetTitle("H/E cut efficiency");
  h1_effHoEvsScPhi->SetTitle("SC #phi");
  h1_effHoEvsScPhi->SetMarkerStyle(20);
  h1_effHoEvsScPhi->SetMarkerColor(4);
  //
  TGraphAsymmErrors *h1_effHoETvsScPhi = new TGraphAsymmErrors(H_scPhi_okHoET,H_scPhi);
  h1_effHoETvsScPhi->SetTitle("H/E cut efficiency");
  h1_effHoETvsScPhi->SetTitle("SC #phi");
  h1_effHoETvsScPhi->SetMarkerStyle(20);
  h1_effHoETvsScPhi->SetMarkerColor(2);
  //
  //
  TGraphAsymmErrors *h1_effHoEtruevsScEta = new TGraphAsymmErrors(H_scEta_okHoEtrue,H_scEta);
  h1_effHoEtruevsScEta->SetTitle("H/Etrue cut efficiency");
  h1_effHoEtruevsScEta->SetTitle("SC #eta");
  h1_effHoEtruevsScEta->SetMarkerStyle(20);
  h1_effHoEtruevsScEta->SetMarkerColor(4);
  //
  TGraphAsymmErrors *h1_effHoEtrueTvsScEta = new TGraphAsymmErrors(H_scEta_okHoEtrueT,H_scEta);
  h1_effHoEtrueTvsScEta->SetTitle("H/Etrue cut efficiency");
  h1_effHoEtrueTvsScEta->SetTitle("SC #eta");
  h1_effHoEtrueTvsScEta->SetMarkerStyle(20);
  h1_effHoEtrueTvsScEta->SetMarkerColor(2);
  //
  //
  TGraphAsymmErrors *h1_effHoEtruevsScPhi = new TGraphAsymmErrors(H_scPhi_okHoEtrue,H_scPhi);
  h1_effHoEtruevsScPhi->SetTitle("H/Etrue cut efficiency");
  h1_effHoEtruevsScPhi->SetTitle("SC #phi");
  h1_effHoEtruevsScPhi->SetMarkerStyle(20);
  h1_effHoEtruevsScPhi->SetMarkerColor(4);
  //
  TGraphAsymmErrors *h1_effHoEtrueTvsScPhi = new TGraphAsymmErrors(H_scPhi_okHoEtrueT,H_scPhi);
  h1_effHoEtrueTvsScPhi->SetTitle("H/Etrue cut efficiency");
  h1_effHoEtrueTvsScPhi->SetTitle("SC #phi");
  h1_effHoEtrueTvsScPhi->SetMarkerStyle(20);
  h1_effHoEtrueTvsScPhi->SetMarkerColor(2);

  // Efficiencies - H/E studies out of gaps
  TGraphAsymmErrors *h1nogap_effHoE_EB = new TGraphAsymmErrors(HnogaptrueET_okHoE_EB,HnogaptrueET_EB);
  TGraphAsymmErrors *h1nogap_effHoE_EE = new TGraphAsymmErrors(HnogaptrueET_okHoE_EE,HnogaptrueET_EE);
  h1nogap_effHoE_EB->SetTitle("H/E cut efficiency");
  h1nogap_effHoE_EE->SetTitle("H/E cut efficiency");
  h1nogap_effHoE_EB->GetXaxis()->SetTitle("ET true");
  h1nogap_effHoE_EE->GetXaxis()->SetTitle("ET true");
  h1nogap_effHoE_EB->SetMarkerStyle(20);
  h1nogap_effHoE_EE->SetMarkerStyle(20);
  h1nogap_effHoE_EB->SetMarkerColor(4);
  h1nogap_effHoE_EE->SetMarkerColor(4);
  // 
  TGraphAsymmErrors *h1nogap_effHoET_EB = new TGraphAsymmErrors(HnogaptrueET_okHoET_EB,HnogaptrueET_EB);
  TGraphAsymmErrors *h1nogap_effHoET_EE = new TGraphAsymmErrors(HnogaptrueET_okHoET_EE,HnogaptrueET_EE);
  h1nogap_effHoET_EB->SetTitle("H/E cut efficiency");
  h1nogap_effHoET_EE->SetTitle("H/E cut efficiency");
  h1nogap_effHoET_EB->GetXaxis()->SetTitle("ET true");
  h1nogap_effHoET_EE->GetXaxis()->SetTitle("ET true");
  h1nogap_effHoET_EB->SetMarkerStyle(20);
  h1nogap_effHoET_EE->SetMarkerStyle(20);
  h1nogap_effHoET_EB->SetMarkerColor(2);
  h1nogap_effHoET_EE->SetMarkerColor(2);
  // 
  TGraphAsymmErrors *h1nogap_effHoEvsScEta = new TGraphAsymmErrors(Hnogap_scEta_okHoE,Hnogap_scEta);
  h1nogap_effHoEvsScEta->SetTitle("H/E cut efficiency");
  h1nogap_effHoEvsScEta->SetTitle("SC #eta");
  h1nogap_effHoEvsScEta->SetMarkerStyle(20);
  h1nogap_effHoEvsScEta->SetMarkerColor(4);
  //
  TGraphAsymmErrors *h1nogap_effHoETvsScEta = new TGraphAsymmErrors(Hnogap_scEta_okHoET,Hnogap_scEta);
  h1nogap_effHoETvsScEta->SetTitle("H/E cut efficiency");
  h1nogap_effHoETvsScEta->SetTitle("SC #eta");
  h1nogap_effHoETvsScEta->SetMarkerStyle(20);
  h1nogap_effHoETvsScEta->SetMarkerColor(2);
  //
  TGraphAsymmErrors *h1nogap_effHoEvsScPhi = new TGraphAsymmErrors(Hnogap_scPhi_okHoE,Hnogap_scPhi);
  h1nogap_effHoEvsScPhi->SetTitle("H/E cut efficiency");
  h1nogap_effHoEvsScPhi->SetTitle("SC #phi");
  h1nogap_effHoEvsScPhi->SetMarkerStyle(20);
  h1nogap_effHoEvsScPhi->SetMarkerColor(4);
  //
  TGraphAsymmErrors *h1nogap_effHoETvsScPhi = new TGraphAsymmErrors(Hnogap_scPhi_okHoET,Hnogap_scPhi);
  h1nogap_effHoETvsScPhi->SetTitle("H/E cut efficiency");
  h1nogap_effHoETvsScPhi->SetTitle("SC #phi");
  h1nogap_effHoETvsScPhi->SetMarkerStyle(20);
  h1nogap_effHoETvsScPhi->SetMarkerColor(2);
  //
  //
  TGraphAsymmErrors *h1nogap_effHoEtruevsScEta = new TGraphAsymmErrors(Hnogap_scEta_okHoEtrue,Hnogap_scEta);
  h1nogap_effHoEtruevsScEta->SetTitle("H/Etrue cut efficiency");
  h1nogap_effHoEtruevsScEta->SetTitle("SC #eta");
  h1nogap_effHoEtruevsScEta->SetMarkerStyle(20);
  h1nogap_effHoEtruevsScEta->SetMarkerColor(4);
  //
  TGraphAsymmErrors *h1nogap_effHoEtrueTvsScEta = new TGraphAsymmErrors(Hnogap_scEta_okHoEtrueT,Hnogap_scEta);
  h1nogap_effHoEtrueTvsScEta->SetTitle("H/Etrue cut efficiency");
  h1nogap_effHoEtrueTvsScEta->SetTitle("SC #eta");
  h1nogap_effHoEtrueTvsScEta->SetMarkerStyle(20);
  h1nogap_effHoEtrueTvsScEta->SetMarkerColor(2);
  //
  //
  TGraphAsymmErrors *h1nogap_effHoEtruevsScPhi = new TGraphAsymmErrors(Hnogap_scPhi_okHoEtrue,Hnogap_scPhi);
  h1nogap_effHoEtruevsScPhi->SetTitle("H/Etrue cut efficiency");
  h1nogap_effHoEtruevsScPhi->SetTitle("SC #phi");
  h1nogap_effHoEtruevsScPhi->SetMarkerStyle(20);
  h1nogap_effHoEtruevsScPhi->SetMarkerColor(4);
  //
  TGraphAsymmErrors *h1nogap_effHoEtrueTvsScPhi = new TGraphAsymmErrors(Hnogap_scPhi_okHoEtrueT,Hnogap_scPhi);
  h1nogap_effHoEtrueTvsScPhi->SetTitle("H/Etrue cut efficiency");
  h1nogap_effHoEtrueTvsScPhi->SetTitle("SC #phi");
  h1nogap_effHoEtrueTvsScPhi->SetMarkerStyle(20);
  h1nogap_effHoEtrueTvsScPhi->SetMarkerColor(2);

  // Efficiency studies, 2012 selection
  TGraphAsymmErrors *h1_effSee_EB = new TGraphAsymmErrors(H_trueET_okSee_EB,H_trueET_EB);
  TGraphAsymmErrors *h1_effSee_EE = new TGraphAsymmErrors(H_trueET_okSee_EE,H_trueET_EE);
  h1_effSee_EB->SetTitle("SigmaIeIe cut efficiency");
  h1_effSee_EE->SetTitle("SigmaIeIe cut efficiency");
  h1_effSee_EB->GetXaxis()->SetTitle("ET true");
  h1_effSee_EE->GetXaxis()->SetTitle("ET true");
  h1_effSee_EB->SetMarkerStyle(20);
  h1_effSee_EE->SetMarkerStyle(20);
  h1_effSee_EB->SetMarkerColor(4);
  h1_effSee_EE->SetMarkerColor(2);
  // 
  TGraphAsymmErrors *h1_effChIso_EB = new TGraphAsymmErrors(H_trueET_okChIso_EB,H_trueET_EB);
  TGraphAsymmErrors *h1_effChIso_EE = new TGraphAsymmErrors(H_trueET_okChIso_EE,H_trueET_EE);
  h1_effChIso_EB->SetTitle("chargedH isolation cut efficiency");
  h1_effChIso_EE->SetTitle("chargedH isolation cut efficiency");
  h1_effChIso_EB->GetXaxis()->SetTitle("ET true");
  h1_effChIso_EE->GetXaxis()->SetTitle("ET true");
  h1_effChIso_EB->SetMarkerStyle(20);
  h1_effChIso_EE->SetMarkerStyle(20);
  h1_effChIso_EB->SetMarkerColor(4);
  h1_effChIso_EE->SetMarkerColor(2);
  // 
  TGraphAsymmErrors *h1_effNhIso_EB = new TGraphAsymmErrors(H_trueET_okNhIso_EB,H_trueET_EB);
  TGraphAsymmErrors *h1_effNhIso_EE = new TGraphAsymmErrors(H_trueET_okNhIso_EE,H_trueET_EE);
  h1_effNhIso_EB->SetTitle("neutralH isolation cut efficiency");
  h1_effNhIso_EE->SetTitle("neutralH isolation cut efficiency");
  h1_effNhIso_EB->GetXaxis()->SetTitle("ET true");
  h1_effNhIso_EE->GetXaxis()->SetTitle("ET true");
  h1_effNhIso_EB->SetMarkerStyle(20);
  h1_effNhIso_EE->SetMarkerStyle(20);
  h1_effNhIso_EB->SetMarkerColor(4);
  h1_effNhIso_EE->SetMarkerColor(2);


  // Plots
  gStyle->SetOptStat(0);

  TCanvas c1a("c1a","c1a",1);
  H_pt_EB -> Draw();
  c1a.SaveAs("pT_EB.png");
  TCanvas c1b("c1b","c1b",1);
  H_pt_EE -> Draw();
  c1b.SaveAs("pT_EE.png");

  TCanvas c2("c2","c2",1);
  H_eta -> Draw();
  c2.SaveAs("eta.png");

  TCanvas c4a("c4a","c4a",1);
  H_E25OverEtrue_EB -> Draw("hist");
  c4a.SaveAs("E25OverEtrue_EB.png");
  TCanvas c4b("c4b","c4b",1);
  H_E25OverEtrue_EE -> Draw("hist");
  c4b.SaveAs("E25OverEtrue_EE.png");

  TCanvas c5a("c5a","c5a",1);
  H_EneOverEtrue_EB -> Draw("hist");
  c5a.SaveAs("EneOverEtrue_EB.png");
  TCanvas c5b("c5b","c5b",1);
  H_EneOverEtrue_EE -> Draw("hist");
  c5b.SaveAs("EneOverEtrue_EE.png");


  // cosmetics
  H_E25OverEtrue_EB->SetLineColor(2);
  H_EneOverEtrue_EB->SetLineColor(4);
  H_E25OverEtrue_EE->SetLineColor(2);
  H_EneOverEtrue_EE->SetLineColor(4);
  H_E25OverEtrue_EB->SetLineWidth(2);
  H_EneOverEtrue_EB->SetLineWidth(2);
  H_E25OverEtrue_EE->SetLineWidth(2);
  H_EneOverEtrue_EE->SetLineWidth(2);

  TLegend *leg;
  leg = new TLegend(0.2,0.2,0.45,0.45);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  leg->AddEntry(P_E25OverEtrueVsETtrue_EB,  "E25/Etrue", "p");
  leg->AddEntry(P_EneOverEtrueVsETtrue_EB,  "Ene/Etrue", "p");

  TCanvas c6a("c6a","c6a",1);  
  H_EneOverEtrue_EB  -> Draw("hist");
  H_E25OverEtrue_EB  -> Draw("samehist");
  leg->Draw();
  c6a.SaveAs("ExOverEtrue_EB.png");
  TCanvas c6b("c6b","c6b",1);  
  H_EneOverEtrue_EE  -> Draw("hist");
  H_E25OverEtrue_EE  -> Draw("samehist");
  leg->Draw();
  c6b.SaveAs("ExOverEtrue_EE.png");

  /*
  TCanvas c8a("c8a","c8a",1);
  P_E25OverEtrueVsETtrue_EB -> Draw();
  c8a.SaveAs("E25OverEtrueVsETtrue_EB.png");
  TCanvas c8b("c8b","c8b",1);
  P_E25OverEtrueVsETtrue_EE -> Draw();
  c8b.SaveAs("E25OverEtrueVsETtrue_EE.png");

  TCanvas c9a("c9a","c9a",1);
  P_EneOverEtrueVsETtrue_EB -> Draw();
  c9a.SaveAs("EneOverEtrueVsETtrue_EB.png");
  TCanvas c9b("c9b","c9b",1);
  P_EneOverEtrueVsETtrue_EE -> Draw();
  c9b.SaveAs("EneOverEtrueVsETtrue_EE.png");
  */

  // cosmetics
  P_E25OverEtrueVsETtrue_EB->SetLineColor(2);
  P_EneOverEtrueVsETtrue_EB->SetLineColor(4);
  P_E25OverEtrueVsETtrue_EE->SetLineColor(2);
  P_EneOverEtrueVsETtrue_EE->SetLineColor(4);
  
  TH2F *myH = new TH2F("myH","",100,0,4500,100,0.8,1.1);
  myH->GetXaxis()->SetTitle("ET true");
  myH->GetYaxis()->SetTitle("Ex / Etrue");

  TCanvas c10a("c10a","c10a",1);  
  myH->Draw();
  P_E25OverEtrueVsETtrue_EB  -> Draw("sameP");
  P_EneOverEtrueVsETtrue_EB  -> Draw("sameP");
  leg->Draw();
  c10a.SaveAs("ExOverEtrueVsETtrue_EB.png");
  TCanvas c10b("c10b","c10b",1);  
  myH->Draw();
  P_E25OverEtrueVsETtrue_EE  -> Draw("sameP");
  P_EneOverEtrueVsETtrue_EE  -> Draw("sameP");
  leg->Draw();
  c10b.SaveAs("ExOverEtrueVsETtrue_EE.png");



  // ------------------------------------------------------------------
  TLegend *leg2;
  leg2 = new TLegend(0.25,0.15,0.5,0.4);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.05);
  leg2->SetFillColor(0);
  leg2->AddEntry(h1_effHoE_EB,  "with cone", "p");
  leg2->AddEntry(h1_effHoET_EB, "with tower", "p");

  TH2F *myH2 = new TH2F("myH2","",100,0,4500,100,0.8,1.);
  myH2->GetXaxis()->SetTitle("ET true");
  myH2->GetYaxis()->SetTitle("H/E efficiency");

  TCanvas c11b("c11b","c11b",1);
  myH2->Draw();
  h1_effHoE_EB->Draw("sameP");
  h1_effHoET_EB->Draw("sameP");
  leg2->Draw();
  c11b.SaveAs("HoEcutEffVsET_barrel.png");

  TCanvas c11c("c11c","c11c",1);
  myH2->Draw();
  h1_effHoE_EE->Draw("sameP");
  h1_effHoET_EE->Draw("sameP");
  leg2->Draw();
  c11c.SaveAs("HoEcutEffVsET_endcap.png");

  TCanvas cc11b("cc11b","cc11b",1);
  myH2->Draw();
  h1nogap_effHoE_EB->Draw("sameP");
  h1nogap_effHoET_EB->Draw("sameP");
  leg2->Draw();
  cc11b.SaveAs("HoEcutEffVsET_gapsRemoved_barrel.png");

  TCanvas cc11c("cc11c","cc11c",1);
  myH2->Draw();
  h1nogap_effHoE_EE->Draw("sameP");
  h1nogap_effHoET_EE->Draw("sameP");
  leg2->Draw();
  cc11c.SaveAs("HoEcutEffVsET_gapsRemoved_endcap.png");


  TLegend *leg3;
  leg3 = new TLegend(0.25,0.15,0.5,0.4);
  leg3->SetFillStyle(0);
  leg3->SetBorderSize(0);
  leg3->SetTextSize(0.05);
  leg3->SetFillColor(0);
  leg3->AddEntry(h1_effHoEvsScEta,  "with cone",  "p");
  leg3->AddEntry(h1_effHoETvsScEta, "with tower", "p");

  TH2F *myH3 = new TH2F("myH3","",100,-2.6,2.6,100,0.8,1.);
  myH3->GetXaxis()->SetTitle("SC #eta");
  myH3->GetYaxis()->SetTitle("H/E efficiency");
  TH2F *myH3b = new TH2F("myH3b","",100,-2.6,2.6,100,0.8,1.);
  myH3b->GetXaxis()->SetTitle("SC #eta");
  myH3b->GetYaxis()->SetTitle("H/Etrue efficiency");

  TCanvas c13a("c13a","c13a",1);
  myH3->Draw();
  h1_effHoEvsScEta->Draw("sameP");
  h1_effHoETvsScEta->Draw("sameP");
  leg3->Draw();
  c13a.SaveAs("HoEcutEffvsScEta_ETgt200.png");

  TCanvas cc13a("cc13a","cc13a",1);
  myH3->Draw();
  h1nogap_effHoEvsScEta->Draw("sameP");
  h1nogap_effHoETvsScEta->Draw("sameP");
  leg3->Draw();
  cc13a.SaveAs("HoEcutEffvsScEta_ETgt200_gapsRemoved.png");

  TCanvas c13b("c13b","c13b",1);
  myH3b->Draw();
  h1_effHoEtruevsScEta->Draw("sameP");
  h1_effHoEtrueTvsScEta->Draw("sameP");
  leg3->Draw();
  c13b.SaveAs("HoEtruecutEffvsScEta_ETgt200.png");

  TCanvas cc13b("cc13b","cc13b",1);
  myH3b->Draw();
  h1nogap_effHoEtruevsScEta->Draw("sameP");
  h1nogap_effHoEtrueTvsScEta->Draw("sameP");
  leg3->Draw();
  cc13b.SaveAs("HoEtruecutEffvsScEta_ETgt200_gapsRemoved.png");
  
  TH2F *myH4 = new TH2F("myH4","",100,-3.15,3.15,100,0.8,1.);
  myH4->GetXaxis()->SetTitle("SC #phi");
  myH4->GetYaxis()->SetTitle("H/E efficiency");
  TH2F *myH4b = new TH2F("myH4b","",100,-3.15,3.15,100,0.8,1.);
  myH4b->GetXaxis()->SetTitle("SC #phi");
  myH4b->GetYaxis()->SetTitle("H/Etrue efficiency");

  TCanvas c14a("c14a","c14a",1);
  myH4->Draw();
  h1_effHoEvsScPhi->Draw("sameP");
  h1_effHoETvsScPhi->Draw("sameP");
  leg3->Draw();
  c14a.SaveAs("HoEcutEffvsScPhi_ETgt200.png");

  TCanvas cc14a("cc14a","cc14a",1);
  myH4->Draw();
  h1nogap_effHoEvsScPhi->Draw("sameP");
  h1nogap_effHoETvsScPhi->Draw("sameP");
  leg3->Draw();
  cc14a.SaveAs("HoEcutEffvsScPhi_ETgt200_gapsRemoved.png");

  TCanvas c14b("c14b","c14b",1);
  myH4b->Draw();
  h1_effHoEtruevsScPhi->Draw("sameP");
  h1_effHoEtrueTvsScPhi->Draw("sameP");
  leg3->Draw();
  c14b.SaveAs("HoEtruecutEffvsScPhi_ETgt200.png");

  TCanvas cc14b("cc14b","cc14b",1);
  myH4b->Draw();
  h1nogap_effHoEtruevsScPhi->Draw("sameP");
  h1nogap_effHoEtrueTvsScPhi->Draw("sameP");
  leg3->Draw();
  cc14b.SaveAs("HoEtruecutEffvsScPhi_ETgt200_gapsRemoved.png");

  TCanvas c14c("c14c","c14c",1);
  TH2F *myH5 = new TH2F("myH5","",100,-2.5,2.5,100,0.,25);
  myH5->GetXaxis()->SetTitle("SC #eta");
  myH5->GetYaxis()->SetTitle("H");
  myH5->Draw();
  P_HvsScEta->Draw("sameP");
  c14c.SaveAs("HvsScEta_ETgt200.png");

  TCanvas c14d("c14d","c14d",1);
  TH2F *myH5b = new TH2F("myH5b","",100,-3.14,3.14,100,0.,25);
  myH5b->GetXaxis()->SetTitle("SC #phi");
  myH5b->GetYaxis()->SetTitle("H");
  myH5b->Draw();
  P_HvsScPhi->Draw("sameP");
  c14d.SaveAs("HvsScPhi_ETgt200.png");

  TLegend *leg4;
  leg4 = new TLegend(0.25,0.15,0.5,0.4);
  leg4->SetFillStyle(0);
  leg4->SetBorderSize(0);
  leg4->SetTextSize(0.05);
  leg4->SetFillColor(0);
  leg4->AddEntry(H_EneOverEtrue_okHoE_EB, "H/E<0.05", "l");
  leg4->AddEntry(H_EneOverEtrue_noHoE_EB, "H/E>0.05", "l");
  
  TCanvas c15a("c15a","c15a",1);
  H_EneOverEtrue_okHoE_EB ->SetLineColor(2);
  H_EneOverEtrue_okHoE_EB ->SetLineWidth(2);
  H_EneOverEtrue_noHoE_EB ->SetLineColor(4);
  H_EneOverEtrue_noHoE_EB ->SetLineWidth(2);
  H_EneOverEtrue_okHoE_EB -> DrawNormalized("hist");
  H_EneOverEtrue_noHoE_EB -> DrawNormalized("samehist");
  leg4->Draw();
  c15a.SetLogy();
  c15a.SaveAs("EoverEtrue_passFailHoE_EB.png");

  TCanvas cc15a("cc15a","cc15a",1);
  Hnogap_EneOverEtrue_okHoE_EB ->SetLineColor(2);
  Hnogap_EneOverEtrue_okHoE_EB ->SetLineWidth(2);
  Hnogap_EneOverEtrue_noHoE_EB ->SetLineColor(4);
  Hnogap_EneOverEtrue_noHoE_EB ->SetLineWidth(2);
  Hnogap_EneOverEtrue_okHoE_EB -> DrawNormalized("hist");
  Hnogap_EneOverEtrue_noHoE_EB -> DrawNormalized("samehist");
  leg4->Draw();
  cc15a.SetLogy();
  cc15a.SaveAs("EoverEtrue_passFailHoE_EB_gapsRemoved.png");

  TCanvas c15b("c15b","c15b",1);
  H_EneOverEtrue_okHoE_EE ->SetLineColor(2);
  H_EneOverEtrue_okHoE_EE ->SetLineWidth(2);
  H_EneOverEtrue_noHoE_EE ->SetLineColor(4);
  H_EneOverEtrue_noHoE_EE ->SetLineWidth(2);
  H_EneOverEtrue_okHoE_EE -> DrawNormalized("hist");
  H_EneOverEtrue_noHoE_EE -> DrawNormalized("samehist");
  leg4->Draw();
  c15b.SetLogy();
  c15b.SaveAs("EoverEtrue_passFailHoE_EE.png");

  TCanvas cc15b("cc15b","cc15b",1);
  Hnogap_EneOverEtrue_okHoE_EE ->SetLineColor(2);
  Hnogap_EneOverEtrue_okHoE_EE ->SetLineWidth(2);
  Hnogap_EneOverEtrue_noHoE_EE ->SetLineColor(4);
  Hnogap_EneOverEtrue_noHoE_EE ->SetLineWidth(2);
  Hnogap_EneOverEtrue_okHoE_EE -> DrawNormalized("hist");
  Hnogap_EneOverEtrue_noHoE_EE -> DrawNormalized("samehist");
  leg4->Draw();
  cc15b.SetLogy();
  cc15b.SaveAs("EoverEtrue_passFailHoE_EE_gapsRemoved.png");

  TCanvas c15c("c15c","c15c",1);
  H_EneOverEtrue_okHoET_EB ->SetLineColor(2);
  H_EneOverEtrue_okHoET_EB ->SetLineWidth(2);
  H_EneOverEtrue_noHoET_EB ->SetLineColor(4);
  H_EneOverEtrue_noHoET_EB ->SetLineWidth(2);
  H_EneOverEtrue_okHoET_EB -> DrawNormalized("hist");
  H_EneOverEtrue_noHoET_EB -> DrawNormalized("samehist");
  leg4->Draw();
  c15c.SetLogy();
  c15c.SaveAs("EoverEtrue_passFailHoEwithTower_EB.png");

  TCanvas cc15c("cc15c","cc15c",1);
  Hnogap_EneOverEtrue_okHoET_EB ->SetLineColor(2);
  Hnogap_EneOverEtrue_okHoET_EB ->SetLineWidth(2);
  Hnogap_EneOverEtrue_noHoET_EB ->SetLineColor(4);
  Hnogap_EneOverEtrue_noHoET_EB ->SetLineWidth(2);
  Hnogap_EneOverEtrue_okHoET_EB -> DrawNormalized("hist");
  Hnogap_EneOverEtrue_noHoET_EB -> DrawNormalized("samehist");
  leg4->Draw();
  cc15c.SetLogy();
  cc15c.SaveAs("EoverEtrue_passFailHoEwithTower_EB_gapsRemoved.png");

  TCanvas c15d("c15d","c15d",1);
  H_EneOverEtrue_okHoET_EE ->SetLineColor(2);
  H_EneOverEtrue_okHoET_EE ->SetLineWidth(2);
  H_EneOverEtrue_noHoET_EE ->SetLineColor(4);
  H_EneOverEtrue_noHoET_EE ->SetLineWidth(2);
  H_EneOverEtrue_okHoET_EE -> DrawNormalized("hist");
  H_EneOverEtrue_noHoET_EE -> DrawNormalized("samehist");
  leg4->Draw();
  c15d.SetLogy();
  c15d.SaveAs("EoverEtrue_passFailHoEwithTower_EE.png");

  TCanvas cc15d("cc15d","cc15d",1);
  Hnogap_EneOverEtrue_okHoET_EE ->SetLineColor(2);
  Hnogap_EneOverEtrue_okHoET_EE ->SetLineWidth(2);
  Hnogap_EneOverEtrue_noHoET_EE ->SetLineColor(4);
  Hnogap_EneOverEtrue_noHoET_EE ->SetLineWidth(2);
  Hnogap_EneOverEtrue_okHoET_EE -> DrawNormalized("hist");
  Hnogap_EneOverEtrue_noHoET_EE -> DrawNormalized("samehist");
  leg4->Draw();
  cc15d.SetLogy();
  cc15d.SaveAs("EoverEtrue_passFailHoEwithTower_EE_gapsRemoved.png");


  
  /*
  TLegend *leg5;
  leg5 = new TLegend(0.25,0.15,0.5,0.4);
  leg5->SetFillStyle(0);
  leg5->SetBorderSize(0);
  leg5->SetTextSize(0.05);
  leg5->SetFillColor(0);
  leg5->AddEntry(h1_effSee_EB, "EB", "p");
  leg5->AddEntry(h1_effSee_EE, "EE", "p");

  TH2F *myH2b = new TH2F("myH2b","",100,0,4500,100,0.9,1.);
  myH2b->GetXaxis()->SetTitle("ET true");
  myH2b->GetYaxis()->SetTitle("SigmaIeIe cut efficiency");

  TCanvas c16a("c16a","c16a",1);
  myH2b->Draw();
  h1_effSee_EB->Draw("sameP");
  h1_effSee_EE->Draw("sameP");
  leg5->Draw();
  c16a.SaveAs("SeecutEffVsET.png");

  TH2F *myH2c = new TH2F("myH2c","",100,0,4500,100,0.8,1.);
  myH2c->GetXaxis()->SetTitle("ET true");
  myH2c->GetYaxis()->SetTitle("chargedH isolation cut efficiency");

  TCanvas c16b("c16b","c16b",1);
  myH2c->Draw();
  h1_effChIso_EB->Draw("sameP");
  h1_effChIso_EE->Draw("sameP");
  leg5->Draw();
  c16b.SaveAs("ChHadronIsocutEffVsET.png");

  TH2F *myH2d = new TH2F("myH2d","",100,0,4500,100,0.9,1.);
  myH2d->GetXaxis()->SetTitle("ET true");
  myH2d->GetYaxis()->SetTitle("neutralH isolation cut efficiency");

  TCanvas c16c("c16c","c16c",1);
  myH2d->Draw();
  h1_effNhIso_EB->Draw("sameP");
  h1_effNhIso_EE->Draw("sameP");
  leg5->Draw();
  c16c.SaveAs("NHadronIsocutEffVsET.png");
  */

  // summary
  h1_effSee_EB  ->SetMarkerColor(2);
  h1_effChIso_EB->SetMarkerColor(3);
  h1_effNhIso_EB->SetMarkerColor(4);
  h1_effHoET_EB ->SetMarkerColor(7);
  h1_effSee_EE  ->SetMarkerColor(2);
  h1_effChIso_EE->SetMarkerColor(3);
  h1_effNhIso_EE->SetMarkerColor(4);
  h1_effHoET_EE ->SetMarkerColor(7);

  TLegend *leg6;
  leg6 = new TLegend(0.5,0.15,0.75,0.35);
  leg6->SetFillStyle(0);
  leg6->SetBorderSize(0);
  leg6->SetTextSize(0.05);
  leg6->SetFillColor(0);
  leg6->AddEntry(h1_effSee_EB, "sIeIe", "p");
  leg6->AddEntry(h1_effChIso_EB,"chH iso", "p"); 
  leg6->AddEntry(h1_effNhIso_EB,"nH iso", "p"); 
  leg6->AddEntry(h1_effHoET_EB,"H/E", "p"); 

  TH2F *myH2r = new TH2F("myH2r","",100,0,4500,100,0.8,1.);
  myH2r->GetXaxis()->SetTitle("ET true");
  myH2r->GetYaxis()->SetTitle("2012 cut efficiency");

  TCanvas c17("c17","c17",1);
  myH2r->Draw();
  h1_effSee_EB->Draw("sameP");
  h1_effChIso_EB->Draw("sameP");
  h1_effNhIso_EB->Draw("sameP");  
  h1_effHoET_EB->Draw("sameP");
  leg6->Draw();
  c17.SaveAs("2012EffVsET_EB.png");

  TCanvas c17b("c17b","c17b",1);
  myH2r->Draw();
  h1_effSee_EE->Draw("sameP");
  h1_effChIso_EE->Draw("sameP");
  h1_effNhIso_EE->Draw("sameP");  
  h1_effHoET_EE->Draw("sameP");
  leg6->Draw();
  c17b.SaveAs("2012EffVsET_EE.png");
}
