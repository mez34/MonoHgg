#include "TStyle.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include <TGraphAsymmErrors.h>
#include "TCanvas.h"
#include "TLegend.h"

#include <iostream>

void selectionEfficiencyPlots() {

  // taking inputs
  // TFile *infile = new TFile("data/fullSel/mergedFinal/GGJets.root");   
  // TFile *infile = new TFile("data/fullSel/mergedFinal/GJets.root");
  // TFile *infile = new TFile("data/fullSel/mergedFinal/RSGravToGG_kMpl-01_M-1500.root");   
  TFile *infile = new TFile("data/fullSel/mergedFinal/RSGravToGG_kMpl-01_M-3000.root");   
  // 
  // TFile *infile = new TFile("data/noNHiso/mergedFinal/GGJets.root");   

  TH1F *h_denomPlusPres0 = (TH1F*)infile->Get("h_denomPlusPres0");
  TH1F *h_denomPlusPres1 = (TH1F*)infile->Get("h_denomPlusPres1");
  TH1F *h_denomPlusPres2 = (TH1F*)infile->Get("h_denomPlusPres2");
  TH1F *h_denomPlusPres3 = (TH1F*)infile->Get("h_denomPlusPres3");

  TH1F *h_num0 = (TH1F*)infile->Get("h_num0");
  TH1F *h_num1 = (TH1F*)infile->Get("h_num1");
  TH1F *h_num2 = (TH1F*)infile->Get("h_num2");
  TH1F *h_num3 = (TH1F*)infile->Get("h_num3");  
  
  TH1F *h_selection = (TH1F*)infile->Get("h_selection");

  // bins
  h_denomPlusPres0->Rebin(4);
  h_denomPlusPres1->Rebin(4);
  h_denomPlusPres2->Rebin(4);
  h_denomPlusPres3->Rebin(4);
  h_num0->Rebin(4);
  h_num1->Rebin(4);
  h_num2->Rebin(4);
  h_num3->Rebin(4);

  // efficiency distributions
  TGraphAsymmErrors *h_effVsPres0 = new TGraphAsymmErrors(h_num0,h_denomPlusPres0);
  TGraphAsymmErrors *h_effVsPres1 = new TGraphAsymmErrors(h_num1,h_denomPlusPres1);
  TGraphAsymmErrors *h_effVsPres2 = new TGraphAsymmErrors(h_num2,h_denomPlusPres2);
  TGraphAsymmErrors *h_effVsPres3 = new TGraphAsymmErrors(h_num3,h_denomPlusPres3);

  h_effVsPres0->SetTitle("EB, highR9");
  h_effVsPres1->SetTitle("EB, lowR9");
  h_effVsPres2->SetTitle("EE, highR9");
  h_effVsPres3->SetTitle("EE, lowR9");

  h_effVsPres0->GetXaxis()->SetTitle("m(#gamma#gamma)");
  h_effVsPres1->GetXaxis()->SetTitle("m(#gamma#gamma)");
  h_effVsPres2->GetXaxis()->SetTitle("m(#gamma#gamma)");
  h_effVsPres3->GetXaxis()->SetTitle("m(#gamma#gamma)");

  // cosmetics
  h_effVsPres0->SetMarkerStyle(20);
  h_effVsPres1->SetMarkerStyle(20);
  h_effVsPres2->SetMarkerStyle(20);
  h_effVsPres3->SetMarkerStyle(20);

  h_effVsPres0->SetMarkerColor(2);
  h_effVsPres1->SetMarkerColor(2);
  h_effVsPres2->SetMarkerColor(2);
  h_effVsPres3->SetMarkerColor(2);


  // breakdown
  h_selection->SetLineColor(4);
  h_selection->SetLineWidth(2);
  h_selection->GetXaxis()->SetBinLabel(1,"no cut");
  h_selection->GetXaxis()->SetBinLabel(2,"2 preselected #gamma");
  h_selection->GetXaxis()->SetBinLabel(3,"2 selected #gamma");
  h_selection->GetXaxis()->SetBinLabel(4,"pT cuts");
  h_selection->GetXaxis()->SetBinLabel(5,"good vertex");
  h_selection->GetXaxis()->SetBinLabel(6,"m(#gamma#gamma) cut");

  // plots
  gStyle->SetOptStat(0);

  TH2F *myH2  = new TH2F("myH2", "",100,0.,6000.,100,0.5,1.);
  myH2->GetXaxis()->SetTitle("m(#gamma#gamma)");
  myH2->GetYaxis()->SetTitle("efficiency");

  TCanvas ca("ca","",1);
  myH2->SetTitle("EB, high R9");
  myH2->Draw();
  h_effVsPres0->Draw("sameP");
  ca.SaveAs("effVsMass_class0.png");

  TCanvas cb("cb","",1);
  myH2->SetTitle("EB, low R9");
  myH2->Draw();
  h_effVsPres1->Draw("sameP");
  cb.SaveAs("effVsMass_class1.png");

  TCanvas cc("cc","",1);
  myH2->SetTitle("EE, high R9");
  myH2->Draw();
  h_effVsPres2->Draw("sameP");
  cc.SaveAs("effVsMass_class2.png");

  TCanvas cd("cd","",1);
  myH2->SetTitle("EE, low R9");
  myH2->Draw();
  h_effVsPres3->Draw("sameP");
  cd.SaveAs("effVsMass_class3.png");
  
  TCanvas ce("ce","",1);
  h_selection->Draw("hist");
  h_selection->SetTitle("");
  // h_selection->SetMinimum(55000);
  // h_selection->SetMaximum(105000);
  ce.SaveAs("breakdown.png");
  ce.SetLogy();
  ce.SaveAs("breakdown_log.png");
}
