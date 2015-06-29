#include <iostream>
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TString.h"
#include "TLegend.h"

void combiner( char * in ){
  TString in_name = in;
  TString fIn1 = "DMHgg50k"; //signal
  TString fIn2 = "WZHtoGG"; 
  TString fIn3 = "QCD1";
  TString fIn4 = "GJets1";
  
  overlay(fIn1, fIn2, fIn3, fIn4, in_name);
}
void overlay( const TString In1, const TString In2, const TString In3, const TString In4, const TString name ){
  TString location = "diPhoPlots/";

  gROOT->SetStyle("Plain");
  gROOT->ProcessLine(".x ./DiphotonStyle.C");
  
  TFile *f1 = new TFile(location+In1+"/diPhotHistos_"+In1+".root");
  TFile *f2 = new TFile(location+In2+"/diPhotHistos_"+In2+".root");
  TFile *f3 = new TFile(location+In3+"/diPhotHistos_"+In3+".root");
  TFile *f4 = new TFile(location+In4+"/diPhotHistos_"+In4+".root");

  TH1D *h_1_1 = (TH1D*) f1->Get(name+"_"+In1);
  TH1D *h_1_2 = (TH1D*) f2->Get(name+"_"+In2);
  TH1D *h_1_3 = (TH1D*) f3->Get(name+"_"+In3);
  TH1D *h_1_4 = (TH1D*) f4->Get(name+"_"+In4);

  TH1D *h1[4];
  h1[0]=h_1_1;
  h1[1]=h_1_2;
  h1[2]=h_1_3;
  h1[3]=h_1_4;

  TCanvas *c1 = new TCanvas("c1","",1200,800);
  c1->cd();
  c1->SetLogy(1);

  h1[0]->SetTitle(name);
  h1[0]->DrawNormalized();
  h1[1]->SetLineColor(kRed);
  h1[1]->DrawNormalized("same");
  h1[2]->SetLineColor(kGreen+2);
  h1[2]->DrawNormalized("same");
  h1[3]->SetLineColor(kBlue);
  h1[3]->DrawNormalized("same");

  TLegend *l1 = new TLegend(0.75,0.75,0.95,0.92);
  l1->SetTextFont(42);
  l1->SetFillColor(0);
  l1->SetBorderSize(2);
  l1->AddEntry(h1[0],"DMHgg","lp");
  l1->AddEntry(h1[1],"WZHtoGG","lp");
  l1->AddEntry(h1[2],"QCD","lp");
  l1->AddEntry(h1[3],"GJets","lp");
  l1->Draw();

  c1->SaveAs("diPhoPlots/"+name+"_comb.png");

}
