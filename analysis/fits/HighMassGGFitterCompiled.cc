#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooCBShape.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooBreitWigner.h"
#include "RooFFTConvPdf.h"
#include "RooHistFunc.h"
#include "RooHistPdf.h"
#include "RooFitResult.h"
#include "RooStats/HLFactory.h"

#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphAsymmErrors.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPaveText.h"

#include <fstream>

using namespace RooFit;
using namespace RooStats;

// ============================================
// to be modified:
static const Int_t NCAT = 4;  
Int_t MINmass= 500;
Int_t MAXmass= 6000;
Float_t Lum = 19500.0;    
bool wantResponse = 0;
bool wantGenLevel = 0;
// ============================================


// Definition of the variables in the input ntuple
RooArgSet* defineVariables() {

  RooRealVar* mgg        = new RooRealVar("mgg",        "M(gg)",       MINmass, MAXmass, "GeV");
  RooRealVar* mggGen     = new RooRealVar("mggGen",     "M(gg) gen",   MINmass, MAXmass, "GeV");
  RooRealVar* eventClass = new RooRealVar("eventClass", "eventClass",    -10,      10,   "");
  // RooRealVar* weight     = new RooRealVar("weight",     "weightings",      0,     1000,  "");   // chiara

  // RooArgSet* ntplVars = new RooArgSet(*mgg, *mggGen, *eventClass, *weight);      
  RooArgSet* ntplVars = new RooArgSet(*mgg, *mggGen, *eventClass);
  
  return ntplVars;
}

void SetConstantParams(const RooArgSet* params) {
  
  cout << endl; cout << "Entering SetConstantParams" << endl;
  TIterator* iter(params->createIterator());
  for (TObject *a = iter->Next(); a != 0; a = iter->Next()) {
    RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
    if (rrv) { rrv->setConstant(true); std::cout << " " << rrv->GetName(); }
  }  
}

// CMS stuffs
TPaveText* get_labelCMS( int legendQuadrant, std::string year, bool sim) {
  
  if( legendQuadrant!=0 && legendQuadrant!=1 && legendQuadrant!=2 && legendQuadrant!=3 ) {
    std::cout << "WARNING! Legend quadrant '" << legendQuadrant << "' not yet implemented for CMS label. Using 2." << std::endl;
    legendQuadrant = 2;
  }
  
  float x1=0.;
  float y1=0.;
  float x2=0.;
  float y2=0.;
  if( legendQuadrant==1 ) {
    x1 = 0.63;
    y1 = 0.83;
    x2 = 0.8;
    y2 = 0.87;
  } else if( legendQuadrant==0 ) {
    x1 = 0.175;
    y1 = 0.953;
    x2 = 0.6;
    y2 = 0.975;
  } else if( legendQuadrant==2 ) {
    x1 = 0.25;
    y1 = 0.2;
    x2 = 0.42;
  }
 
  TPaveText* cmslabel = new TPaveText( x1, y1, x2, y2, "" );
  cmslabel->SetFillColor(kWhite);
  cmslabel->SetTextSize(0.038);
  if( legendQuadrant==0 ) cmslabel->SetTextAlign(11);
  cmslabel->SetTextSize(0.038);
  cmslabel->SetTextFont(42);
  
  std::string leftText;
   
  if (sim)  leftText = "CMS Simulation"; 
  else {
    leftText = "CMS Preliminary, xxx fb^{-1}";
  }
  cmslabel->AddText(leftText.c_str());
  return cmslabel;
}

TPaveText* get_labelSqrt( int legendQuadrant ) {

  if( legendQuadrant!=0 && legendQuadrant!=1 && legendQuadrant!=2 ) {
    std::cout << "WARNING! Legend quadrant '" << legendQuadrant << "' not yet implemented for Sqrt label. Using 2." << std::endl;
    legendQuadrant = 2;
  }

  float x1=0.;
  float y1=0.;
  float x2=0.;
  float y2=0.;
  if( legendQuadrant==1 ) {
    x1 = 0.63;
    y1 = 0.78;
    x2 = 0.8;
    y2 = 0.82;
  } else if( legendQuadrant==2 ) {
    x1 = 0.25;
    y1 = 0.78;
    x2 = 0.42;
    y2 = 0.82;
  } else if( legendQuadrant==0 ) {
    x1 = 0.65;
    y1 = 0.953;
    x2 = 0.87;
    y2 = 0.975;
  }

  TPaveText* label_sqrt = new TPaveText(x1,y1,x2,y2, "");
  label_sqrt->SetFillColor(kWhite);
  label_sqrt->SetTextSize(0.038);
  label_sqrt->SetTextFont(42);
  label_sqrt->SetTextAlign(31); 
  label_sqrt->AddText("#sqrt{s} = 13 TeV");
  return label_sqrt;
}

// loading signal data and making roodatasets
void AddSigData(RooWorkspace* w, Float_t mass, TString coupling) {
  
  Int_t ncat = NCAT;
  
  // Variables
  RooArgSet* ntplVars = defineVariables();

  // -------------------------  
  // Files
  int iMass = abs(mass);   
  TString inDir = "data/newSelection/mergedFinal/";
  TChain* sigTree = new TChain();
  cout << "reading file " << inDir+TString(Form("RSGravToGG_kMpl-"))+coupling+TString(Form("_M-%d.root/DiPhotonTree", iMass)) << endl;
  sigTree->Add(inDir+TString(Form("RSGravToGG_kMpl-"))+coupling+TString(Form("_M-%d.root/DiPhotonTree", iMass)));
  sigTree->SetTitle("sigTree");
  sigTree->SetName("sigTree");


  // -------------------------
  // common preselection cut on mgg and mggGen
  TString mainCut1 = TString::Format("mgg>=500 && mgg<=6000 && mggGen>=500 && mggGen<=6000");   
  // RooDataSet sigWeighted("sigWeighted","dataset",sigTree,*ntplVars,mainCut1,"weight");   // chiara!
  RooDataSet sigWeighted("sigWeighted","dataset",sigTree,*ntplVars,mainCut1);


  // -------------------------
  // reduced mass
  RooFormulaVar *massReduced_formula = new RooFormulaVar("massReduced_formula","","@0/@1",RooArgList(*w->var("mgg"),*w->var("mggGen")));
  RooRealVar* massReduced = (RooRealVar*) sigWeighted.addColumn(*massReduced_formula);
  massReduced->SetName("massReduced");
  massReduced->SetTitle("massReduced");
  w->import(*massReduced);  
  massReduced->setRange(0.5, 1.5);
  
  // common preselection cut on the reduced mass
  TString mainCut = TString::Format("massReduced>0.5 && massReduced <1.5"); 

  
  // -------------------------
  // split in categories, wrt mgg - this is the dataset to be used for the convolution
  cout << endl;
  cout << "preparing dataset with observable mgg" << endl;
  RooDataSet* signal[NCAT];
  for (int c=0; c<ncat; ++c) {
    if (c==0) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("mgg"),mainCut+TString::Format("&& eventClass==0"));
    if (c==1) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("mgg"),mainCut+TString::Format("&& eventClass==1"));
    if (c==2) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("mgg"),mainCut+TString::Format("&& eventClass==2"));
    if (c==3) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("mgg"),mainCut+TString::Format("&& eventClass==3"));

    w->import(*signal[c],Rename(TString::Format("SigWeight_cat%d",c)));
    
    cout << "cat " << c << ", signal[c]: " << endl;
    signal[c]->Print("v");
    cout << "---- for category " << c << ", nX for signal[c]:  " << signal[c]->sumEntries() << endl; 
    cout << endl;
  }

  // Create full weighted signal data set without categorization
  RooDataSet* signalAll = (RooDataSet*) sigWeighted.reduce(*w->var("mgg"),mainCut);
  w->import(*signalAll, Rename("SigWeight"));
  cout << "now signalAll" << endl;
  signalAll->Print("v");
  cout << "---- nX for signalAll:  " << signalAll->sumEntries() << endl; 
  cout << endl;


  // -------------------------
  // split in categories, wrt massReduced - to study the detector response
  if (wantResponse) {
    cout << endl;
    cout << endl;
    cout << "preparing dataset with observable massReduced" << endl;
    RooDataSet* signalR[NCAT];
    for (int c=0; c<ncat; ++c) {
      if (c==0) signalR[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& eventClass==0"));
      if (c==1) signalR[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& eventClass==1"));
      if (c==2) signalR[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& eventClass==2"));
      if (c==3) signalR[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& eventClass==3"));
      w->import(*signalR[c],Rename(TString::Format("SigWeightReduced_cat%d",c)));
    }
    cout << endl;
  }


  // -------------------------
  // split in categories, wrt genMass - to study the theory width
  if (wantGenLevel) {   // chiara: se lo usi nei toys sistema questo switch
    cout << endl;
    cout << endl;
    cout << "preparing dataset with observable mggGen, no split in categories since they're all the same" << endl;
    RooDataSet* signalGAll = (RooDataSet*) sigWeighted.reduce(*w->var("mggGen"),mainCut);
    w->import(*signalGAll, Rename("SigWeightGen"));
    cout << endl;
  }

  cout << "workspace summary" << endl;
  w->Print();
}

// Signal model: doubleCB. To describe the detector resolution, reco mgg used
void SigModelResponseDoubleCBFit(RooWorkspace* w, Float_t mass, TString coupling) {
  
  int iMass = abs(mass);   
  
  // dataset
  RooDataSet* signal[NCAT];
  RooRealVar* mgg = w->var("mgg");     
  
  // fit function
  RooDoubleCB* ResponseDoubleCB[NCAT];

  TCanvas* c1 = new TCanvas("c1", "c1", 1);
  c1->cd();

  for(int c = 0; c<NCAT; c++){
    
    // loading the dataset
    signal[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c)); 

    // DoubleCB
    RooFormulaVar CBmean(TString::Format("CB_sig_mean_cat%d",c),"", "@0", *w->var(TString::Format("Sig_mean_cat%d",c)));
    RooFormulaVar CBsigma(TString::Format("CB_sig_sigma_cat%d",c), "", "sqrt(@0*@0)", *w->var(TString::Format("Sig_sigma_cat%d",c)));
    RooFormulaVar CBalpha1(TString::Format("CB_sig_alpha1_cat%d",c),"", "@0", *w->var( TString::Format("Sig_alpha1_cat%d",c)));
    RooFormulaVar CBn1(TString::Format("CB_sig_n1cat%d",c),"", "@0", *w->var( TString::Format("Sig_n1_cat%d",c)));
    RooFormulaVar CBalpha2(TString::Format("CB_sig_alpha2_cat%d",c),"", "@0", *w->var( TString::Format("Sig_alpha2_cat%d",c)));
    RooFormulaVar CBn2(TString::Format("CB_sig_n2cat%d",c),"", "@0", *w->var( TString::Format("Sig_n2_cat%d",c)));
    ResponseDoubleCB[c] = new RooDoubleCB(TString::Format("ResponseDoubleCB_cat%d",c),TString::Format("ResponseDoubleCB_cat%d",c) , *mgg, CBmean, CBsigma, CBalpha1, CBn1, CBalpha2, CBn2) ;    
    w->import(*ResponseDoubleCB[c]);

    // Fit with ResponseDoubleCB
    RooFitResult* fitresults = 0;
    if (mass==1500) 
      fitresults = (RooFitResult* ) ResponseDoubleCB[c]->fitTo(*signal[c], SumW2Error(kFALSE), Range(1250, 1700), RooFit::Save(kTRUE));
    else if (mass==750) 
      fitresults = (RooFitResult* ) ResponseDoubleCB[c]->fitTo(*signal[c], SumW2Error(kFALSE), Range(500, 1200), RooFit::Save(kTRUE));
    else if (mass==5000) 
      fitresults = (RooFitResult* ) ResponseDoubleCB[c]->fitTo(*signal[c], SumW2Error(kFALSE), Range(4000, 5500), RooFit::Save(kTRUE));
    
    std::cout<<TString::Format("******************************** Signal Fit results DoubleCB mass %f cat %d***********************************", mass, c)<<std::endl;
    fitresults->Print("V");
   


    // Plot
    RooPlot* plotG = mgg->frame(Range(1250, 1650),Title("Mgg, response"),Bins(60));
    plotG->GetXaxis()->SetTitle("m_{reco}");
    plotG->GetXaxis()->SetTitleFont(42);
    plotG->GetXaxis()->SetTitleSize(0.04);
    plotG->GetXaxis()->SetTitleOffset(1.40);
    plotG->GetYaxis()->SetTitle("Events");
    plotG->GetYaxis()->SetTitleFont(42);
    plotG->GetYaxis()->SetTitleSize(0.04);

    TLegend* legmc = new TLegend(0.6, 0.58, 0.91, 0.91, "", "");
    legmc->SetTextSize(0.0286044);  
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->AddEntry(plotG->getObject(0),"Simulation","LP");    
    legmc->AddEntry(plotG->getObject(1),"DoubleCB ","L");
    
    TLatex* latex = new TLatex(0.21, 0.76, TString::Format("#splitline{m_{X}=%d GeV}{#splitline{}{Class %d}}",iMass,c));
    latex->SetTextSize(0.038);
    latex->SetTextAlign(11);
    latex->SetTextFont(42); 
    latex->SetNDC();

    // Lin scale
    if (mass==1500) {
      plotG = mgg->frame(Range(1400, 1550),Title("Mgg, response"),Bins(60));           
    } else if (mass==750) {
      plotG = mgg->frame(Range(710, 780),Title("Mgg, response"),Bins(60));          
    } else if (mass==5000) {
      plotG = mgg->frame(Range(4500, 5300),Title("Mgg, response"),Bins(60));          
    }
    signal[c]->plotOn(plotG);
    ResponseDoubleCB[c]->plotOn(plotG, LineColor(kBlue));
    plotG->Draw();
    latex->Draw("same");
    legmc->Draw("same");
    c1->SetLogy(0);
    c1->SaveAs(TString::Format("plots/responseDoubleCB_cat%d.png",c));

    // The Log scale
    if (mass==1500) {
      plotG = mgg->frame(Range(1250, 1650),Title("Mgg, response"),Bins(60));        
    } else if (mass==750) {
      plotG = mgg->frame(Range(550, 850),Title("Mgg, response"),Bins(60));             
    } else if (mass==5000) {
      plotG = mgg->frame(Range(4000, 5500),Title("Mgg, response"),Bins(60));           
    }
    signal[c]->plotOn(plotG);
    ResponseDoubleCB[c]->plotOn(plotG, LineColor(kBlue));
    plotG->Draw();
    latex->Draw("same");
    legmc->Draw("same");
    c1->SetLogy();
    c1->SaveAs(TString::Format("plots/responseDoubleCB_cat%d_LOG.png",c));


    // saving as constant in the WS    
    w->defineSet(TString::Format("ResponseDoubleCBPdfParam_cat%d",c),RooArgSet(*w->var(TString::Format("Sig_sigma_cat%d",c)), 
									       *w->var(TString::Format("Sig_alpha1_cat%d",c)),
									       *w->var(TString::Format("Sig_alpha2_cat%d",c)),
									       *w->var(TString::Format("Sig_n1_cat%d",c)),
									       *w->var(TString::Format("Sig_n2_cat%d",c)),	   
									       *w->var(TString::Format("Sig_mean_cat%d",c))));
    SetConstantParams(w->set(TString::Format("ResponseDoubleCBPdfParam_cat%d",c)));

    w->Print();
  }
}

// Signal model: sum of 2 CBs. To describe the detector resolution, reco mgg used
void SigModelResponseCBCBFit(RooWorkspace* w, Float_t mass, TString coupling) {
  
  int iMass = abs(mass);   
  
  // Dataset
  RooDataSet* signal[NCAT];
  RooRealVar* mgg = w->var("mgg");     

  // fit functions
  RooCBShape* ResponseCBpos[NCAT];
  RooCBShape* ResponseCBneg[NCAT];
  RooAddPdf* ResponseAdd[NCAT];
  
  TCanvas* c1 = new TCanvas("c1", "c1", 1);
  c1->cd();

  for(int c = 0; c<NCAT; c++){
    
    // taking the dataset
    signal[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c)); 

    // cb pos                                                               
    RooFormulaVar CBpos_mean(TString::Format("MassCB_sig_mean_cat%d",c),"", "@0", *w->var(TString::Format("Mass_sig_mean_cat%d",c)));
    RooFormulaVar CBpos_sigma(TString::Format("MassCB_sig_sigma_cat%d",c), "", "sqrt(@0*@0)", *w->var(TString::Format("Mass_sig_sigma_cat%d",c)));
    RooFormulaVar CBpos_alphaCB(TString::Format("MassCB_sig_alphaCBpos_cat%d",c),"", "@0", *w->var( TString::Format("Mass_sig_alphaCBpos_cat%d",c)));
    RooFormulaVar CBpos_n(TString::Format("MassCB_sig_nCBpos_cat%d",c),"", "@0", *w->var( TString::Format("Mass_sig_nCBpos_cat%d",c)));
    ResponseCBpos[c] = new RooCBShape(TString::Format("ResponseCBpos_cat%d",c),TString::Format("ResponseCBpos_cat%d",c) , *mgg, CBpos_mean, CBpos_sigma,  CBpos_alphaCB, CBpos_n) ;
    
    // cb neg (same mean and sigma)
    RooFormulaVar CBneg_n(TString::Format("MassCB_sig_nCBneg_cat%d",c),"", "@0", *w->var( TString::Format("Mass_sig_nCBneg_cat%d",c)));
    RooFormulaVar CBneg_alphaCB(TString::Format("MassCB_sig_alphaCBneg_cat%d",c),"", "@0", *w->var( TString::Format("Mass_sig_alphaCBneg_cat%d",c)));
    ResponseCBneg[c] = new RooCBShape(TString::Format("ResponseCBneg_cat%d",c),TString::Format("ResponseCBneg_cat%d",c) , *mgg, CBpos_mean, CBpos_sigma, CBneg_alphaCB, CBneg_n) ;
   
    // combination pos and neg
    RooFormulaVar CB_frac(TString::Format("MassCB_sig_frac_cat%d",c),"", "@0", *w->var(TString::Format("Mass_sig_frac_cat%d",c)));
    ResponseAdd[c]= new RooAddPdf(TString::Format("ResponseAddPdf_cat%d",c),TString::Format("ResponseAddPdf_cat%d",c) , RooArgList(*ResponseCBneg[c], *ResponseCBpos[c]), CB_frac);

    w->import(*ResponseAdd[c]);
   
    // Fit with ResponseAdd
    RooFitResult* fitresults = 0;
    if (mass==1500) 
      fitresults = (RooFitResult* ) ResponseAdd[c]->fitTo(*signal[c], SumW2Error(kFALSE), Range(1250, 1700), RooFit::Save(kTRUE));
    else if (mass==750) 
      fitresults = (RooFitResult* ) ResponseAdd[c]->fitTo(*signal[c], SumW2Error(kFALSE), Range(500, 1200), RooFit::Save(kTRUE));
    else if (mass==5000) 
      fitresults = (RooFitResult* ) ResponseAdd[c]->fitTo(*signal[c], SumW2Error(kFALSE), Range(4000, 5500), RooFit::Save(kTRUE));
    
    std::cout<<TString::Format("******************************** Signal Fit results CB+CB  mass %f cat %d***********************************", mass, c)<<std::endl;
    fitresults->Print("V");
   
    
    
    // Plot
    RooPlot* plotG = mgg->frame(Range(1250, 1650),Title("Mgg, response"),Bins(60));
    plotG->GetXaxis()->SetTitle("m_{reco}");
    plotG->GetXaxis()->SetTitleFont(42);
    plotG->GetXaxis()->SetTitleSize(0.04);
    plotG->GetXaxis()->SetTitleOffset(1.40);
    plotG->GetYaxis()->SetTitle("Events");
    plotG->GetYaxis()->SetTitleFont(42);
    plotG->GetYaxis()->SetTitleSize(0.04);
    TLegend* legmc = new TLegend(0.6, 0.58, 0.91, 0.91, "", "");
    legmc->SetTextSize(0.0286044);  
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->AddEntry(plotG->getObject(0),"Simulation","LP");    
    legmc->AddEntry(plotG->getObject(1),"Sum of two CB ","L");
    legmc->AddEntry(plotG->getObject(2),"CB 1","L");   
    legmc->AddEntry(plotG->getObject(3),"CB 2","L");   
    
    TLatex* latex = new TLatex(0.21, 0.76, TString::Format("#splitline{m_{X}=%d GeV}{#splitline{}{Class %d}}",iMass,c));
    latex->SetTextSize(0.038);
    latex->SetTextAlign(11);
    latex->SetTextFont(42); 
    latex->SetNDC();


    // First lin scale
    if (mass==1500) {
      plotG = mgg->frame(Range(1400, 1550),Title("Mgg, response"),Bins(60));           
    } else if (mass==750) {
      plotG = mgg->frame(Range(710, 780),Title("Mgg, response"),Bins(60));               
    } else if (mass==5000) {
      plotG = mgg->frame(Range(4500, 5300),Title("Mgg, response"),Bins(60));          
    }
    signal[c]->plotOn(plotG);
    ResponseAdd[c]->plotOn(plotG, LineColor(kBlue));
    ResponseAdd[c]->plotOn(plotG,Components(TString::Format("ResponseCBneg_cat%d",c)), LineColor(kGreen), LineStyle(kDotted));
    ResponseAdd[c]->plotOn(plotG,Components(TString::Format("ResponseCBpos_cat%d",c)), LineColor(kRed), LineStyle(kDashed));
    plotG->Draw();
    latex->Draw("same");
    legmc->Draw("same");
    c1->SetLogy(0);
    c1->SaveAs(TString::Format("plots/responseAbsoluteFitCBCB_cat%d.png",c));

    // The log scale
    if (mass==1500) {
      plotG = mgg->frame(Range(1250, 1650),Title("Mgg, response"),Bins(60));              
    } else if (mass==750) {
      plotG = mgg->frame(Range(550, 850),Title("Mgg, response"),Bins(60));             
    } else if (mass==5000) {
      plotG = mgg->frame(Range(4000, 5500),Title("Mgg, response"),Bins(60));           
    }
    signal[c]->plotOn(plotG);
    ResponseAdd[c]->plotOn(plotG, LineColor(kBlue));
    ResponseAdd[c]->plotOn(plotG,Components(TString::Format("ResponseCBneg_cat%d",c)), LineColor(kGreen), LineStyle(kDotted));
    ResponseAdd[c]->plotOn(plotG,Components(TString::Format("ResponseCBpos_cat%d",c)), LineColor(kRed), LineStyle(kDashed));
    plotG->Draw();
    latex->Draw("same");
    legmc->Draw("same");
    c1->SetLogy();
    c1->SaveAs(TString::Format("plots/responseAbsoluteFitCBCB_cat%d_LOG.png",c));


    // saving as constant in the WS
    w->defineSet(TString::Format("ResponseAddPdfParam_cat%d",c),RooArgSet(*w->var(TString::Format("Mass_sig_sigma_cat%d",c)), 
									  *w->var(TString::Format("Mass_sig_alphaCBpos_cat%d",c)),
									  *w->var(TString::Format("Mass_sig_alphaCBneg_cat%d",c)),
									  *w->var(TString::Format("Mass_sig_nCBpos_cat%d",c)),
									  *w->var(TString::Format("Mass_sig_nCBneg_cat%d",c)),	   
									  *w->var(TString::Format("Mass_sig_frac_cat%d",c)),  
									  *w->var(TString::Format("Mass_sig_mean_cat%d",c))));
    SetConstantParams(w->set(TString::Format("ResponseAddPdfParam_cat%d",c)));

    w->Print();
  }
}

// Signal model: sum of two CBs. Detector resolution function based on mreco/mgen -1
void SigModelResponseReducedCBCBFit(RooWorkspace* w, Float_t mass, TString coupling) {
  
  int iMass = abs(mass);   
  
  // Dataset 
  RooDataSet* signal[NCAT];
  RooRealVar* mgg = w->var("mgg");     
  RooRealVar* massReduced = w->var("massReduced");
    
  // fit functions
  RooCBShape* ResponseCBpos[NCAT];
  RooCBShape* ResponseCBneg[NCAT];
  RooAddPdf* ResponseAdd[NCAT];
  
  TCanvas* c1 = new TCanvas("c1", "c1", 1);
  c1->cd();
  
  for(int c = 0; c<NCAT; c++){

    // taking the dataset  
    signal[c] = (RooDataSet*) w->data(TString::Format("SigWeightReduced_cat%d",c));
    
    // cb pos                                                               
    RooFormulaVar CBpos_mean(TString::Format("ReducedMassCB_sig_mean_cat%d",c),"", "@0", *w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)));
    RooFormulaVar CBpos_sigma(TString::Format("ReducedMassCB_sig_sigma_cat%d",c), "", "sqrt(@0*@0)", *w->var(TString::Format("ReducedMass_sig_sigma_cat%d",c)));
    RooFormulaVar CBpos_alphaCB(TString::Format("ReducedMassCB_sig_alphaCBpos_cat%d",c),"", "@0", *w->var( TString::Format("ReducedMass_sig_alphaCBpos_cat%d",c)));
    RooFormulaVar CBpos_n(TString::Format("ReducedMassCB_sig_nCBpos_cat%d",c),"", "@0", *w->var( TString::Format("ReducedMass_sig_nCBpos_cat%d",c)));
    ResponseCBpos[c] = new RooCBShape(TString::Format("ResponseCBpos_cat%d",c),TString::Format("ResponseCBpos_cat%d",c) , *massReduced, CBpos_mean, CBpos_sigma,  CBpos_alphaCB, CBpos_n) ;
    
    // cb neg (same mean and sigma)
    RooFormulaVar CBneg_n(TString::Format("ReducedMassCB_sig_nCBneg_cat%d",c),"", "@0", *w->var( TString::Format("ReducedMass_sig_nCBneg_cat%d",c)));
    RooFormulaVar CBneg_alphaCB(TString::Format("ReducedMassCB_sig_alphaCBneg_cat%d",c),"", "@0", *w->var( TString::Format("ReducedMass_sig_alphaCBneg_cat%d",c)));
    ResponseCBneg[c] = new RooCBShape(TString::Format("ResponseCBneg_cat%d",c),TString::Format("ResponseCBneg_cat%d",c) , *massReduced, CBpos_mean, CBpos_sigma, CBneg_alphaCB, CBneg_n) ;
   
    // combination pos and neg
    RooFormulaVar CB_frac(TString::Format("ReducedMassCB_sig_frac_cat%d",c),"", "@0", *w->var(TString::Format("ReducedMass_sig_frac_cat%d",c)));
    ResponseAdd[c]= new RooAddPdf(TString::Format("ResponseAddPdf_cat%d",c),TString::Format("ResponseAddPdf_cat%d",c) , RooArgList(*ResponseCBneg[c], *ResponseCBpos[c]), CB_frac);

    w->import(*ResponseAdd[c]);
   
    // Fit with ResponseAdd
    RooFitResult* fitresults = (RooFitResult* ) ResponseAdd[c]->fitTo(*signal[c], SumW2Error(kFALSE), Range(0.95, 1.05), RooFit::Save(kTRUE));
    std::cout<<TString::Format("******************************** Signal Fit results CB+CB  mass %f cat %d***********************************", mass, c)<<std::endl;
    fitresults->Print("V");
   

    // Plot
    RooPlot* plotG = massReduced->frame(Range(0.95, 1.05),Title("Mass Reduced"),Bins(60));
    signal[c]->plotOn(plotG);

    ResponseAdd[c]->plotOn(plotG, LineColor(kBlue));
    ResponseAdd[c]->plotOn(plotG,Components(TString::Format("ResponseCBneg_cat%d",c)), LineColor(kGreen), LineStyle(kDotted));
    ResponseAdd[c]->plotOn(plotG,Components(TString::Format("ResponseCBpos_cat%d",c)), LineColor(kRed), LineStyle(kDashed));
  
    plotG->GetXaxis()->SetTitle("#frac{m_{reco}}{m_{true}}");
    plotG->GetXaxis()->SetTitleFont(42);
    plotG->GetXaxis()->SetTitleSize(0.04);
    plotG->GetXaxis()->SetTitleOffset(1.40);
    plotG->GetYaxis()->SetTitle("Events/0.0024 units");
    plotG->GetYaxis()->SetTitleFont(42);
    plotG->GetYaxis()->SetTitleSize(0.04);

    TLegend* legmc = new TLegend(0.6, 0.58, 0.91, 0.91, "", "");
    legmc->SetTextSize(0.0286044);  
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->AddEntry(plotG->getObject(0),"Simulation","LP");    
    legmc->AddEntry(plotG->getObject(1),"Sum of two CB ","L");
    legmc->AddEntry(plotG->getObject(2),"CB 1","L");   
    legmc->AddEntry(plotG->getObject(3),"CB 2","L");   
    
    TLatex* latex = new TLatex(0.21, 0.76, TString::Format("#splitline{m_{X}=%d GeV}{#splitline{}{Class %d}}",iMass,c));
    latex->SetTextSize(0.038);
    latex->SetTextAlign(11);
    latex->SetTextFont(42); 
    latex->SetNDC();
   
    plotG->Draw();
    
    latex->Draw("same");
    legmc->Draw("same");
    int iPos=11 ;

    c1->SetLogy(0);
    c1->SaveAs(TString::Format("plots/responseFitCBCB_cat%d.png",c));
    c1->SetLogy();
    c1->SaveAs(TString::Format("plots/responseFitCBCB_cat%d_LOG.png",c));
    c1->SaveAs(TString::Format("plots/responseFitCBCB_cat%d_LOG.pdf",c));
    

    // saving as constant in the WS  
    w->defineSet(TString::Format("ResponseAddPdfParam_cat%d",c),RooArgSet(*w->var(TString::Format("ReducedMass_sig_sigma_cat%d",c)), 
									  *w->var(TString::Format("ReducedMass_sig_alphaCBpos_cat%d",c)),
									  *w->var(TString::Format("ReducedMass_sig_alphaCBneg_cat%d",c)),
									  *w->var(TString::Format("ReducedMass_sig_nCBpos_cat%d",c)),
									  *w->var(TString::Format("ReducedMass_sig_nCBneg_cat%d",c)),	   
									  *w->var(TString::Format("ReducedMass_sig_frac_cat%d",c)),  
									  *w->var(TString::Format("ReducedMass_sig_mean_cat%d",c))));
    SetConstantParams(w->set(TString::Format("ResponseAddPdfParam_cat%d",c)));

    w->Print();
  }
}

// Signal model: doubleCB. To describe the detector resolution, mgg/mggGen -1 used
void SigModelResponseReducedDoubleCBFit(RooWorkspace* w, Float_t mass, TString coupling) {
  
  int iMass = abs(mass);   
  
  // Dataset 
  RooDataSet* signal[NCAT];
  RooRealVar* mgg = w->var("mgg");     
  RooRealVar* massReduced = w->var("massReduced");
  
  // fit function
  RooDoubleCB* ResponseDoubleCB[NCAT];

  // to write down the results in a file
  ofstream *outTxtFile = new ofstream("HighMassGG_temp.rs",ios::out);
  *outTxtFile << "mgg[500,6000];" << endl;
  *outTxtFile << "mggGen[500,6000];" << endl;
  *outTxtFile << "weight[0, 1000];" << endl;
  *outTxtFile << "eventClass[-10,10];" << endl;
  *outTxtFile << endl;
  
  TCanvas* c1 = new TCanvas("c1", "c1", 1);
  c1->cd();
  
  for(int c = 0; c<NCAT; c++){

    // loading the dataset  
    signal[c] = (RooDataSet*) w->data(TString::Format("SigWeightReduced_cat%d",c));
    
    // DoubleCB 
    RooFormulaVar CBmean(TString::Format("ReducedMassCB_sig_mean_cat%d",c),"", "@0", *w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)));
    RooFormulaVar CBsigma(TString::Format("ReducedMassCB_sig_sigma_cat%d",c), "", "sqrt(@0*@0)", *w->var(TString::Format("ReducedMass_sig_sigma_cat%d",c)));
    RooFormulaVar CBalpha1(TString::Format("ReducedMassCB_sig_alpha1_cat%d",c),"", "@0", *w->var( TString::Format("ReducedMass_sig_alpha1_cat%d",c)));
    RooFormulaVar CBn1(TString::Format("ReducedMassCB_sig_n1cat%d",c),"", "@0", *w->var( TString::Format("ReducedMass_sig_n1_cat%d",c)));
    RooFormulaVar CBalpha2(TString::Format("ReducedMassCB_sig_alpha2_cat%d",c),"", "@0", *w->var( TString::Format("ReducedMass_sig_alpha2_cat%d",c)));
    RooFormulaVar CBn2(TString::Format("ReducedMassCB_sig_n2cat%d",c),"", "@0", *w->var( TString::Format("ReducedMass_sig_n2_cat%d",c)));
    ResponseDoubleCB[c] = new RooDoubleCB(TString::Format("ResponseDoubleCB_cat%d",c),TString::Format("ResponseDoubleCB_cat%d",c) , *massReduced, CBmean, CBsigma, CBalpha1, CBn1, CBalpha2, CBn2) ;
    w->import(*ResponseDoubleCB[c]);
   
    // Fit with ResponseDoubleCB      
    RooFitResult* fitresults = (RooFitResult* ) ResponseDoubleCB[c]->fitTo(*signal[c], SumW2Error(kFALSE), Range(0.95, 1.05), RooFit::Save(kTRUE));
    std::cout<<TString::Format("******************************** Signal Fit results doubleCB  mass %f cat %d***********************************", mass, c)<<std::endl;
    fitresults->Print("V");

    // saving fit parameters
    *outTxtFile << TString::Format("Sig_sigma_cat%d[",c)  << CBsigma.getVal()  << "];" << endl;
    *outTxtFile << TString::Format("Sig_alpha1_cat%d[",c) << CBalpha1.getVal() << "];" << endl;
    *outTxtFile << TString::Format("Sig_alpha2_cat%d[",c) << CBalpha2.getVal() << "];" << endl;
    *outTxtFile << TString::Format("Sig_n1_cat%d[",c)     << CBn1.getVal()     << "];" << endl;
    *outTxtFile << TString::Format("Sig_n2_cat%d[",c)     << CBn2.getVal()     << "];" << endl;
    *outTxtFile << TString::Format("Sig_mean_cat%d[",c)   << CBmean.getVal()   << "];" << endl;
    *outTxtFile << endl;

    // Plot
    RooPlot* plotG = massReduced->frame(Range(0.95, 1.05),Title("Mass Reduced"),Bins(60));
    plotG->GetXaxis()->SetTitle("#frac{m_{reco}}{m_{true}}");
    plotG->GetXaxis()->SetTitleFont(42);
    plotG->GetXaxis()->SetTitleSize(0.04);
    plotG->GetXaxis()->SetTitleOffset(1.40);
    plotG->GetYaxis()->SetTitle("Events/0.0024 units");
    plotG->GetYaxis()->SetTitleFont(42);
    plotG->GetYaxis()->SetTitleSize(0.04);

    TLegend* legmc = new TLegend(0.6, 0.58, 0.91, 0.91, "", "");
    legmc->SetTextSize(0.0286044);  
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->AddEntry(plotG->getObject(0),"Simulation","LP");    
    legmc->AddEntry(plotG->getObject(1),"DoubleCB ","L");
    
    TLatex* latex = new TLatex(0.21, 0.76, TString::Format("#splitline{m_{X}=%d GeV}{#splitline{}{Class %d}}",iMass,c));
    latex->SetTextSize(0.038);
    latex->SetTextAlign(11);
    latex->SetTextFont(42); 
    latex->SetNDC();

    signal[c]->plotOn(plotG);
    ResponseDoubleCB[c]->plotOn(plotG, LineColor(kBlue));
    plotG->Draw();
    latex->Draw("same");
    legmc->Draw("same");

    c1->SetLogy(0);
    c1->SaveAs(TString::Format("plots/responseDoubleCB_cat%d.png",c));
    c1->SetLogy();
    c1->SaveAs(TString::Format("plots/responseDoubleCB_cat%d_LOG.png",c));
    c1->SaveAs(TString::Format("plots/responseDoubleCB_cat%d_LOG.pdf",c));
    
    // saving as constant in the WS  
    w->defineSet(TString::Format("ResponseDoubleCBPdfParam_cat%d",c),RooArgSet(*w->var(TString::Format("ReducedMass_sig_sigma_cat%d",c)), 
									       *w->var(TString::Format("ReducedMass_sig_alpha1_cat%d",c)),
									       *w->var(TString::Format("ReducedMass_sig_alpha2_cat%d",c)),
									       *w->var(TString::Format("ReducedMass_sig_n1_cat%d",c)),
									       *w->var(TString::Format("ReducedMass_sig_n2_cat%d",c)),	   
									       *w->var(TString::Format("ReducedMass_sig_mean_cat%d",c))));
    SetConstantParams(w->set(TString::Format("ResponseDoubleCBPdfParam_cat%d",c)));

    w->Print();
  }
}

//-------------------------------------------------------------------------
// Signal model: BW only fit to the gen level mass - to check it is doable or not due to the PDFs
// Fit all categories together
void SigModelBWFit(RooWorkspace* w, Float_t mass, TString coupling) {

  int iMass = abs(mass);   
  
  // Variables
  RooRealVar* mggGen = w->var("mggGen");     

  // dataset
  RooDataSet* signal = (RooDataSet*) w->data(TString::Format("SigWeightGen"));

  // fit function
  RooFormulaVar meanBW("meanBWgen","","@0",*w->var("meanBW"));   
  RooFormulaVar sigmaBW("sigmaBWgen","","@0",*w->var("sigmaBW"));     
  RooBreitWigner *genMassBW = new RooBreitWigner("genMassBW","genMassBW",*mggGen,meanBW,sigmaBW);  
  w->import(*genMassBW);
  
  TCanvas* c1 = new TCanvas("c1", "c1", 1);
  c1->cd();

  std::cout<<TString::Format("******************************** gen level mass fit with BW, %f ***********************************",mass)<<std::endl;
  float fitMin = 1350.;
  float fitMax = 1650.;
  if (coupling=="001") {
    if (mass==750)  { fitMin = 740.;  fitMax = 760.; }
    if (mass==1500) { fitMin = 1490.; fitMax = 1510.; }
    if (mass==5000) { fitMin = 4980.; fitMax = 5020.; }
  } else if (coupling=="01") {
    if (mass==1500) { fitMin = 1200.;  fitMax = 1800.; }
    if (mass==3000) { fitMin = 2750.;  fitMax = 3250.; }
  } else if (coupling=="02") {
    if (mass==1500) { fitMin = 1000.;  fitMax = 2000.; }
    if (mass==3000) { fitMin = 2000.;  fitMax = 4000.; }
    if (mass==5000) { fitMin = 4000.;  fitMax = 6000.; }
  }
  RooFitResult* fitresults = (RooFitResult* ) genMassBW->fitTo(*signal, SumW2Error(kFALSE), Range(fitMin,fitMax), RooFit::Save(kFALSE));    

  // Plot
  RooPlot* plotG = mggGen->frame(Range(fitMin,fitMax),Title("Gen Level mgg"),Bins(60));
  plotG->GetXaxis()->SetTitle("m_{true}");
  plotG->GetXaxis()->SetTitleFont(42);
  plotG->GetXaxis()->SetTitleSize(0.04);
  plotG->GetXaxis()->SetTitleOffset(1.40);
  plotG->GetYaxis()->SetTitleFont(42);
  plotG->GetYaxis()->SetTitleSize(0.04);

  TLegend* legmc = new TLegend(0.6, 0.58, 0.91, 0.91, "", "");
  legmc->SetTextSize(0.0286044);  
  legmc->SetTextFont(42);
  legmc->SetBorderSize(0);
  legmc->SetFillStyle(0);
  legmc->AddEntry(plotG->getObject(0),"Simulation","LP");    
  legmc->AddEntry(plotG->getObject(1),"BW fit","L");
  
  TLatex* latex = new TLatex(0.21, 0.76, TString::Format("m_{X}=%d GeV",iMass));
  latex->SetTextSize(0.038);
  latex->SetTextAlign(11);
  latex->SetTextFont(42); 
  latex->SetNDC();

  signal->plotOn(plotG);
  genMassBW->plotOn(plotG, LineColor(kBlue));
  plotG->Draw();
  latex->Draw("same");
  legmc->Draw("same");

  c1->SetLogy(0);
  c1->SaveAs("plots/mggGenFitBW.png");
  c1->SetLogy();
  c1->SaveAs("plots/mggGenFitBW_LOG.png");
  
  w->defineSet("genMassBWPdfParam",RooArgSet(*w->var("meanBW"),*w->var("sigmaBW")));  
  SetConstantParams(w->set("genMassBWPdfParam"));

  // now put the dataset into a roodatahist and save
  TString myName= TString::Format("bSigWeightGen_mass%d",iMass)+TString("_")+TString(coupling);
  mggGen->setRange(fitMin,fitMax);
  mggGen->setBins(60);
  RooDataHist bSigWeightGen("bSigWeightGen","bSigWeightGen",RooArgList(*mggGen),*signal);
  bSigWeightGen.Print();
  w->import(bSigWeightGen,Rename(myName));
  w->Print();
  
  TFile myOutFile("outputHits.root","UPDATE");
  myOutFile.cd();
  bSigWeightGen.Write(myName);
}

//-------------------------------------------------------------------------

// Fit signal with model with CBCB convoluted with BW
void SigModelFitConvBW(RooWorkspace* w, Float_t mass) {

  Float_t MASS(mass);  

  // Dataset
  RooDataSet* sigToFit[NCAT];
  RooRealVar* mgg = w->var("mgg"); 

  // Fit to Signal 
  for (int c=0; c<NCAT; ++c) {
    cout << "---------- Category = " << c << endl;
    
    // CB - chiara: tornare alla versione di Livia che prende dal file quando avremo modello interpolato
    RooFormulaVar CBpos_mean(TString::Format("CBpos_mean_cat%d",c),"","@0",*w->var(TString::Format("ConvMass_sig_mean_cat%d",c)) );
    RooFormulaVar CBpos_sigma(TString::Format("CBpos_sigma_cat%d",c),"","@0",*w->var(TString::Format("Mass_sig_sigma_cat%d",c)) );
    RooFormulaVar CBpos_alphaCB(TString::Format("CBpos_alphaCB_cat%d",c),"","@0",*w->var(TString::Format("Mass_sig_alphaCBpos_cat%d",c)) );
    RooFormulaVar CBneg_alphaCB(TString::Format("CBneg_alphaCB_cat%d",c),"","@0",*w->var(TString::Format("Mass_sig_alphaCBneg_cat%d",c)) );
    RooFormulaVar CBpos_n(TString::Format("CBpos_n_cat%d",c),"","@0",*w->var(TString::Format("Mass_sig_nCBpos_cat%d",c)) );
    RooFormulaVar CBneg_n(TString::Format("CBneg_n_cat%d",c),"","@0",*w->var(TString::Format("Mass_sig_nCBneg_cat%d",c)) );
    RooFormulaVar CBpos_frac(TString::Format("CBpos_frac_cat%d",c),"","@0",*w->var(TString::Format("Mass_sig_frac_cat%d",c)) );

    RooCBShape ResCBpos(TString::Format("ResCBpos_cat%d",c),TString::Format("ResCBpos_cat%d",c) , *mgg, CBpos_mean, CBpos_sigma,CBpos_alphaCB, CBpos_n) ;
    RooCBShape ResCBneg(TString::Format("ResCBneg_cat%d",c),TString::Format("ResCBneg_cat%d",c) , *mgg, CBpos_mean, CBpos_sigma,CBneg_alphaCB, CBneg_n) ;
    mgg->setBins(5000, "cache");  
    
    RooAddPdf ResAddPdf(TString::Format("ResAddPdf_cat%d",c),TString::Format("ResAddPdf_cat%d",c) , RooArgList(ResCBneg, ResCBpos), CBpos_frac);

    // BW 
    RooFormulaVar meanBW(TString::Format("massBW_cat%d",c),"","@0",*w->var(TString::Format("meanBW_cat%d",c)) );  
    RooFormulaVar sigmaBW(TString::Format("widthBW_cat%d",c),"","@0",*w->var(TString::Format("sigmaBW_cat%d",c)) );
    RooBreitWigner SigModelBW(TString::Format("BW_cat%d",c),TString::Format("BW_cat%d",c), *mgg, meanBW, sigmaBW);
      
    // Convolution
    RooFFTConvPdf* ConvolutedRes_CB;
    ConvolutedRes_CB = new RooFFTConvPdf(TString::Format("mggSig_cat%d",c),TString::Format("mggSig_cat%d",c), *mgg,SigModelBW, ResAddPdf);
    w->import(*ConvolutedRes_CB);


    // Fit and Plot
    sigToFit[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c)); 
    
    float fitMin = 1250.;
    float fitMax = 1700.;
    if (mass==750)  {fitMin = 500.;   fitMax = 1200.; }
    if (mass==1500) {fitMin = 1250.;  fitMax = 1700.; }
    if (mass==5000) {fitMin = 4000.;  fitMax = 5500.; }

    RooFitResult* fitresults_CB = (RooFitResult* ) ConvolutedRes_CB->fitTo(*sigToFit[c], SumW2Error(kFALSE), Range(fitMin,fitMax), RooFit::Save(kTRUE));
    fitresults_CB->Print("V");
    
    RooPlot* plotOnlyResPdf = w->var("mgg")->frame(Range(fitMin,fitMax),Bins(100));
    sigToFit[c]->plotOn(plotOnlyResPdf, LineColor(kRed), LineStyle(kDashed));
    double max = plotOnlyResPdf->GetMaximum();
    
    RooPlot* plotPhotonsMassAll = w->var("mgg")->frame(Range(fitMin,fitMax),Bins(100));
    sigToFit[c]->plotOn(plotPhotonsMassAll);
    SigModelBW.plotOn(plotPhotonsMassAll, LineColor(kGreen), LineStyle(kDashed));
    ResAddPdf.plotOn(plotPhotonsMassAll, LineColor(kRed), LineStyle(kDashed));
    ConvolutedRes_CB->plotOn(plotPhotonsMassAll, LineColor(kBlue));
    
    TCanvas* c1 = new TCanvas("c1","PhotonsMass",0,0,800,800);
    c1->cd(1);
    plotPhotonsMassAll->Draw();  
    plotPhotonsMassAll->GetYaxis()->SetRangeUser(0.01, max*1.2);
    plotPhotonsMassAll->GetXaxis()->SetRangeUser(210, 290);
    
    TLatex *lat  = new TLatex(0.55,0.9,TString::Format("Cat: %d", c));  
    lat->SetTextSize(0.038);
    lat->SetTextAlign(11);
    lat->SetTextFont(42); 
    lat->SetNDC();
    
    TLegend *legmc = new TLegend(0.55, 0.6, 0.87, 0.88, "");
    legmc->AddEntry(plotPhotonsMassAll->getObject(0),"Simulation","LPE");
    legmc->AddEntry(plotPhotonsMassAll->getObject(1),"BW","L");
    legmc->AddEntry(plotPhotonsMassAll->getObject(2)," CB + CB ","L");
    legmc->AddEntry(plotPhotonsMassAll->getObject(3),"BW #otimes Resolution","L");
    legmc->SetTextSize(0.0206044);
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();
    lat->Draw("same");
    
    int massI(mass);
    c1->SetLogy();
    plotPhotonsMassAll->GetXaxis()->SetTitle("m_{#gamma #gamma}[GeV]");
    
    c1->SetLogy(0);
    c1->SaveAs("plots/signalCBCconvBW"+TString::Format(("_M%d_cat%d_free.png"),massI, c));
    
    c1->SetLogy();
    plotPhotonsMassAll->GetYaxis()->SetRangeUser(0.01,max*10. );
    plotPhotonsMassAll->GetXaxis()->SetRangeUser(210, 290);
    c1->SaveAs("plots/signalCBCconvBW"+TString::Format(("_M%d_cat%d_LOG_free.png"),massI,c));
    
    
    // IMPORTANT: fix all pdf parameters to constant
    w->defineSet(TString::Format("ConvolutedPdfParam_cat%d",c),RooArgSet( *w->var(TString::Format("Mass_sig_sigma_cat%d",c)), 
									  *w->var(TString::Format("Mass_sig_alphaCBpos_cat%d",c)),
									  *w->var(TString::Format("Mass_sig_alphaCBneg_cat%d",c)),
									  *w->var(TString::Format("Mass_sig_nCBpos_cat%d",c)),
									  *w->var(TString::Format("Mass_sig_nCBneg_cat%d",c)),	   
									  *w->var(TString::Format("Mass_sig_frac_cat%d",c)),  
									  *w->var(TString::Format("ConvMass_sig_mean_cat%d",c)),
									  *w->var(TString::Format("meanBW_cat%d",c)),
									  *w->var(TString::Format("sigmaBW_cat%d",c))));

    SetConstantParams(w->set(TString::Format("ConvolutedPdfParam_cat%d",c)));
    
    w->Print("V");
  }
}

// Fit signal with model with DoubleCB convoluted with BW
void SigModelFitDoubleCBConvBW(RooWorkspace* w, Float_t mass) {

  Float_t MASS(mass);  
  
  // Dataset
  RooDataSet* sigToFit[NCAT];
  RooRealVar* mgg = w->var("mgg"); 

  // Fit to Signal 
  for (int c=0; c<NCAT; ++c) {
    cout << "---------- Category = " << c << endl;
    
    // DoubleCB
    RooFormulaVar CBmean(TString::Format("CBmean_cat%d",c),"","@0",*w->var(TString::Format("ConvSig_mean_cat%d",c)) );
    RooFormulaVar CBsigma(TString::Format("CBsigma_cat%d",c),"","@0",*w->var(TString::Format("Sig_sigma_cat%d",c)) );
    RooFormulaVar CBalpha1(TString::Format("CBalpha1_cat%d",c),"","@0",*w->var(TString::Format("Sig_alpha1_cat%d",c)) );
    RooFormulaVar CBn1(TString::Format("CBn1_cat%d",c),"","@0",*w->var(TString::Format("Sig_n1_cat%d",c)) );
    RooFormulaVar CBalpha2(TString::Format("CBalpha2_cat%d",c),"","@0",*w->var(TString::Format("Sig_alpha2_cat%d",c)) );
    RooFormulaVar CBn2(TString::Format("CBn2_cat%d",c),"","@0",*w->var(TString::Format("Sig_n2_cat%d",c)) );
    RooDoubleCB ResponseDoubleCB(TString::Format("ResponseDoubleCB_cat%d",c),TString::Format("ResponseDoubleCB_cat%d",c) , *mgg, CBmean, CBsigma, CBalpha1, CBn1, CBalpha2, CBn2);
    mgg->setBins(5000, "cache");  

    // BW 
    RooFormulaVar meanBW(TString::Format("massBW_cat%d",c),"","@0",*w->var(TString::Format("meanBW_cat%d",c)) );  
    RooFormulaVar sigmaBW(TString::Format("widthBW_cat%d",c),"","@0",*w->var(TString::Format("sigmaBW_cat%d",c)) );
    RooBreitWigner SigModelBW(TString::Format("BW_cat%d",c),TString::Format("BW_cat%d",c), *mgg, meanBW, sigmaBW);
      
    // Convolution
    RooFFTConvPdf* ConvolutedRes_CB;
    ConvolutedRes_CB = new RooFFTConvPdf(TString::Format("mggSig_cat%d",c),TString::Format("mggSig_cat%d",c), *mgg,SigModelBW, ResponseDoubleCB);
    w->import(*ConvolutedRes_CB);


    // Fit and Plot
    sigToFit[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c)); 
    
    float fitMin = 1250.;
    float fitMax = 1700.;
    if (mass==750)  {fitMin = 500.;   fitMax = 1200.; }
    if (mass==1500) {fitMin = 1250.;  fitMax = 1700.; }
    if (mass==5000) {fitMin = 4000.;  fitMax = 5500.; }

    RooFitResult* fitresults_CB = (RooFitResult* ) ConvolutedRes_CB->fitTo(*sigToFit[c], SumW2Error(kFALSE), Range(fitMin,fitMax), RooFit::Save(kTRUE));
    fitresults_CB->Print("V");
    
    RooPlot* plotOnlyResPdf = w->var("mgg")->frame(Range(fitMin,fitMax),Bins(100));
    sigToFit[c]->plotOn(plotOnlyResPdf, LineColor(kRed), LineStyle(kDashed));
    double max = plotOnlyResPdf->GetMaximum();
    
    RooPlot* plotPhotonsMassAll = w->var("mgg")->frame(Range(fitMin,fitMax),Bins(100));
    plotPhotonsMassAll->SetTitle(TString::Format("DoubleCB conv BW, cat%d",c));
    sigToFit[c]->plotOn(plotPhotonsMassAll);
    SigModelBW.plotOn(plotPhotonsMassAll, LineColor(kGreen), LineStyle(kDashed));
    ResponseDoubleCB.plotOn(plotPhotonsMassAll, LineColor(kRed), LineStyle(kDashed));
    ConvolutedRes_CB->plotOn(plotPhotonsMassAll, LineColor(kBlue));
    
    TCanvas* c1 = new TCanvas("c1","PhotonsMass",0,0,800,800);
    c1->cd(1);
    plotPhotonsMassAll->Draw();  
    plotPhotonsMassAll->GetYaxis()->SetRangeUser(0.01, max*1.2);
    plotPhotonsMassAll->GetXaxis()->SetRangeUser(210, 290);
    
    TLegend *legmc = new TLegend(0.6, 0.58, 0.091, 0.91, "");
    legmc->AddEntry(plotPhotonsMassAll->getObject(0),"Simulation","LPE");
    legmc->AddEntry(plotPhotonsMassAll->getObject(1),"BW","L");
    legmc->AddEntry(plotPhotonsMassAll->getObject(2),"DoubleCB","L");
    legmc->AddEntry(plotPhotonsMassAll->getObject(3),"BW #otimes Resolution","L");
    legmc->SetTextSize(0.0206044);
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();
    
    int massI(mass);
    c1->SetLogy();
    plotPhotonsMassAll->GetXaxis()->SetTitle("m_{#gamma #gamma}[GeV]");
    
    c1->SetLogy(0);
    c1->SaveAs("plots/signalCBCconvBW"+TString::Format(("_M%d_cat%d_free.png"),massI, c));
    
    c1->SetLogy();
    plotPhotonsMassAll->GetYaxis()->SetRangeUser(0.01,max*10. );
    plotPhotonsMassAll->GetXaxis()->SetRangeUser(210, 290);
    c1->SaveAs("plots/signalCBCconvBW"+TString::Format(("_M%d_cat%d_LOG_free.png"),massI,c));
    
    
    // IMPORTANT: fix all pdf parameters to constant
    w->defineSet(TString::Format("ConvolutedPdfParam_cat%d",c),RooArgSet( *w->var(TString::Format("Sig_sigma_cat%d",c)), 
									  *w->var(TString::Format("Sig_alpha1_cat%d",c)),
									  *w->var(TString::Format("Sig_alpha2_cat%d",c)),
									  *w->var(TString::Format("Sig_n1_cat%d",c)),
									  *w->var(TString::Format("Sig_n2_cat%d",c)),	   
									  *w->var(TString::Format("ConvSig_mean_cat%d",c)),
									  *w->var(TString::Format("sigmaBW_cat%d",c)),
									  *w->var(TString::Format("meanBW_cat%d",c))));
    SetConstantParams(w->set(TString::Format("ConvolutedPdfParam_cat%d",c)));
    
    w->Print("V");
  }
}

// Fit signal using toys from fitted doubleCB and nominal BW
void SigModelFromToys(RooWorkspace* w, Float_t mass) {

  Float_t MASS(mass);  
  
  // Original dataset and variable
  RooDataSet* sigToFit[NCAT];
  RooRealVar* mgg = w->var("mgg"); 
  
  // Toy dataset
  RooDataSet* pseudoData[NCAT];
  RooRealVar* mggGenCB = new RooRealVar("mggGenCB", "M(gg)", -500., 500., "GeV");
  
  // Fit to Signal 
  for (int c=0; c<NCAT; ++c) {
    cout << "---------- Category = " << c << endl;
    
    // DoubleCB - must be centred at 1 with all parameters fixed to the result of the fit for this mass at k=0.01 (or no resonance)
    RooFormulaVar CBmean(TString::Format("CBmean_cat%d",c),"","@0",*w->var(TString::Format("Sig_mean_cat%d",c)) );
    RooFormulaVar CBsigma(TString::Format("CBsigma_cat%d",c),"","@0",*w->var(TString::Format("Sig_sigma_cat%d",c)) );
    RooFormulaVar CBalpha1(TString::Format("CBalpha1_cat%d",c),"","@0",*w->var(TString::Format("Sig_alpha1_cat%d",c)) );
    RooFormulaVar CBn1(TString::Format("CBn1_cat%d",c),"","@0",*w->var(TString::Format("Sig_n1_cat%d",c)) );
    RooFormulaVar CBalpha2(TString::Format("CBalpha2_cat%d",c),"","@0",*w->var(TString::Format("Sig_alpha2_cat%d",c)) );
    RooFormulaVar CBn2(TString::Format("CBn2_cat%d",c),"","@0",*w->var(TString::Format("Sig_n2_cat%d",c)) );
    RooDoubleCB ResponseDoubleCB(TString::Format("ResponseDoubleCB_cat%d",c),TString::Format("ResponseDoubleCB_cat%d",c) , *mggGenCB, CBmean, CBsigma, CBalpha1, CBn1, CBalpha2, CBn2);

    // BW - all parameters fixed to the expected values
    RooFormulaVar meanBW(TString::Format("massBW_cat%d",c),"","@0",*w->var(TString::Format("meanBW_cat%d",c)) );  
    RooFormulaVar sigmaBW(TString::Format("widthBW_cat%d",c),"","@0",*w->var(TString::Format("sigmaBW_cat%d",c)) );
    RooBreitWigner SigModelBW(TString::Format("BW_cat%d",c),TString::Format("BW_cat%d",c), *mgg, meanBW, sigmaBW);

    // Original dataset
    sigToFit[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c)); 
    int origEntries = sigToFit[c]->numEntries();

    // Pseudodata from our generation
    pseudoData[c] = new RooDataSet("pseudoData","pseudoData", RooArgSet(*mgg));
      
    // Generation from the two pdfs
    cout << "generating " << origEntries << " entries" << endl;
    for (int ii=0; ii<origEntries; ii++) {
      // for (int ii=0; ii<2000; ii++) {
      if (ii%100==0) cout << ii << endl;
      RooDataSet* dataCB = ResponseDoubleCB.generate(*mggGenCB,1);  
      RooDataSet* dataBW = SigModelBW.generate(*mgg,1);     
      RooArgSet setCB = *dataCB->get(0);
      RooArgSet setBW = *dataBW->get(0);
      RooRealVar* varCB = (RooRealVar*)setCB.find("mggGenCB");      
      RooRealVar* varBW = (RooRealVar*)setBW.find("mgg");      
      float mggCB = varCB->getVal();
      float mggBW = varBW->getVal();

      *mgg = mggCB*mggBW;  
      pseudoData[c]->add(RooArgSet(*mgg));  
    }
    pseudoData[c]->Print("V");
    w->import(*pseudoData[c],Rename(TString::Format("SigPseudodata_cat%d",c)));

    // Now build a new CB centred in 0 and with free parameters
    RooFormulaVar CBmeanF(TString::Format("CBmeanF_cat%d",c),"","@0",*w->var(TString::Format("SigF_mean_cat%d",c)) );
    RooFormulaVar CBsigmaF(TString::Format("CBsigmaF_cat%d",c),"","@0",*w->var(TString::Format("SigF_sigma_cat%d",c)) );
    RooFormulaVar CBalpha1F(TString::Format("CBalpha1F_cat%d",c),"","@0",*w->var(TString::Format("SigF_alpha1_cat%d",c)) );
    RooFormulaVar CBn1F(TString::Format("CBn1F_cat%d",c),"","@0",*w->var(TString::Format("SigF_n1_cat%d",c)) );
    RooFormulaVar CBalpha2F(TString::Format("CBalpha2F_cat%d",c),"","@0",*w->var(TString::Format("SigF_alpha2_cat%d",c)) );
    RooFormulaVar CBn2F(TString::Format("CBn2F_cat%d",c),"","@0",*w->var(TString::Format("SigF_n2_cat%d",c)) );
    RooDoubleCB ResponseDoubleCBF(TString::Format("ResponseDoubleCBF_cat%d",c),TString::Format("ResponseDoubleCBF_cat%d",c) , *mgg, CBmeanF, CBsigmaF, CBalpha1F, CBn1F, CBalpha2F, CBn2F);    
    mgg->setBins(5000, "cache");  

    // And convolve with the nominal BW - everything fixed
    RooFFTConvPdf* ConvolutedRes;
    ConvolutedRes = new RooFFTConvPdf(TString::Format("mggSig_cat%d",c),TString::Format("mggSig_cat%d",c), *mgg,SigModelBW, ResponseDoubleCBF);
    w->import(*ConvolutedRes);
    



    // Fit and Plot    
    float fitMin = 1250.;
    float fitMax = 1700.;
    if (mass==750)  {fitMin = 500.;   fitMax = 1200.; }
    if (mass==1500) {fitMin = 1250.;  fitMax = 1700.; }
    if (mass==5000) {fitMin = 4000.;  fitMax = 5500.; }

    RooFitResult* fitresults_CB = (RooFitResult* ) ConvolutedRes->fitTo(*pseudoData[c], SumW2Error(kFALSE), Range(fitMin,fitMax), RooFit::Save(kTRUE));
    fitresults_CB->Print("V");

    RooPlot* plotOnlyResPdf = w->var("mgg")->frame(Range(fitMin,fitMax),Bins(100));
    pseudoData[c]->plotOn(plotOnlyResPdf, LineColor(kRed), LineStyle(kDashed));
    double max = plotOnlyResPdf->GetMaximum();

    RooPlot* plotPhotonsMassAll = w->var("mgg")->frame(Range(fitMin,fitMax),Bins(100));
    plotPhotonsMassAll->GetXaxis()->SetTitle("m_{#gamma #gamma}[GeV]");
    plotPhotonsMassAll->SetTitle(TString::Format("DoubleCB conv BW from toys, cat%d",c));
    pseudoData[c]->plotOn(plotPhotonsMassAll);
    SigModelBW.plotOn(plotPhotonsMassAll, LineColor(kGreen), LineStyle(kDashed));
    ResponseDoubleCBF.plotOn(plotPhotonsMassAll, LineColor(kRed), LineStyle(kDashed));
    ConvolutedRes->plotOn(plotPhotonsMassAll, LineColor(kBlue));

    TCanvas* c1 = new TCanvas("c1","PhotonsMass",0,0,800,800);
    c1->cd(1);
    int massI(mass);    
    c1->SetLogy();
    plotPhotonsMassAll->GetYaxis()->SetRangeUser(0.01, max*10.);
    plotPhotonsMassAll->Draw();

    TLegend *legmc = new TLegend(0.6, 0.58, 0.091, 0.91, "");
    legmc->AddEntry(plotPhotonsMassAll->getObject(0),"Pseudodata","LPE");
    legmc->AddEntry(plotPhotonsMassAll->getObject(1),"BW","L");
    legmc->AddEntry(plotPhotonsMassAll->getObject(2),"DoubleCB","L");
    legmc->AddEntry(plotPhotonsMassAll->getObject(3),"BW #otimes Resolution","L");
    legmc->SetTextSize(0.0206044);
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();
    c1->SaveAs("plots/signalCBCconvBW"+TString::Format(("_M%d_cat%d_LOG_fromPseudoData.png"),massI,c));

    if(mass==750) plotPhotonsMassAll = w->var("mgg")->frame(Range(600,900),Bins(100));
    plotPhotonsMassAll->GetXaxis()->SetTitle("m_{#gamma #gamma}[GeV]");
    plotPhotonsMassAll->SetTitle(TString::Format("DoubleCB conv BW from toys, cat%d",c));
    pseudoData[c]->plotOn(plotPhotonsMassAll);
    SigModelBW.plotOn(plotPhotonsMassAll, LineColor(kGreen), LineStyle(kDashed));
    ResponseDoubleCBF.plotOn(plotPhotonsMassAll, LineColor(kRed), LineStyle(kDashed));
    ConvolutedRes->plotOn(plotPhotonsMassAll, LineColor(kBlue));

    c1->SetLogy(0);
    plotPhotonsMassAll->GetYaxis()->SetRangeUser(0.01, max*1.2);
    plotPhotonsMassAll->Draw();  
    legmc->Draw();
    c1->SaveAs("plots/signalCBCconvBW"+TString::Format(("_M%d_cat%d_fromPseudoData.png"),massI, c));

    // check that the generation was ok
    RooPlot* plotDatasets = w->var("mgg")->frame(Range(fitMin,fitMax),Bins(100));
    if(mass==750) plotDatasets = w->var("mgg")->frame(Range(500,900),Bins(100));
    plotDatasets->SetTitle(TString::Format("Test, cat%d",c));
    pseudoData[c]->plotOn(plotDatasets,MarkerColor(kRed));
    sigToFit[c]->plotOn(plotDatasets);
    c1->cd(1);
    c1->SetLogy(0);
    plotDatasets->Draw();  
    c1->SaveAs("plots/compareDataset"+TString::Format(("_cat%d.png"),c));  
    c1->SetLogy();
    plotDatasets->Draw();
    c1->SaveAs("plots/compareDataset"+TString::Format(("_cat%d_LOG.png"),c));


    // IMPORTANT: fix all pdf parameters to constant
    w->defineSet(TString::Format("ConvolutedPdfParam_cat%d",c),RooArgSet( *w->var(TString::Format("SigF_sigma_cat%d",c)), 
									  *w->var(TString::Format("SigF_alpha1_cat%d",c)),
									  *w->var(TString::Format("SigF_alpha2_cat%d",c)),
									  *w->var(TString::Format("SigF_n1_cat%d",c)),
									  *w->var(TString::Format("SigF_n2_cat%d",c)),	   
									  *w->var(TString::Format("SigF_mean_cat%d",c)),
									  *w->var(TString::Format("sigmaBW_cat%d",c)),
									  *w->var(TString::Format("meanBW_cat%d",c))));
    
    SetConstantParams(w->set(TString::Format("ConvolutedPdfParam_cat%d",c)));
    w->Print("V");
  }
}

//------------------------------------------------------------------------------------------------------
// Fit signal using toys from fitted doubleCB and binned theory distribution
void SigModelFromBinnedToys(RooWorkspace* w, Float_t mass, TString coupling) {

  Float_t MASS(mass);  
  int iMass = (int)mass;
  
  // Original dataset 
  RooDataSet* sigToFitReco[NCAT];

  // RooRealVars from the WS
  RooRealVar* mgg    = w->var("mgg"); 
  RooRealVar* mggGen = w->var("mggGen");

  // RooRealVars to generate
  RooRealVar* mggGenCB = new RooRealVar("mggGenCB", "M(gg)", -500., 500., "GeV");  
  
  // Toy dataset to describe the convolution
  RooDataSet* pseudoData[NCAT];
  
  // RooDataHist with mgg at gen level
  TFile *myRDHistFile = new TFile("outputHits.root","READ");
  TString myName= TString::Format("bSigWeightGen_mass%d",iMass)+TString("_")+coupling;
  cout << "loading " << myName << endl;
  RooDataHist *bSigWeightGen = (RooDataHist*)myRDHistFile->Get(myName);  
  RooHistPdf hSigWeightGen("hSigWeightGen","hSigWeightGen",*mggGen,*bSigWeightGen,0) ;

  // Fit to Signal 
  for (int c=0; c<NCAT; ++c) {
    cout << "---------- Category = " << c << endl;

    // Original dataset
    sigToFitReco[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c)); 
    int origEntries = sigToFitReco[c]->numEntries();

    // Pseudodata from our generation
    pseudoData[c] = new RooDataSet("pseudoData","pseudoData", RooArgSet(*mgg));  
    
    // DoubleCB - must be centred at 1 with all parameters fixed to the result of the fit for this mass at k=0.01 (or no resonance)
    RooFormulaVar CBmean(TString::Format("CBmean_cat%d",c),"","@0",*w->var(TString::Format("Sig_mean_cat%d",c)) );
    RooFormulaVar CBsigma(TString::Format("CBsigma_cat%d",c),"","@0",*w->var(TString::Format("Sig_sigma_cat%d",c)) );
    RooFormulaVar CBalpha1(TString::Format("CBalpha1_cat%d",c),"","@0",*w->var(TString::Format("Sig_alpha1_cat%d",c)) );
    RooFormulaVar CBn1(TString::Format("CBn1_cat%d",c),"","@0",*w->var(TString::Format("Sig_n1_cat%d",c)) );
    RooFormulaVar CBalpha2(TString::Format("CBalpha2_cat%d",c),"","@0",*w->var(TString::Format("Sig_alpha2_cat%d",c)) );
    RooFormulaVar CBn2(TString::Format("CBn2_cat%d",c),"","@0",*w->var(TString::Format("Sig_n2_cat%d",c)) );
    RooDoubleCB ResponseDoubleCB(TString::Format("ResponseDoubleCB_cat%d",c),TString::Format("ResponseDoubleCB_cat%d",c) , *mggGenCB, CBmean, CBsigma, CBalpha1, CBn1, CBalpha2, CBn2);

    // Generation from the two pdfs
    cout << "generating " << origEntries << " entries" << endl;
    for (int ii=0; ii<origEntries; ii++) {
    // for (int ii=0; ii<5000; ii++) {
      if (ii%500==0) cout << ii << endl;
      
      RooDataSet* dataCB = ResponseDoubleCB.generate(*mggGenCB,1);  
      RooDataSet* dataBW = hSigWeightGen.generate(*mggGen,1,AutoBinned(false));
      RooArgSet setCB = *dataCB->get(0);
      RooArgSet setBW = *dataBW->get(0);
      RooRealVar* varCB = (RooRealVar*)setCB.find("mggGenCB");      
      RooRealVar* varBW = (RooRealVar*)setBW.find("mggGen");      
      float mggCB = varCB->getVal();
      float mggBW = varBW->getVal();
      *mgg = mggCB*mggBW;  
      pseudoData[c]->add(RooArgSet(*mgg));   
    }
    cout << "my pseudodata" << endl;
    pseudoData[c]->Print("V");
    w->import(*pseudoData[c],Rename(TString::Format("SigPseudodata_cat%d",c)));
  
    TCanvas* c1 = new TCanvas("c1","PhotonsMass",0,0,800,800);
    RooPlot* plotDatasets = w->var("mgg")->frame(Range(1250,1700),Bins(100));
    if(mass==750)  plotDatasets = w->var("mgg")->frame(Range(500,900),Bins(100));
    if(mass==1500) plotDatasets = w->var("mgg")->frame(Range(1250,1700),Bins(100));
    plotDatasets->SetTitle(TString::Format("Toys vs data, cat%d",c));
    pseudoData[c]->plotOn(plotDatasets,MarkerColor(kRed));
    sigToFitReco[c]->plotOn(plotDatasets);
    c1->cd(1);
    c1->SetLogy(0);
    plotDatasets->Draw();  
    c1->SaveAs("plots/compareDataset"+TString::Format(("_cat%d.png"),c));  
    c1->SetLogy();
    plotDatasets->Draw();
    c1->SaveAs("plots/compareDataset"+TString::Format(("_cat%d_LOG.png"),c));

    /*
    // Now build a new CB centred in 0 and with free parameters
    RooFormulaVar CBmeanF(TString::Format("CBmeanF_cat%d",c),"","@0",*w->var(TString::Format("SigF_mean_cat%d",c)) );
    RooFormulaVar CBsigmaF(TString::Format("CBsigmaF_cat%d",c),"","@0",*w->var(TString::Format("SigF_sigma_cat%d",c)) );
    RooFormulaVar CBalpha1F(TString::Format("CBalpha1F_cat%d",c),"","@0",*w->var(TString::Format("SigF_alpha1_cat%d",c)) );
    RooFormulaVar CBn1F(TString::Format("CBn1F_cat%d",c),"","@0",*w->var(TString::Format("SigF_n1_cat%d",c)) );
    RooFormulaVar CBalpha2F(TString::Format("CBalpha2F_cat%d",c),"","@0",*w->var(TString::Format("SigF_alpha2_cat%d",c)) );
    RooFormulaVar CBn2F(TString::Format("CBn2F_cat%d",c),"","@0",*w->var(TString::Format("SigF_n2_cat%d",c)) );
    RooDoubleCB ResponseDoubleCBF(TString::Format("ResponseDoubleCBF_cat%d",c),TString::Format("ResponseDoubleCBF_cat%d",c) , *mgg, CBmeanF, CBsigmaF, CBalpha1F, CBn1F, CBalpha2F, CBn2F);    
    mgg->setBins(5000, "cache");  

    // And convolve with the nominal BW - everything fixed
    //RooFormulaVar meanBW(TString::Format("massBW_cat%d",c),"","@0",*w->var(TString::Format("meanBW_cat%d",c)) );  
    //RooFormulaVar sigmaBW(TString::Format("widthBW_cat%d",c),"","@0",*w->var(TString::Format("sigmaBW_cat%d",c)) );
    //RooBreitWigner SigModelBW(TString::Format("BW_cat%d",c),TString::Format("BW_cat%d",c), *mgg, meanBW, sigmaBW);



    RooFFTConvPdf* ConvolutedRes;
    ConvolutedRes = new RooFFTConvPdf(TString::Format("mggSig_cat%d",c),TString::Format("mggSig_cat%d",c), *mgg,SigModelBW, ResponseDoubleCBF);
    w->import(*ConvolutedRes);
    
    // Fit and Plot    
    float fitMin = 1250.;
    float fitMax = 1700.;
    if (mass==750)  {fitMin = 500.;   fitMax = 1200.; }
    if (mass==1500) {fitMin = 1250.;  fitMax = 1700.; }
    if (mass==5000) {fitMin = 4000.;  fitMax = 5500.; }

    RooFitResult* fitresults_CB = (RooFitResult* ) ConvolutedRes->fitTo(*pseudoData[c], SumW2Error(kFALSE), Range(fitMin,fitMax), RooFit::Save(kTRUE));
    fitresults_CB->Print("V");

    RooPlot* plotOnlyResPdf = w->var("mgg")->frame(Range(fitMin,fitMax),Bins(100));
    pseudoData[c]->plotOn(plotOnlyResPdf, LineColor(kRed), LineStyle(kDashed));
    double max = plotOnlyResPdf->GetMaximum();
    

    RooPlot* plotPhotonsMassAll = w->var("mgg")->frame(Range(fitMin,fitMax),Bins(100));
    plotPhotonsMassAll->SetTitle(TString::Format("DoubleCB conv BW from toys, cat%d",c));
    pseudoData[c]->plotOn(plotPhotonsMassAll);
    SigModelBW.plotOn(plotPhotonsMassAll, LineColor(kGreen), LineStyle(kDashed));
    ResponseDoubleCBF.plotOn(plotPhotonsMassAll, LineColor(kRed), LineStyle(kDashed));
    ConvolutedRes->plotOn(plotPhotonsMassAll, LineColor(kBlue));
    
    TCanvas* c1 = new TCanvas("c1","PhotonsMass",0,0,800,800);
    c1->cd(1);
    plotPhotonsMassAll->Draw();  
    plotPhotonsMassAll->GetYaxis()->SetRangeUser(0.01, max*1.2);
    plotPhotonsMassAll->GetXaxis()->SetRangeUser(210, 290);
    
    TLegend *legmc = new TLegend(0.6, 0.58, 0.091, 0.91, "");
    legmc->AddEntry(plotPhotonsMassAll->getObject(0),"Pseudodata","LPE");
    legmc->AddEntry(plotPhotonsMassAll->getObject(1),"BW","L");
    legmc->AddEntry(plotPhotonsMassAll->getObject(2),"DoubleCB","L");
    legmc->AddEntry(plotPhotonsMassAll->getObject(3),"BW #otimes Resolution","L");
    legmc->SetTextSize(0.0206044);
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();
    
    int massI(mass);
    c1->SetLogy();
    plotPhotonsMassAll->GetXaxis()->SetTitle("m_{#gamma #gamma}[GeV]");

    c1->SetLogy(0);
    c1->SaveAs("plots/signalCBCconvBW"+TString::Format(("_M%d_cat%d_fromPseudoData.png"),massI, c));
    
    c1->SetLogy();
    plotPhotonsMassAll->GetYaxis()->SetRangeUser(0.01,max*10. );
    plotPhotonsMassAll->GetXaxis()->SetRangeUser(210, 290);
    c1->SaveAs("plots/signalCBCconvBW"+TString::Format(("_M%d_cat%d_LOG_fromPseudoData.png"),massI,c));


    // check that the generation was ok
    RooPlot* plotDatasets = w->var("mgg")->frame(Range(fitMin,fitMax),Bins(100));
    plotDatasets->SetTitle(TString::Format("Test, cat%d",c));
    pseudoData[c]->plotOn(plotDatasets,MarkerColor(kRed));
    sigToFitReco[c]->plotOn(plotDatasets);
    c1->cd(1);
    c1->SetLogy(0);
    plotDatasets->Draw();  
    c1->SaveAs("plots/compareDataset"+TString::Format(("_cat%d.png"),c));  
    c1->SetLogy();
    plotDatasets->Draw();
    c1->SaveAs("plots/compareDataset"+TString::Format(("_cat%d_LOG.png"),c));


    // IMPORTANT: fix all pdf parameters to constant
    w->defineSet(TString::Format("ConvolutedPdfParam_cat%d",c),RooArgSet( *w->var(TString::Format("SigF_sigma_cat%d",c)), 
									  *w->var(TString::Format("SigF_alpha1_cat%d",c)),
									  *w->var(TString::Format("SigF_alpha2_cat%d",c)),
									  *w->var(TString::Format("SigF_n1_cat%d",c)),
									  *w->var(TString::Format("SigF_n2_cat%d",c)),	   
									  *w->var(TString::Format("SigF_mean_cat%d",c)),
									  *w->var(TString::Format("sigmaBW_cat%d",c)),
									  *w->var(TString::Format("meanBW_cat%d",c))));
    
    SetConstantParams(w->set(TString::Format("ConvolutedPdfParam_cat%d",c)));
    w->Print("V");
    */

  }
}
//------------------------------------------------------------------------------------------------------


// Write signal pdfs and datasets into the workspace 
void MakeSigWS(RooWorkspace* w, const char* fileBaseName, Float_t mass, std::string coupling){
  
  TString wsDir = "workspaces/";
  Int_t ncat = NCAT;
  
  RooWorkspace *wAll = new RooWorkspace("w_all","w_all");  

  // Retrieve stuff
  w->Print("V");
  for (int c=0; c<ncat; ++c) {

    // convolution
    wAll->import(*w->pdf("mggSig"+TString::Format("_cat%d",c)));  // this is the convolution
    
    // dataset
    wAll->import(*w->data(TString::Format("SigWeight_cat%d",c)));            // real dataset
    // wAll->import(*w->data(TString::Format("SigPseudodata_cat%d",c)));     // dataset from toys
  }
  std::cout << "done with importing signal pdf and datasets" << std::endl;

  wAll->import(*w->var("massReduced"));
  wAll->import(*w->var("mggGen"));

  int imass = (int)mass;
  TString filename(wsDir+TString(fileBaseName)+TString::Format("_m%d",imass)+"_"+TString(coupling)+".root");
  TFile fileWs(filename,"RECREATE");
  fileWs.cd();
  wAll->writeToFile(filename);
  cout << "Write signal workspace in: " << filename << " file" << endl;
  
  return;
}

void runfits(const Float_t mass=1500, string coupling="001") {

  //******************************************************************//
  //  Steps:
  //     - create signal and background data sets 
  //     - make and fit signal and background  models 
  //     - write signal and background workspaces in root files
  //     - write data card
  //*******************************************************************//
  
  TString fileBaseName("HighMassGG");    
  TString fileBkgName("HighMassGG.inputbkg");
  HLFactory hlf("HLFactory", "HighMassGG.rs", false);
  RooWorkspace* w = hlf.GetWs();
 
  // import luminosity in the ws
  RooRealVar lumi("lumi","lumi",Lum);
  w->import(lumi); 
  
  // range for the variables
  w->var("mgg")->setMin(MINmass);
  w->var("mgg")->setMax(MAXmass);
  w->var("mggGen")->setMin(MINmass);
  w->var("mggGen")->setMax(MAXmass);
  w->Print("v");
  
  cout << endl; 
  cout << "Now add signal data" << endl;
  AddSigData(w, mass, coupling);   

  cout << endl; 
  if (wantResponse) {
    cout << "Now prepare signal model fit - resolution function" << endl;  
    // SigModelResponseDoubleCBFit(w, mass, coupling); 
    SigModelResponseReducedDoubleCBFit(w, mass, coupling); 
  }    

  cout << endl;
  if (wantGenLevel) {
    cout << "Now try BW only on gen level mgg" << endl;
    SigModelBWFit(w, mass, coupling);     
  }
  
  cout << endl;
  cout << "Now prepare signal model fit - resolution function x BW" << endl;  
  // SigModelFitDoubleCBConvBW(w, mass);

  cout << endl;
  cout << "Alternative: generate following nominal BW and fitted doubleCB and build another dataset" << endl;
  // SigModelFromToys(w, mass);  
  SigModelFromBinnedToys(w, mass,coupling); 

  cout << endl;
  cout << "Now preparing signal WS" << endl;
  // MakeSigWS(w, fileBaseName, mass, coupling);

  return;
}
