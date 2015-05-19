#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooCBShape.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooBreitWigner.h"
#include "RooFFTConvPdf.h"
#include "RooHistFunc.h"
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

using namespace RooFit;
using namespace RooStats;

// ============================================
// to be modified:
static const Int_t NCAT = 4;  
Int_t MINmass= 500;
Int_t MAXmass= 6000;
std::string filePOSTfix="";
double signalScaler=1.00;
Float_t Lum = 19500.0;    
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
 
  TPaveText* cmslabel = new TPaveText( x1, y1, x2, y2, "brNDC" );
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

  TPaveText* label_sqrt = new TPaveText(x1,y1,x2,y2, "brNDC");
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
  RooFormulaVar *massReduced_formula = new RooFormulaVar("massReduced_formula","","@0/@1 -1",RooArgList(*w->var("mgg"),*w->var("mggGen")));
  RooRealVar* massReduced = (RooRealVar*) sigWeighted.addColumn(*massReduced_formula);
  massReduced->SetName("massReduced");
  massReduced->SetTitle("massReduced");
  w->import(*massReduced);  
  massReduced->setRange(-0.5, 0.5);
  
  // common preselection cut on the reduced mass
  TString mainCut = TString::Format("massReduced>-0.5 && massReduced <0.5"); 

  
  // -------------------------
  // split in categories, wrt mgg
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
  // split in categories, wrt massReduced
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
      fitresults = (RooFitResult* ) ResponseDoubleCB[c]->fitTo(*signal[c], SumW2Error(kFALSE), Range(400, 1200), RooFit::Save(kTRUE));
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

    TLegend* legmc = new TLegend(0.6, 0.58, 0.91, 0.91, "", "brNDC");
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
      plotG = mgg->frame(Range(1250, 1650),Title("Mgg, response"),Bins(60));        
    } else if (mass==750) {
      plotG = mgg->frame(Range(710, 780),Title("Mgg, response"),Bins(60));          
    } else if (mass==5000) {
      plotG = mgg->frame(Range(4000, 5500),Title("Mgg, response"),Bins(60));        
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
      plotG = mgg->frame(Range(1400, 1550),Title("Mgg, response"),Bins(60));           
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
      fitresults = (RooFitResult* ) ResponseAdd[c]->fitTo(*signal[c], SumW2Error(kFALSE), Range(400, 1200), RooFit::Save(kTRUE));
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
    TLegend* legmc = new TLegend(0.6, 0.58, 0.91, 0.91, "", "brNDC");
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
      plotG = mgg->frame(Range(1250, 1650),Title("Mgg, response"),Bins(60));              
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
      plotG = mgg->frame(Range(1400, 1550),Title("Mgg, response"),Bins(60));           
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
    RooFitResult* fitresults = (RooFitResult* ) ResponseAdd[c]->fitTo(*signal[c], SumW2Error(kFALSE), Range(-0.05, 0.05), RooFit::Save(kTRUE));
    std::cout<<TString::Format("******************************** Signal Fit results CB+CB  mass %f cat %d***********************************", mass, c)<<std::endl;
    fitresults->Print("V");
   

    // Plot
    RooPlot* plotG = massReduced->frame(Range(-0.05, 0.05),Title("Mass Reduced"),Bins(60));
    signal[c]->plotOn(plotG);

    ResponseAdd[c]->plotOn(plotG, LineColor(kBlue));
    ResponseAdd[c]->plotOn(plotG,Components(TString::Format("ResponseCBneg_cat%d",c)), LineColor(kGreen), LineStyle(kDotted));
    ResponseAdd[c]->plotOn(plotG,Components(TString::Format("ResponseCBpos_cat%d",c)), LineColor(kRed), LineStyle(kDashed));
  
    plotG->GetXaxis()->SetTitle("#frac{m_{reco}-m_{true}}{m_{true}}");
    plotG->GetXaxis()->SetTitleFont(42);
    plotG->GetXaxis()->SetTitleSize(0.04);
    plotG->GetXaxis()->SetTitleOffset(1.40);
    plotG->GetYaxis()->SetTitle("Events/0.0024 units");
    plotG->GetYaxis()->SetTitleFont(42);
    plotG->GetYaxis()->SetTitleSize(0.04);

    TLegend* legmc = new TLegend(0.6, 0.58, 0.91, 0.91, "", "brNDC");
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
    RooFitResult* fitresults = (RooFitResult* ) ResponseDoubleCB[c]->fitTo(*signal[c], SumW2Error(kFALSE), Range(-0.05, 0.05), RooFit::Save(kTRUE));
    std::cout<<TString::Format("******************************** Signal Fit results doubleCB  mass %f cat %d***********************************", mass, c)<<std::endl;
    fitresults->Print("V");
   
    // Plot
    RooPlot* plotG = massReduced->frame(Range(-0.05, 0.05),Title("Mass Reduced"),Bins(60));
    plotG->GetXaxis()->SetTitle("#frac{m_{reco}-m_{true}}{m_{true}}");
    plotG->GetXaxis()->SetTitleFont(42);
    plotG->GetXaxis()->SetTitleSize(0.04);
    plotG->GetXaxis()->SetTitleOffset(1.40);
    plotG->GetYaxis()->SetTitle("Events/0.0024 units");
    plotG->GetYaxis()->SetTitleFont(42);
    plotG->GetYaxis()->SetTitleSize(0.04);

    TLegend* legmc = new TLegend(0.6, 0.58, 0.91, 0.91, "", "brNDC");
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
void SigModelBWFit(RooWorkspace* w, Float_t mass, TString coupling) {

  int iMass = abs(mass);   

  // Files
  TString inDir = "data/newSelection/mergedFinal/";
  TChain* sigTree = new TChain();
  sigTree->Add(inDir+TString(Form("RSGravToGG_kMpl-"))+coupling+TString(Form("_M-%d.root/DiPhotonTree", iMass)));
  sigTree->SetTitle("sigTree");
  sigTree->SetName("sigTree");

  // Variables
  RooArgSet* ntplVars = defineVariables();
  
  // common preselection cut
  TString mainCut = TString::Format("mgg>=500 && mgg<=6000 && mggGen>=500 && mggGen<=6000");   
  // RooDataSet sigWeightedGen("sigWeightedGen","dataset",sigTree,*ntplVars,mainCut,"weight");
  RooDataSet sigWeightedGen("sigWeightedGen","dataset",sigTree,*ntplVars,mainCut);
  RooRealVar* mggGen = w->var("mggGen");     

  // fit functions
  RooDataSet* signal[NCAT];
  RooBreitWigner *genMassBW[NCAT]; 
  
  TCanvas* c1 = new TCanvas("c1", "c1", 1);
  c1->cd();
  TPaveText* label_cms  = get_labelCMS(0, "2015", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
  
  for(int c = 0; c<NCAT; c++){
    
    // splitting in categories
    if (c==0) signal[c] = (RooDataSet*) sigWeightedGen.reduce(*w->var("mggGen"),mainCut+TString::Format("&& eventClass==0"));
    if (c==1) signal[c] = (RooDataSet*) sigWeightedGen.reduce(*w->var("mggGen"),mainCut+TString::Format("&& eventClass==1"));
    if (c==2) signal[c] = (RooDataSet*) sigWeightedGen.reduce(*w->var("mggGen"),mainCut+TString::Format("&& eventClass==2"));
    if (c==3) signal[c] = (RooDataSet*) sigWeightedGen.reduce(*w->var("mggGen"),mainCut+TString::Format("&& eventClass==3"));
    w->import(*signal[c],Rename(TString::Format("SigWeightGen_cat%d",c))); 
    
    // BW
    RooFormulaVar meanBW(TString::Format("meanBWgen_cat%d",c),"","@0",*w->var(TString::Format("meanBW_cat%d",c)));   
    RooFormulaVar sigmaBW(TString::Format("sigmaBWgen_cat%d",c),"","@0",*w->var(TString::Format("sigmaBW_cat%d",c)));     
    genMassBW[c] = new RooBreitWigner(TString::Format("genMassBW_cat%d",c),TString::Format("genMassBW_cat%d",c),*mggGen,meanBW,sigmaBW);  

    w->import(*genMassBW[c]);
   
    // Fit with this BW
    // chiara
    RooFitResult* fitresults = 0;
    if (coupling=="001") {
      if (mass==750)  fitresults = (RooFitResult* ) genMassBW[c]->fitTo(*signal[c], SumW2Error(kFALSE), Range(740, 760), RooFit::Save(kFALSE));    
      if (mass==1500) fitresults = (RooFitResult* ) genMassBW[c]->fitTo(*signal[c], SumW2Error(kFALSE), Range(1490, 1510), RooFit::Save(kFALSE));
      if (mass==5000) fitresults = (RooFitResult* ) genMassBW[c]->fitTo(*signal[c], SumW2Error(kFALSE), Range(4980, 5020), RooFit::Save(kFALSE));
    } else if (coupling=="01") {
      if (mass==1500) fitresults = (RooFitResult* ) genMassBW[c]->fitTo(*signal[c], SumW2Error(kFALSE), Range(1400, 1600), RooFit::Save(kFALSE));
      if (mass==3000) fitresults = (RooFitResult* ) genMassBW[c]->fitTo(*signal[c], SumW2Error(kFALSE), Range(2750, 3250), RooFit::Save(kFALSE));
    } else if (coupling=="02") {
      if (mass==1500) fitresults = (RooFitResult* ) genMassBW[c]->fitTo(*signal[c], SumW2Error(kFALSE), Range(1000, 2000), RooFit::Save(kFALSE));
      if (mass==3000) fitresults = (RooFitResult* ) genMassBW[c]->fitTo(*signal[c], SumW2Error(kFALSE), Range(2000, 4000), RooFit::Save(kFALSE));
      if (mass==5000) fitresults = (RooFitResult* ) genMassBW[c]->fitTo(*signal[c], SumW2Error(kFALSE), Range(4000, 6000), RooFit::Save(kFALSE));
    }
    std::cout<<TString::Format("******************************** gen level mass fit with BW, %f cat %d***********************************", mass, c)<<std::endl;
    fitresults->Print("V");
   
    // Plot
    // chiara
    RooPlot* plotG = mggGen->frame(Range(740,760),Title("Gen Level mgg"),Bins(60));
    if (coupling=="001") {
      if (mass==750)  plotG = mggGen->frame(Range(740,760),Title("Gen Level mgg"),Bins(60));
      if (mass==1500) plotG = mggGen->frame(Range(1490,1510),Title("Gen Level mgg"),Bins(60));
      if (mass==5000) plotG = mggGen->frame(Range(4980,5020),Title("Gen Level mgg"),Bins(60));
    } else if (coupling=="01") {
      if (mass==1500) plotG = mggGen->frame(Range(1400,1600),Title("Gen Level mgg"),Bins(60));
      if (mass==3000) plotG = mggGen->frame(Range(2750,3250),Title("Gen Level mgg"),Bins(60));
    } else if (coupling=="02") {
      if (mass==1500) plotG = mggGen->frame(Range(1000,2000),Title("Gen Level mgg"),Bins(60));
      if (mass==3000) plotG = mggGen->frame(Range(2000,4000),Title("Gen Level mgg"),Bins(60));
      if (mass==5000) plotG = mggGen->frame(Range(4000,6000),Title("Gen Level mgg"),Bins(60));
    }
    signal[c]->plotOn(plotG);

    genMassBW[c]->plotOn(plotG, LineColor(kBlue));
  
    plotG->GetXaxis()->SetTitle("m_{true}");
    plotG->GetXaxis()->SetTitleFont(42);
    plotG->GetXaxis()->SetTitleSize(0.04);
    plotG->GetXaxis()->SetTitleOffset(1.40);
    plotG->GetYaxis()->SetTitleFont(42);
    plotG->GetYaxis()->SetTitleSize(0.04);

    TLegend* legmc = new TLegend(0.6, 0.58, 0.91, 0.91, "", "brNDC");
    legmc->SetTextSize(0.0286044);  
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->AddEntry(plotG->getObject(0),"Simulation","LP");    
    legmc->AddEntry(plotG->getObject(1),"BW fit","L");
    
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
    c1->SaveAs(TString::Format("plots/mggGenFitBW_cat%d.png",c));
    c1->SetLogy();
    c1->SaveAs(TString::Format("plots/mggGenFitBW_cat%d_LOG.png",c));
	       
    w->defineSet(TString::Format("genMassBWPdfParam_cat%d",c),RooArgSet(*w->var(TString::Format("meanBW_cat%d",c)), 
									*w->var(TString::Format("sigmaBW_cat%d",c))));
    SetConstantParams(w->set(TString::Format("genMassBWPdfParam_cat%d",c)));

    w->Print();
  }
}
//-------------------------------------------------------------------------

// Fit signal with model with CB convoluted with BW
void SigModelFitConvBW(RooWorkspace* w, Float_t mass, Double_t width) {

  Float_t MASS(mass);  
  Float_t minMassFit(mass*0.8);
  Float_t maxMassFit(mass*1.2); 

  // chiara
  RooDataSet* sigToFit[NCAT];
  RooRealVar* mgg = w->var("mgg"); 

  // systematics removed for the moment

  TPaveText* label_cms  = get_labelCMS(0, "2015", true);
  TPaveText* label_sqrt = get_labelSqrt(0);

  // chiara: to be created
  // TFile* f = new TFile("sigShapeCorrections.root", "READ"); 

  // Fit to Signal 
  for (int c=0; c<NCAT; ++c) {
    cout << "---------- Category = " << c << endl;
    
    // chiara: to be uncommented when we have the parametric model
    // get sigma from TF1:     
    // TF1* fcn = (TF1*)f->Get(TString::Format("f%d",c));   
    // Float_t massF = (Float_t) (*w->var("MG")).getVal();   
    // Float_t sigmaCorr = fcn->Eval(massF);
    // if(massF==1500)sigmaCorr=1;   
    // Float_t sigmaCorr=1;   
    // RooRealVar rooSigmaCorr (TString::Format("rooSigmaCorr_cat%d",c), TString::Format("rooSigmaCorr_cat%d",c), sigmaCorr, "");
    // rooSigmaCorr.setConstant(); 
    // w->import(rooSigmaCorr);

    // chiara: to be uncommented for syst.
    // ( *w->var(TString::Format("mShift_cat%d",c))).setConstant();
    // ( *w->var(TString::Format("mSmear_cat%d",c))).setConstant(); 

    // CB - chiara: tornare alla versione di Livia che prende dal file quando avremo modello interpolato
    RooFormulaVar CBpos_mean(TString::Format("CBpos_mean_cat%d",c),"","@0",*w->var(TString::Format("ConvMass_sig_mean_cat%d",c)) );
    RooFormulaVar CBpos_sigma(TString::Format("CBpos_sigma_cat%d",c),"","@0",*w->var(TString::Format("Mass_sig_sigma_cat%d",c)) );
    ////RooFormulaVar CBpos_sigma(TString::Format("CBpos_sigma_cat%d",c),"","(sqrt(@0*@0*@2*@2)*@1)",RooArgList(*w->var(TString::Format("Mass_sig_sigma_cat%d",c)),*w->var("MG"),*w->var(TString::Format("rooSigmaCorr_cat%d",c)) ) );   
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

    // chiara: provo libero
    RooFormulaVar meanBW(TString::Format("massBW_cat%d",c),"","@0",*w->var(TString::Format("meanBW_cat%d",c)) );  
    RooFormulaVar sigmaBW(TString::Format("widthBW_cat%d",c),"","@0",*w->var(TString::Format("sigmaBW_cat%d",c)) );
    RooBreitWigner SigModelBW(TString::Format("BW_cat%d",c),TString::Format("BW_cat%d",c), *mgg, meanBW, sigmaBW);
   /*
    // chiara era cosi
    // RooFormulaVar meanBW(TString::Format("meanBW_cat%d",c),"","@0",*w->var("MG"));    
    // RooRealVar sigmaBW_var(TString::Format("sigmaBW_var_cat%d",c), TString::Format("sigmaBW_var_cat%d",c), width);
    // sigmaBW_var.setConstant();  
    // cout << "import width = " << width << endl;
    // w->import(sigmaBW_var);
    // RooFormulaVar* sigmaBW;
    // if(width<1) sigmaBW = new RooFormulaVar(TString::Format("sigmaBW_cat%d",c),"","@0",*w->var(TString::Format("sigmaBW_var_cat%d",c))); 

    // chiara: se voglio allargare alla Livia
    // else if(width==2)  sigmaBW = new RooFormulaVar(TString::Format("sigmaBW_cat%d",c),"","@0*0.02",*w->var("MG"));   
    // else if(width==5)  sigmaBW = new RooFormulaVar(TString::Format("sigmaBW_cat%d",c),"","@0*0.05",*w->var("MG"));   
    // else if(width==7)  sigmaBW = new RooFormulaVar(TString::Format("sigmaBW_cat%d",c),"","@0*0.07",*w->var("MG"));   
    // else if(width==10) sigmaBW = new RooFormulaVar(TString::Format("sigmaBW_cat%d",c),"","@0*0.10",*w->var("MG"));   
        
    // BW with mean and sigma as above
    // RooBreitWigner SigModelBW(TString::Format("BW_cat%d",c),TString::Format("BW_cat%d",c), *mgg, meanBW, *sigmaBW);
    */
      
    // Convolution
    RooFFTConvPdf* ConvolutedRes_CB;
    ConvolutedRes_CB = new RooFFTConvPdf(TString::Format("mggSig_cat%d",c),TString::Format("mggSig_cat%d",c), *mgg,SigModelBW, ResAddPdf);
    w->import(*ConvolutedRes_CB);

    // chiara: serve?
    // RooHistFunc* rooFunc_norm = getRooHistFunc(c,w->var("MG"));
    // w->import(*rooFunc_norm);
    // std::cout<<"SIG NORM ----->"<<rooFunc_norm->getVal(*w->var("MG"));

    // make plot for some cases
    if(1){  // chiara      
    // if(width < 1){       

      sigToFit[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c)); 
      //sigToFit[c].Print("v");
      
      RooFitResult* fitresults_CB = 0;
      if (mass==1500)
	fitresults_CB = (RooFitResult* ) ConvolutedRes_CB->fitTo(*sigToFit[c], SumW2Error(kFALSE), Range(1250, 1700), RooFit::Save(kTRUE));
      else if (mass==750)
	fitresults_CB = (RooFitResult* ) ConvolutedRes_CB->fitTo(*sigToFit[c], SumW2Error(kFALSE), Range(400, 1200), RooFit::Save(kTRUE));
      else if (mass==5000)
	fitresults_CB = (RooFitResult* ) ConvolutedRes_CB->fitTo(*sigToFit[c], SumW2Error(kFALSE), Range(4000, 5500), RooFit::Save(kTRUE));
      // RooFitResult* fitresults_CB = (RooFitResult* ) ResAddPdf.fitTo(*sigToFit[c], SumW2Error(kFALSE), Range(1200, 1700), RooFit::Save(kTRUE));
      // RooFitResult* fitresults_CB = (RooFitResult* ) SigModelBW.fitTo(*sigToFit[c], SumW2Error(kFALSE), Range(1200, 1700), RooFit::Save(kTRUE));
      fitresults_CB->Print("V");
      
      RooPlot* plotOnlyResPdf = w->var("mgg")->frame(Range(1250,1700),Bins(100));
      if (mass==1500) plotOnlyResPdf = w->var("mgg")->frame(Range(1250,1700),Bins(100));
      else if (mass==750) plotOnlyResPdf = w->var("mgg")->frame(Range(400,1200),Bins(100));
      else if (mass==5000) plotOnlyResPdf = w->var("mgg")->frame(Range(4000,5500),Bins(100));
      sigToFit[c]->plotOn(plotOnlyResPdf, LineColor(kRed), LineStyle(kDashed));
      double max = plotOnlyResPdf->GetMaximum();

      RooPlot* plotPhotonsMassAll = w->var("mgg")->frame(Range(1250,1700),Bins(100));
      if (mass==1500) plotPhotonsMassAll = w->var("mgg")->frame(Range(1250,1700),Bins(100));
      else if (mass==750) plotPhotonsMassAll = w->var("mgg")->frame(Range(400,1200),Bins(100));
      else if (mass==5000) plotPhotonsMassAll = w->var("mgg")->frame(Range(4000,5500),Bins(100));
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

      TLegend *legmc = new TLegend(0.55, 0.6, 0.87, 0.88, "brNDC");
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
      // label_cms->Draw("same");
      // label_sqrt->Draw("same");
      
      int massI(mass);
      c1->SetLogy();
      plotPhotonsMassAll->GetXaxis()->SetTitle("m_{#gamma #gamma}[GeV]");
     
      c1->SetLogy(0);
      c1->SaveAs("plots/signalCBCconvBW"+TString::Format(("_M%d_cat%d_free.png"),massI, c));
      
      c1->SetLogy();
      plotPhotonsMassAll->GetYaxis()->SetRangeUser(0.01,max*10. );
      plotPhotonsMassAll->GetXaxis()->SetRangeUser(210, 290);
      c1->SaveAs("plots/signalCBCconvBW"+TString::Format(("_M%d_cat%d_LOG_free.png"),massI,c));
    }  

    /*

    // plot signal model at different widths
    bool plotW = true;
    if(plotW && c==0){
      RooRealVar var_01("var_w01", "var_w01", 0.1);
      var_01.setConstant();	
      RooFormulaVar sigmaBW_01("w01", "w01","@0", var_01);     
      RooBreitWigner SiBW_01("sigBW_01","sigBW_01" , *mgg, meanBW, sigmaBW_01);
      RooFFTConvPdf  ConvolutedRes_01("conv01", "conv01", *mgg,SiBW_01, ResAddPdf);
      
      RooRealVar var_3("var_w3", "var_w3",3);
      var_3.setConstant();	
      RooFormulaVar sigmaBW_3("w3", "w3","@0",  var_3);     
      RooBreitWigner SiBW_3("sigBW_3","sigBW_3" , *mgg, meanBW, sigmaBW_3);
      RooFFTConvPdf  ConvolutedRes_3("conv3", "conv3", *mgg,SiBW_3, ResAddPdf);
      
      RooRealVar var_6("var_w6", "var_w6", 6);
      var_6.setConstant();	
      RooFormulaVar sigmaBW_6("w6", "w6","@0", var_6);     
      RooBreitWigner SiBW_6("sigBW_6","sigBW_6" , *mgg, meanBW, sigmaBW_6);
      RooFFTConvPdf  ConvolutedRes_6("conv6", "conv6", *mgg,SiBW_6, ResAddPdf);
      
      RooRealVar var_10("var_w10", "var_w10", 10);
      var_10.setConstant();	
      RooFormulaVar sigmaBW_10("w10", "w10","@0", var_10);     
      RooBreitWigner SiBW_10("sigBW_10","sigBW_10" , *mgg, meanBW, sigmaBW_10);
      RooFFTConvPdf  ConvolutedRes_10("conv10", "conv10", *mgg,SiBW_10, ResAddPdf);
      
      RooRealVar var_15("var_w15", "var_w15",15);
      var_15.setConstant();	
      RooFormulaVar sigmaBW_15("w15", "w15","@0",  var_15);     
      RooBreitWigner SiBW_15("sigBW_15","sigBW_15" , *mgg, meanBW, sigmaBW_15);
      RooFFTConvPdf  ConvolutedRes_15("conv15", "conv15", *mgg,SiBW_15, ResAddPdf);
      
      RooPlot* plotWidths = w->var("mgg")->frame(Range(minMassFit-20,maxMassFit+20),Bins(100));
      ConvolutedRes_15.plotOn( plotWidths, LineColor(kAzure+3));
      ConvolutedRes_10.plotOn( plotWidths, LineColor(kAzure+2));
      ConvolutedRes_6.plotOn( plotWidths, LineColor(kAzure+1));
      ConvolutedRes_3.plotOn( plotWidths, LineColor(kViolet+1));
      ConvolutedRes_01.plotOn( plotWidths, LineColor(kViolet-9));
      plotWidths->Draw();
      
      label_cms->Draw("same");
      label_sqrt->Draw("same");
      
      TLegend* leg = new TLegend(0.598851,0.6044755,0.84253,0.928252,"", "brNDC");
      
      leg->SetBorderSize(0.);
      leg->SetFillColor(kWhite);
      leg->SetTextFont(42);
      plotWidths->GetYaxis()->SetRangeUser(0.001, 1.);
      plotWidths->GetXaxis()->SetTitle("m_{#gamma #gamma}[GeV]");
      plotWidths->GetYaxis()->SetTitle(" ");
      leg->AddEntry(plotWidths->getObject(0), "Width = 15 GeV", "L");
      leg->AddEntry(plotWidths->getObject(1), "Width = 10 GeV", "L");
      leg->AddEntry(plotWidths->getObject(2),"Width = 6 GeV", "L");
      leg->AddEntry(plotWidths->getObject(3),"Width = 3 GeV", "L");
      leg->AddEntry(plotWidths->getObject(4), "Width = 0.1 GeV", "L");
      leg->Draw("same");
      
      c1->SaveAs("plots/signalModels_differentWidths.png");
      c1->SaveAs("plots/signalModels_differentWidths.pdf");
    }
    */

    // IMPORTANT: fix all pdf parameters to constant
    w->defineSet(TString::Format("ConvolutedPdfParam_cat%d",c),RooArgSet( *w->var(TString::Format("MassCB_sig_sigma_cat%d",c)), 
									  *w->var(TString::Format("MassCB_sig_alphaCBpos_cat%d",c)),
									  *w->var(TString::Format("MassCB_sig_alphaCBneg_cat%d",c)),
									  *w->var(TString::Format("MassCB_sig_nCBpos_cat%d",c)),
									  *w->var(TString::Format("MassCB_sig_nCBneg_cat%d",c)),	   
									  *w->var(TString::Format("MassCB_sig_frac_cat%d",c)),  
									  *w->var(TString::Format("MassCB_sig_mean_cat%d",c)),
									  *w->var(TString::Format("sigmaBW_var_cat%d",c))));
    
    SetConstantParams(w->set(TString::Format("ConvolutedPdfParam_cat%d",c)));
    
    w->Print("V");
  }
}

// Write signal pdfs and datasets into the workspace 
void MakeSigWS(RooWorkspace* w, const char* fileBaseName, Float_t width, std::string model){
  
  TString wsDir = "workspaces/"+filePOSTfix;
  Int_t ncat = NCAT;

  RooWorkspace *wAll = new RooWorkspace("w_all","w_all");  

  //********************************//
  // Retrieve P.D.F.s
   //w->Print("V");
  for (int c=0; c<ncat; ++c) {
    //  std::cout<<"flag"<<std::endl;
    wAll->import(*w->pdf("mggSig"+TString::Format("_cat%d",c)));//*w->pdf("mggSigCBCExt"+TString::Format("_cat%d",c))
    //  std::cout<<"flag"<<std::endl;    
    //wAll->import(*w->pdf("mggSig"+TString::Format("_Inter_cat%d",c)));//*w->pdf("mggSigCBCExt"+TString::Format("_cat%d",c))
    //std::cout<<"flag"<<std::endl;
    
    wAll->import(*w->data(TString::Format("SigWeight_cat%d",c)));
    //	std::cout<<"flag"<<std::endl;
    wAll->import(*w->function("mggSig"+TString::Format("_cat%d_norm",c)));
                                                 
  }
  std::cout << "done with importing signal pdfs" << std::endl;
  wAll->import(*w->var("massReduced"));
  wAll->import(*w->var("mggGen"));
  // (2) Systematics on energy scale and resolution // chiara: per ora tutte le sistematiche non hanno senso
  // wAll->factory("CMS_hgg_sig_m0_absShift[1,1.0,1.0]");
  // wAll->factory("CMS_hgg_sig_m0_absShift_cat0[1,1.0,1.0]");
  // wAll->factory("CMS_hgg_sig_m0_absShift_cat1[1,1.0,1.0]");
  // wAll->factory("prod::CMS_hgg_sig_m0_cat0(massggnewvtx_sig_m0_cat0, CMS_hgg_sig_m0_absShift)");
  // wAll->factory("prod::CMS_hgg_sig_m0_cat1(massggnewvtx_sig_m0_cat1, CMS_hgg_sig_m0_absShift)");

  // (3) Systematics on resolution: create new sigmas
  // wAll->factory("CMS_hgg_sig_sigmaScale[1,1.0,1.0]");
  // wAll->factory("CMS_hgg_sig_sigmaScale_cat0[1,1.0,1.0]");
  // wAll->factory("CMS_hgg_sig_sigmaScale_cat1[1,1.0,1.0]");
  // wAll->factory("prod::CMS_hgg_sig_sigma_cat0(massggnewvtx_sig_sigma0_cat0, CMS_hgg_sig_sigmaScale)");
  // wAll->factory("prod::CMS_hgg_sig_gsigma_cat0(massggnewvtx_sig_sigma1_cat0, CMS_hgg_sig_sigmaScale)");
  // wAll->factory("prod::CMS_hgg_sig_sigma_cat1(massggnewvtx_sig_sigma0_cat1, CMS_hgg_sig_sigmaScale)");
  // wAll->factory("prod::CMS_hgg_sig_gsigma_cat1(massggnewvtx_sig_sigma1_cat1, CMS_hgg_sig_sigmaScale)")

  TString filename(wsDir+TString(fileBaseName)+TString::Format(("_m%.2f_w%.2f.inputsig_"+model+".root").c_str(),w->var("MH")->getVal(),width));
  wAll->writeToFile(filename);
  cout << "Write signal workspace in: " << filename << " file" << endl;
  
  return;
}

/*
// Write background pdfs and datasets into the workspace 
void MakeBkgWS(RooWorkspace* w, const char* fileBaseName, double mass) {

  TString wsDir = "workspaces/"+filePOSTfix;
  Int_t ncat = NCAT;  

  RooWorkspace *wAll = new RooWorkspace("w_all","w_all");
  
  //********************************
  // Retrieve the datasets and PDFs
  RooDataSet* data[NCAT];
 
  for (int c=0; c<ncat; ++c) {
  
    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    ((RooRealVar*) data[c]->get()->find("mgg"))->setBins(320) ;
 
    RooDataHist* dataBinned = data[c]->binnedClone();
 
    wAll->import(*w->pdf(TString::Format("mggBkg_cat%d",c)));
 
    wAll->import(*dataBinned, Rename(TString::Format("data_obs_cat%d",c)));
 
    wAll->import(*w->data(TString::Format("Data_cat%d",c)), Rename(TString::Format("data_unbinned_obs_cat%d",c)));
 
  }
  std::cout << "done with importing background pdfs" << std::endl;
  

  TString filename;
  filename = (wsDir+TString(fileBaseName)+TString::Format("_m%.2f.root",w->var("MH")->getVal()));


  wAll->writeToFile(filename);
  cout << "Write background workspace in: " << filename << " file" << endl;

  std::cout << std::endl; 
  std::cout << "observation:" << std::endl;
  for (int c=0; c<ncat; ++c) {
    std::cout << "cat " << c << ", " << wAll->data(TString::Format("data_obs_cat%d",c))->sumEntries() << endl;
    wAll->data(TString::Format("data_obs_cat%d",c))->Print();
  }
  std::cout << std::endl;
  
  return;
}
*/

void runfits(const Float_t mass=1500, string coupling="001", Bool_t dobands = false, Float_t width=0.1) {

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
 
  // RooFitResult* fitresults;
  
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
  cout << "Now prepare signal model fit - resolution function" << endl;  
  // SigModelResponseCBCBFit(w, mass, coupling);     
  // SigModelResponseDoubleCBFit(w, mass, coupling); 
  // SigModelResponseReducedCBCBFit(w, mass, coupling);       
  SigModelResponseReducedDoubleCBFit(w, mass, coupling); 

  cout << endl;
  cout << "Now try test with BW only on gen level mgg" << endl;
  // SigModelBWFit(w, mass, coupling);     
  
  cout << endl;
  cout << "Now prepare signal model fit - resolution function x BW" << endl;  
  // SigModelFitConvBW(w, mass, width);

  return;
}
