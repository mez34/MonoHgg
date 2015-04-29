using namespace RooFit;
using namespace RooStats ;

// ============================================
// to be modified:
static const Int_t NCAT = 4; 
Int_t MINmass= 500;
Int_t MAXmass= 6000;
std::string filePOSTfix="";
double signalScaler=1.00;
Float_t Lum = 19500.0;    
// ============================================

// functions
void AddSigData(RooWorkspace*, Float_t);

void AddBkgData(RooWorkspace*, Float_t);

void SigModelResponseCBCBFit(RooWorkspace*);

void SigModelFitConvBW(RooWorkspace*, Float_t, Double_t);

RooFitResult*  BkgModelFitDiJetFunc(RooWorkspace*, Bool_t, Float_t, bool);

RooFitResult*  BkgModelFitExpolFunc(RooWorkspace*, Bool_t, Float_t, bool);

RooHistFunc* getRooHistFunc(int cat, RooRealVar* var);

void SetConstantParams(const RooArgSet* params);

// CMS stuffs, cosmetics
TPaveText* get_labelCMS( int legendQuadrant, std::string year, bool sim);
TPaveText* get_labelSqrt( int legendQuadrant );

// Definition of the variables in the input ntuple
RooArgSet* defineVariables() {

  RooRealVar* mgg        = new RooRealVar("mgg",        "M(gg)",       MINmass, MAXmass, "GeV");
  RooRealVar* mggTrue    = new RooRealVar("mggTrue",    "M(gg) true",  MINmass, MAXmass, "GeV");
  RooRealVar* eventClass = new RooRealVar("eventClass", "eventClass",    -10,      10,   "");
  RooRealVar* weight     = new RooRealVar("weight",     "weightings",      0,     1000,  "");
  RooRealVar* nvtx       = new RooRealVar("nvtx",       "nvtx",            0,      50,   "");

  RooArgSet* ntplVars = new RooArgSet(*mgg, *mggTrue, *eventClass, *weight, *nvtx);
  
  return ntplVars;
}

void runfits(const Float_t mass=1500, Bool_t dobands = false, Float_t width=0.1, std::string model) {

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
 
  RooFitResult* fitresults;
  
  // import luminosity in the ws
  RooRealVar lumi("lumi","lumi",Lum);
  w->import(lumi); 
  
  // compute roohistfunc with min and max for the fit range 
  RooRealVar* var = new RooRealVar("var","var", 500, 6000 );   // maximum range

  // define the range for the fit at each mass
  RooHistFunc* rooFitMin = getRooHistFuncFitMIN(0, var);   
  RooHistFunc* rooFitMax = getRooHistFuncFitMAX(0, var);   
  TF1* fFitMin = getFuncFitMIN(0, var);
  TF1* fFitMax = getFuncFitMAX(0, var);
  makeRooHistPlot(rooFitMin,rooFitMax,fFitMin,fFitMax,var); 
  
  
  cout << endl;
  cout << endl;
  cout << "--------------------------------------------------------" << endl;
  cout << "mass range for mass " << mass << " under study: " << fFitMin->Eval(mass) << " - " << fFitMax->Eval(mass) << std::endl;
  cout << "--------------------------------------------------------" << endl;
  cout << endl;
  cout << endl;  

  var->setVal(mass);
  
  // and now set the range for the fit at this mass
  double newmin = fFitMin->Eval(mass);
  double newmax = fFitMax->Eval(mass);
  MINmass=newmin;
  MAXmass=newmax;
  Double_t MMIN = MINmass; 
  Double_t MMAX = MAXmass; 
  w->var("mgg")->setMin(MMIN);
  w->var("mgg")->setMax(MMAX);
  w->var("mggTrue")->setMin(MMIN);
  w->var("mggTrue")->setMax(MMAX);

  w->Print("v");
  
  
  cout << endl; 
  cout << "Now add signal data" << endl;
  AddSigData(w, mass);   

  cout << endl; 
  cout << "Now prepare signal model fit" << endl;  
  SigModelResponseCBCBFit(w, mass);     

  // chiara: arrivata qui e sistemati i metodi usati sopra



  /*
  cout << endl; 
  cout << "Now add background data" << endl;  
  w->var("mgg")->setMin(MINmass);         
  w->var("mgg")->setMax(MAXmass);
  AddBkgData(w, mass);     

  // fit sig
  //// SigModelFitConvBW(w, mass, width, model);      
  // SigModelFitConvRelBW(w, mass, width, model);     
  
  cout << endl; 
  // cout << "Now prepare background model fit" << endl;
  bool blind = false;                   // chiara: vedi se questi si possono mettere sopra comuni a tutto
  bool dobandsHere= false;              // chiara: vedi se questi si possono mettere sopra comuni a tutto  
  // for(int c=0; c<NCAT;c++){
  // fitresults = BkgModelFitExpPARFunc(w, dobandsHere, mass,c, blind); //MAKE FIG 3  <--
  // }

  // chiara: questi erano commentati da Livia
  // Make statistical treatment. Setup the limit on X production
  cout << endl; 
  cout << "Now make signal workspace" << endl;
  // MakeSigWS(w, fileBaseName+"_8TeV", width, model);     //<---
  cout << endl; cout << "Now make background workspace" << endl;
  // MakeBkgWS(w, fileBkgName, mass); //<--
  
  cout << endl; cout << "Now prepare datacards" << endl;
  // for (int c=0; c<NCAT; c++) MakeDataCard_1Channel(w, fileBaseName, fileBkgName,width,c, model); //<--

  */

  return;
}

void AddSigData(RooWorkspace* w, Float_t mass) {
  
  Int_t ncat = NCAT;
  TString inDir = "/afs/cern.ch/work/c/crovelli/Flashgg/EXO_7_4_0_pre9/src/Diphotons/analysis/python/";
  
  // Variables
  RooArgSet* ntplVars = defineVariables();
  // ntplVars->add(*w->var("mgg"));       // chiara: questo serve?        
  // (*w->var("mgg")).setRange(130,230);  // chiara: questo serve?
  
  // Files
  int iMass = abs(mass);   
  TChain* sigTree = new TChain();
  sigTree->Add(inDir+TString(Form("diPhotons_RS_%d.root/DiPhotonTree", iMass)));
  sigTree->SetTitle("sigTree");
  sigTree->SetName("sigTree");

  // common preselection cut
  TString mainCut = "mgg>=500 && mgg<=6000";   
  
  // Create signal dataset 
  RooDataSet sigWeighted("sigWeighted","dataset",sigTree,*ntplVars,mainCut,"weight");
  cout << endl;
  cout << "sigWeighted" << endl;
  sigWeighted.Print("v");
  cout << "---- nX:  " << sigWeighted.sumEntries() << endl; 
  
  // apply a common preselection cut; split in categories
  cout << endl;
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
}

// Data dataset
void AddBkgData(RooWorkspace* w, Float_t mass) {

  // initializations
  Int_t ncat = NCAT;

  Float_t minMassFit, maxMassFit;
  minMassFit = MINmass;
  maxMassFit = MAXmass;

  // retrieve the data tree
  TFile dataFile("/afs/cern.ch/work/c/crovelli/Flashgg/EXO_7_4_0_pre9/src/Diphotons/analysis/python/diPhotons_GGj.root");  // chiara
  TTree* dataTree = (TTree*) dataFile.Get("DiPhotonTree");

  // Variables
  RooArgSet* ntplVars = defineVariables();
  ntplVars->add(*w->var("mgg"));
  (*w->var("mgg")).setRange(MINmass, MAXmass);
  
  // common preselection cut
  TString mainCut = TString::Format(" mgg>=(%.1f) && mgg<=(%.1f)", minMassFit, maxMassFit);   
  
  // Create dataset
  RooDataSet Data("Data","dataset",dataTree,*ntplVars,mainCut,"weight");
  cout << endl;
  cout << "Data, everything: " << endl;
  Data.Print("v");
  cout << "---- nX:  " << Data.sumEntries() << endl;
  cout << endl;

  // split into NCAT categories; apply common preselection cut
  RooDataSet* dataToFit[NCAT];  
  for (int c=0; c<ncat; ++c) {

    if (c==0) dataToFit[c] = (RooDataSet*) Data.reduce(*w->var("mgg"),mainCut+TString::Format("&& eventClass==0"));
    if (c==1) dataToFit[c] = (RooDataSet*) Data.reduce(*w->var("mgg"),mainCut+TString::Format("&& eventClass==1"));
    if (c==2) dataToFit[c] = (RooDataSet*) Data.reduce(*w->var("mgg"),mainCut+TString::Format("&& eventClass==2"));
    if (c==3) dataToFit[c] = (RooDataSet*) Data.reduce(*w->var("mgg"),mainCut+TString::Format("&& eventClass==3"));

    cout << endl; cout << "for category = " << c << endl;
    dataToFit[c]->Print("v");
    cout << "---- nX:  " << dataToFit[c]->sumEntries() << endl;

     w->import(*dataToFit[c],Rename(TString::Format("Data_cat%d",c)));
  }

  // Create full data set without categorization
  cout << "data, no split" << endl;
  RooDataSet* data = (RooDataSet*) Data.reduce(*w->var("mgg"),mainCut);
  w->import(*data, Rename("Data"));
  cout << endl;
 
  data->Print("v");
  cout << "---- nX:  " << data->sumEntries() << endl; 
}


// ok
// Signal model: sum of two CB with same mass and mean
void SigModelResponseCBCBFit(RooWorkspace* w, Float_t mass) {

  int iMass = abs(mass);   

  // Files
  TString inDir = "/afs/cern.ch/work/c/crovelli/Flashgg/EXO_7_4_0_pre9/src/Diphotons/analysis/python/";
  TChain* sigTree = new TChain();
  sigTree->Add(inDir+TString(Form("diPhotons_RS_%d.root/DiPhotonTree", iMass)));
  sigTree->SetTitle("sigTree");
  sigTree->SetName("sigTree");

  // Variables
  RooArgSet* ntplVars = defineVariables();
  
  // common preselection cut
  TString mainCut1 = TString::Format("mgg>=500 && mgg<=6000");   
  RooDataSet sigWeighted("sigWeighted","dataset",sigTree,*ntplVars,mainCut1,"weight");
  RooRealVar* mgg = w->var("mgg");     

  // reduced mass
  RooFormulaVar *massReduced_formula = new RooFormulaVar("massReduced_formula","","@0/@1 -1",RooArgList(*w->var("mgg"),*w->var("mggTrue")));
  RooRealVar* massReduced = (RooRealVar*) sigWeighted.addColumn(*massReduced_formula);
  massReduced->SetName("massReduced");
  massReduced->SetTitle("massReduced");
  w->import(*massReduced);  
  massReduced->setRange(-0.5, 0.5);
  
  // common preselection cut on the reduced mass
  TString mainCut = TString::Format("massReduced>-0.5 && massReduced <0.5"); 
  
  // fit functions
  RooDataSet* signal[NCAT];
  RooCBShape* ResponseCBpos[NCAT];
  RooCBShape* ResponseCBneg[NCAT];
  RooGaussian* ResponseGauss[NCAT];
  RooAddPdf* ResponseAddGauss[NCAT];
  RooAddPdf* ResponseAdd[NCAT];
  
  TCanvas* c1 = new TCanvas("c1", "c1", 1);
  c1->cd();
  TPaveText* label_cms  = get_labelCMS(0, "2015", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
  
  // chiara: sistematiche, commento per ora
  // Double_t scaleSyst;
  // Double_t smearSyst;
  
  for(int c = 0; c<NCAT; c++){
    
    TLatex *lat  = new TLatex(0.6,0.9,TString::Format("Cat: %d", c));  
    lat->SetTextSize(0.038);
    lat->SetTextAlign(11);
    lat->SetTextFont(42); 
    lat->SetNDC();

    // chiara: sistematiche, commento per ora    
    // RooRealVar *mShift = new RooRealVar(TString::Format("mShift_cat%d", c), TString::Format("mShift_cat%d",c), -10., 10.);
    // mShift->setVal(0.);
    ////  mShift->setConstant();
    // w->import(*mShift);
    
    // RooRealVar *mSmear = new RooRealVar(TString::Format("mSmear_cat%d", c), TString::Format("mSmear_cat%d",c), -10., 10.);
    // mSmear->setVal(0.);
    ////  mSmear->setConstant();
    // w->import(*mSmear);
    
    w->Print();
    
    // splitting in categories
    if (c==0) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& eventClass==0"));
    if (c==1) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& eventClass==1"));
    if (c==2) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& eventClass==2"));
    if (c==3) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& eventClass==3"));
    w->import(*signal[c],Rename(TString::Format("SigWeightReduced_cat%d",c))); 
    
    // chiara: sistematiche, commento per ora  
    // if (c==0 || c==1) scaleSyst = 0.005;
    // if (c==2 || c==3) scaleSyst = 0.007;
    // if (c==0) smearSyst = 0.005;
    // if (c==1) smearSyst = 0.0058;
    // if (c==2 || c==3) smearSyst = 0.01;

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
    RooFitResult* fitresults = (RooFitResult* ) ResponseAdd[c]->fitTo(*signal[c], SumW2Error(kTRUE), Range(-1, 1), RooFit::Save(kTRUE));
    std::cout<<TString::Format("******************************** Signal Fit results CB+CB  mass %f cat %d***********************************", mass, c)<<std::endl;
    fitresults->Print("V");
   
    
    // Plot
    RooPlot* plotG = massReduced->frame(Range(-0.12, 0.12),Title("Mass Reduced"),Bins(60));
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

    TLegend* legmc = new TLegend(0.6, 0.62, 0.91, 0.95, "", "brNDC");
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
    
    c1->SetLogy();
    c1->SaveAs(TString::Format("plots/responseFitCBCB_cat%d_LOG.png",c));
    c1->SaveAs(TString::Format("plots/responseFitCBCB_cat%d_LOG.pdf",c));
	       
    w->defineSet(TString::Format("ResponseAddPdfParam_cat%d",c),RooArgSet(*w->var(TString::Format("ReducedMass_sig_sigma_cat%d",c)), 
									  *w->var(TString::Format("ReducedMass_sig_alphaCBpos_cat%d",c)),
									  *w->var(TString::Format("ReducedMass_sig_alphaCBneg_cat%d",c)),
									  *w->var(TString::Format("ReducedMass_sig_nCBpos_cat%d",c)),
									  *w->var(TString::Format("ReducedMass_sig_nCBneg_cat%d",c)),	   
									  *w->var(TString::Format("ReducedMass_sig_frac_cat%d",c)),  
									  *w->var(TString::Format("ReducedMass_sig_mean_cat%d",c))));
    SetConstantParams(w->set(TString::Format("ResponseAddPdfParam_cat%d",c)));
  }
}


// Fit signal with model with CB convoluted with BW
// chiara: qui ci sono diverse width, da capire se andra' tenuto o no
void SigModelFitConvBW(RooWorkspace* w, Float_t mass, Double_t width, std::string model) {

  Int_t ncat = NCAT;
  RooDataSet* sigToFit[NCAT];
  
  Float_t MASS(mass);  
  Float_t minMassFit(mass*0.8);
  Float_t maxMassFit(mass*1.2); 

  RooRealVar* mgg = w->var("mgg"); 

  // chiara: sistematiche, per il momento commento
  // Double_t scaleSyst;
  // Double_t smearSyst;

  std::cout<<"----------------------------------------------------------------------------------------"<<std::endl;

  //  mgg->setRange("sigrange",minMassFit-20,maxMassFit+20); 

  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
  TFile* f = new TFile("sigShapeCorrections.root", "READ");

  // Fit to Signal 
  for (int c=0; c<ncat; ++c) {
    cout << "---------- Category = " << c << endl;

    // chiara: sicuri che va lasciato commentato?
    // sigToFit[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c));
    // w->import(*sigToFit[c]);

    // chiara: sistematiche, per il momento commento    
    // introduce systs 
    // if (c==0 || c==1) scaleSyst = 0.005;
    // if (c==2 || c==3) scaleSyst = 0.007;
    // if (c==0)         smearSyst = 0.005;
    // if (c==1)         smearSyst = 0.0058;
    // if (c==2 || c==3) smearSyst = 0.01;

    // get sigma from TF1:   
    TF1* fcn = (TF1*)f->Get(TString::Format("f%d",c));
    Float_t massF = (Float_t) (*w->var("MH")).getVal();
    Float_t sigmaCorr = fcn->Eval(massF);
    if (massF==150) sigmaCorr=1;
    RooRealVar rooSigmaCorr (TString::Format("rooSigmaCorr_cat%d",c), TString::Format("rooSigmaCorr_cat%d",c), sigmaCorr, "");
    rooSigmaCorr.setConstant();
    w->import(rooSigmaCorr);
    ( *w->var(TString::Format("mShift_cat%d",c))).setConstant();
    ( *w->var(TString::Format("mSmear_cat%d",c))).setConstant();
    

    // cb
    RooFormulaVar CBpos_mean_draw(TString::Format("CBpos_mean_draw_cat%d",c),"","@0+@1",RooArgList(*w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)),*w->var("MH")));
    RooFormulaVar CBpos_mean(TString::Format("CBpos_mean_cat%d",c),"",TString::Format("@0+%f*@1", scaleSyst),RooArgList(*w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)),*w->var(TString::Format("mShift_cat%d",c))));
    RooFormulaVar CBpos_sigma(TString::Format("CBpos_sigma_cat%d",c),"",TString::Format("(sqrt(@0*@0*@3*@3+%f*%f*@2)*@1)",smearSyst,smearSyst),RooArgList(*w->var(TString::Format("ReducedMass_sig_sigma_cat%d",c)),*w->var("MH"),*w->var(TString::Format("mSmear_cat%d",c)),*w->var(TString::Format("rooSigmaCorr_cat%d",c)) ) );
    // std::cout <<"-------------------> SIGMA: "<<CBpos_sigma->getVal()<<"    MASS: "<<(*w->var("MH")).getVal()<<std::endl;
    
    RooFormulaVar CBpos_alphaCB(TString::Format("CBpos_alphaCB_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_alphaCBpos_cat%d",c)) );
    RooFormulaVar CBneg_alphaCB(TString::Format("CBneg_alphaCB_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_alphaCBneg_cat%d",c)) );
    RooFormulaVar CBpos_n(TString::Format("CBpos_n_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_Npos_cat%d",c)) );
    RooFormulaVar CBneg_n(TString::Format("CBneg_n_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_Nneg_cat%d",c)) );
    RooFormulaVar CBpos_frac(TString::Format("CBpos_frac_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_frac_cat%d",c)) );
    
    RooCBShape ResCBpos_draw(TString::Format("ResCBpos_draw_cat%d",c),TString::Format("ResCBpos_draw_cat%d",c) , *mgg, CBpos_mean_draw, CBpos_sigma, CBpos_alphaCB, CBpos_n) ;
    RooCBShape ResCBneg_draw(TString::Format("ResCBneg_draw_cat%d",c),TString::Format("ResCBneg_draw_cat%d",c) , *mgg, CBpos_mean_draw, CBpos_sigma, CBneg_alphaCB, CBneg_n) ;
    RooCBShape ResCBpos(TString::Format("ResCBpos_cat%d",c),TString::Format("ResCBpos_cat%d",c) , *mgg, CBpos_mean, CBpos_sigma,CBpos_alphaCB, CBpos_n) ;
    RooCBShape ResCBneg(TString::Format("ResCBneg_cat%d",c),TString::Format("ResCBneg_cat%d",c) , *mgg, CBpos_mean, CBpos_sigma,CBneg_alphaCB, CBneg_n) ;
    mgg->setBins(40000, "cache");  

    RooAddPdf ResAddPdf(TString::Format("ResAddPdf_cat%d",c),TString::Format("ResAddPdf_cat%d",c) , RooArgList(ResCBneg, ResCBpos), CBpos_frac);
    RooAddPdf ResAddPdf_draw(TString::Format("ResAddPdf_draw_cat%d",c),TString::Format("ResAddPdf_draw_cat%d",c) , RooArgList(ResCBneg_draw, ResCBpos_draw), CBpos_frac);


    // BW
    RooFormulaVar meanBW(TString::Format("meanBW_cat%d",c),"","@0",*w->var("MH"));  
    RooRealVar sigmaBW_var(TString::Format("sigmaBW_var_cat%d",c), TString::Format("sigmaBW_var_cat%d",c), width);
    // std::cout<<" width:--------> "<<width<<std::endl;
    sigmaBW_var.setConstant();
    w->import(sigmaBW_var);
    
    RooFormulaVar* sigmaBW;
    if(width<1)sigmaBW = new RooFormulaVar(TString::Format("sigmaBW_cat%d",c),"","@0",*w->var(TString::Format("sigmaBW_var_cat%d",c))); 
    else if(width==2)  sigmaBW = new RooFormulaVar(TString::Format("sigmaBW_cat%d",c),"","@0*0.02",*w->var("MH"));   
    else if(width==5)  sigmaBW = new RooFormulaVar(TString::Format("sigmaBW_cat%d",c),"","@0*0.05",*w->var("MH"));   
    else if(width==7)  sigmaBW = new RooFormulaVar(TString::Format("sigmaBW_cat%d",c),"","@0*0.07",*w->var("MH"));   
    else if(width==10) sigmaBW = new RooFormulaVar(TString::Format("sigmaBW_cat%d",c),"","@0*0.10",*w->var("MH"));   
    RooBreitWigner SigModelBW(TString::Format("BW_cat%d",c),TString::Format("BW_cat%d",c), *mgg, meanBW, *sigmaBW);

    /*
    // correction to BW
    TFile* file = new TFile("BW_corrections.root","READ");
    
    TGraph2D* g0 = (TGraph2D*)file->Get(TString::Format("p0_cat%d",c));   
    double var1 = g0->Interpolate(mass, width);  
    RooRealVar v1(TString::Format("v1_cat%d",c),TString::Format("v1_cat%d",c), var1);
    v1.setConstant();
    std::cout<<"++++++++++++++++++++++++++ "<<var1<<std::endl; 
    w->import(v1);
 
    RooFormulaVar f1(TString::Format("f1_cat%d",c), TString::Format("f1_cat%d",c), "@0",*w->var(TString::Format("v1_cat%d",c)));
    
    RooGenericPdf pf(TString::Format("pf_cat%d",c),  "exp(@1*(@0-@2))", RooArgList(*mgg, f1,*w->var("MH")));
    RooProdPdf* prod;
    */

    // Convolution
    RooFFTConvPdf*  ConvolutedRes_CB;
    ConvolutedRes_CB = new RooFFTConvPdf(TString::Format("mggSig_cat%d",c),TString::Format("mggSig_cat%d",c), *mgg,SigModelBW, ResAddPdf);
    w->import(*ConvolutedRes_CB);
    
    RooHistFunc* rooFunc_norm = getRooHistFunc(c,w->var("MH"), model );
    w->import(*rooFunc_norm);

    std::cout<<"SIG NORM ----->"<<rooFunc_norm->getVal(*w->var("MH"));
    // w->Print("V");

    if(width <2. && mass < 150){ //if i want to plot the fit
      sigToFit[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c));
      
      RooFitResult* fitresults_CB = (RooFitResult* ) ConvolutedRes_CB.fitTo(*sigToFit[c], RooFit::Save(kTRUE));
      fitresults_CB->Print("V");
      
      RooPlot* plotOnlyResPdf = w->var("mgg")->frame(Range(minMassFit-20,maxMassFit+20),Bins(100));
      sigToFit[c]->plotOn(plotOnlyResPdf, LineColor(kRed), LineStyle(kDashed));
      double max = plotOnlyResPdf->GetMaximum();
      
      // Plot to verify everything is ok
      RooPlot* plotPhotonsMassAll = w->var("mgg")->frame(Range(minMassFit-20,maxMassFit+20),Bins(100));
      sigToFit[c]->plotOn(plotPhotonsMassAll);
      SigModelBW.plotOn(plotPhotonsMassAll, LineColor(kGreen), LineStyle(kDashed));
      //  ResAddPdf_draw.plotOn(plotPhotonsMassAll, LineColor(kRed), LineStyle(kDashed));
      ConvolutedRes_CB.plotOn(plotPhotonsMassAll, LineColor(kBlue));

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

      TLegend *legmc = new TLegend(0.55, 0.6, 0.87, 0.88, ("Model: "+model).c_str(), "brNDC");
      legmc->AddEntry(plotPhotonsMassAll->getObject(0),"Simulation","LPE");
      legmc->AddEntry(plotPhotonsMassAll->getObject(1),"BW","L");
      // legmc->AddEntry(plotPhotonsMassAll->getObject(2)," CB + CB ","L");
      legmc->AddEntry(plotPhotonsMassAll->getObject(2),"BW #otimes Resolution","L");
      legmc->SetTextSize(0.0206044);
      legmc->SetTextFont(42);
      legmc->SetBorderSize(0);
      legmc->SetFillStyle(0);
      legmc->Draw();
      lat->Draw("same");
      label_cms->Draw("same");
      label_sqrt->Draw("same");
      
      int massI(mass);
      c1->SetLogy();
      
      plotPhotonsMassAll->GetXaxis()->SetTitle("m_{#gamma #gamma}[GeV]");
     
      c1->SetLogy(0);
      c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_"+model+".png").c_str(),massI, c));
      c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_"+model+".pdf").c_str(),massI, c));
      
      c1->SetLogy();
      plotPhotonsMassAll->GetYaxis()->SetRangeUser(0.01,max*10. );
      plotPhotonsMassAll->GetXaxis()->SetRangeUser(210, 290);
      c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_LOG_"+model+".png").c_str(),massI,c));
      c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_LOG_"+model+".pdf").c_str(),massI,c));
    }  

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
      
      c1->SaveAs("plots/SignalModels_differentWidths.png");
      c1->SaveAs("plots/SignalModels_differentWidths.pdf");
    }
    
    // IMPORTANT: fix all pdf parameters to constant
    w->defineSet(TString::Format("ConvolutedPdfParam_cat%d",c),RooArgSet( *w->var(TString::Format("ReducedMass_sig_sigma_cat%d",c)), 
									  *w->var(TString::Format("ReducedMass_sig_alphaCBpos_cat%d",c)),
									  *w->var(TString::Format("ReducedMass_sig_alphaCBneg_cat%d",c)),
									  *w->var(TString::Format("ReducedMass_sig_Npos_cat%d",c)),
									  *w->var(TString::Format("ReducedMass_sig_Nneg_cat%d",c)),	   
									  *w->var(TString::Format("ReducedMass_sig_frac_cat%d",c)),  
									  *w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)),
									  *w->var(TString::Format("sigmaBW_var_cat%d",c))));
    
    SetConstantParams(w->set(TString::Format("ConvolutedPdfParam_cat%d",c)));
    
    //w->Print("V");
  }
}

// Fit signal with model with CB convoluted with BW 
void SigModelFitConvRelBW(RooWorkspace* w, Float_t mass, Double_t width, std::string model) {

  Int_t ncat = NCAT;
  RooDataSet* sigToFit[NCAT];
  
  Float_t MASS(mass);  
  Float_t minMassFit(mass*0.8);
  Float_t maxMassFit(mass*1.2); 

  RooRealVar* mgg = w->var("mgg"); 

  // chiara: sistematiche, per il momento commento     
  // Double_t scaleSyst;
  // Double_t smearSyst;

  std::cout<<"----------------------------------------------------------------------------------------"<<std::endl;

  //  mgg->setRange("sigrange",minMassFit-20,maxMassFit+20); 

  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
  TFile* f = new TFile("sigShapeCorrections.root", "READ");

  // Fit Signal 
  for (int c=0; c<ncat; ++c) {
    cout << "---------- Category = " << c << endl;

    // chiara: sicuri che va lasciato commentato?        
    // sigToFit[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c));
    // w->import(*sigToFit[c]);

    // chiara: sistematiche, per il momento commento  
    // introduce systs 
    // if (c==0 || c==1) scaleSyst = 0.005;
    // if (c==2 || c==3) scaleSyst = 0.007;
    // if (c==0)         smearSyst = 0.005;
    // if (c==1)         smearSyst = 0.0058;
    // if (c==2 || c==3) smearSyst = 0.01;

    // get sigma from TF1:   - chiara: da capire
    TF1* fcn = (TF1*)f->Get(TString::Format("f%d",c));
    Float_t massF = (Float_t) (*w->var("MH")).getVal();
    Float_t sigmaCorr = fcn->Eval(massF);
    if(massF==150) sigmaCorr=1;
    RooRealVar rooSigmaCorr (TString::Format("rooSigmaCorr_cat%d",c), TString::Format("rooSigmaCorr_cat%d",c), sigmaCorr, "");
    rooSigmaCorr.setConstant();
    w->import(rooSigmaCorr);
    ( *w->var(TString::Format("mShift_cat%d",c))).setConstant();
    ( *w->var(TString::Format("mSmear_cat%d",c))).setConstant();

    // cb
    RooFormulaVar CBpos_mean_draw(TString::Format("CBpos_mean_draw_cat%d",c),"","@0+@1",RooArgList(*w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)),*w->var("MH")));
    RooFormulaVar CBpos_mean(TString::Format("CBpos_mean_cat%d",c),"",TString::Format("@0+%f*@1", scaleSyst),RooArgList(*w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)),*w->var(TString::Format("mShift_cat%d",c))));
    RooFormulaVar CBpos_sigma(TString::Format("CBpos_sigma_cat%d",c),"",TString::Format("(sqrt(@0*@0*@3*@3+%f*%f*@2)*@1)",smearSyst,smearSyst),RooArgList(*w->var(TString::Format("ReducedMass_sig_sigma_cat%d",c)),*w->var("MH"),*w->var(TString::Format("mSmear_cat%d",c)),*w->var(TString::Format("rooSigmaCorr_cat%d",c)) ) );
    //    std::cout<<"-------------------> SIGMA: "<<CBpos_sigma->getVal()<<"    MASS: "<<(*w->var("MH")).getVal()<<std::endl;

    RooFormulaVar CBpos_alphaCB(TString::Format("CBpos_alphaCB_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_alphaCBpos_cat%d",c)) );
    RooFormulaVar CBneg_alphaCB(TString::Format("CBneg_alphaCB_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_alphaCBneg_cat%d",c)) );
    RooFormulaVar CBpos_n(TString::Format("CBpos_n_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_Npos_cat%d",c)) );
    RooFormulaVar CBneg_n(TString::Format("CBneg_n_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_Nneg_cat%d",c)) );
    RooFormulaVar CBpos_frac(TString::Format("CBpos_frac_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_frac_cat%d",c)) );
    
    RooCBShape ResCBpos_draw(TString::Format("ResCBpos_draw_cat%d",c),TString::Format("ResCBpos_draw_cat%d",c) , *mgg, CBpos_mean_draw, CBpos_sigma, CBpos_alphaCB, CBpos_n) ;
    RooCBShape ResCBneg_draw(TString::Format("ResCBneg_draw_cat%d",c),TString::Format("ResCBneg_draw_cat%d",c) , *mgg, CBpos_mean_draw, CBpos_sigma, CBneg_alphaCB, CBneg_n) ;
    RooCBShape ResCBpos(TString::Format("ResCBpos_cat%d",c),TString::Format("ResCBpos_cat%d",c) , *mgg, CBpos_mean, CBpos_sigma,CBpos_alphaCB, CBpos_n) ;
    RooCBShape ResCBneg(TString::Format("ResCBneg_cat%d",c),TString::Format("ResCBneg_cat%d",c) , *mgg, CBpos_mean, CBpos_sigma,CBneg_alphaCB, CBneg_n) ;
    mgg->setBins(40000, "cache");  

    RooAddPdf ResAddPdf(TString::Format("ResAddPdf_cat%d",c),TString::Format("ResAddPdf_cat%d",c) , RooArgList(ResCBneg, ResCBpos), CBpos_frac);
    RooAddPdf ResAddPdf_Inter(TString::Format("ResAddPdf_Inter_cat%d",c),TString::Format("ResAddPdf_Inter_cat%d",c) , RooArgList(ResCBneg, ResCBpos), CBpos_frac);
    RooAddPdf ResAddPdf_draw(TString::Format("ResAddPdf_draw_cat%d",c),TString::Format("ResAddPdf_draw_cat%d",c) , RooArgList(ResCBneg_draw, ResCBpos_draw), CBpos_frac);


    // BW
    RooFormulaVar meanBW(TString::Format("meanBW_cat%d",c),"","@0",*w->var("MH"));  
    RooRealVar sigmaBW_var(TString::Format("sigmaBW_var_cat%d",c), TString::Format("sigmaBW_var_cat%d",c), width);
    std::cout<<" width:--------> "<<width<<std::endl;
    sigmaBW_var.setConstant();
    w->import(sigmaBW_var);
    
    RooFormulaVar* sigmaBW;
    if(width<1)        sigmaBW = new RooFormulaVar(TString::Format("sigmaBW_cat%d",c),"","@0",*w->var(TString::Format("sigmaBW_var_cat%d",c))); 
    else if(width==2)  sigmaBW = new RooFormulaVar(TString::Format("sigmaBW_cat%d",c),"","@0*0.02",*w->var("MH"));   
    else if(width==5)  sigmaBW = new RooFormulaVar(TString::Format("sigmaBW_cat%d",c),"","@0*0.05",*w->var("MH"));   
    else if(width==7)  sigmaBW = new RooFormulaVar(TString::Format("sigmaBW_cat%d",c),"","@0*0.07",*w->var("MH"));   
    else if(width==10) sigmaBW = new RooFormulaVar(TString::Format("sigmaBW_cat%d",c),"","@0*0.10",*w->var("MH"));   


    // Convolution
    RooGenericPdf SigModelBW(TString::Format("BW_cat%d",c),"1 / ( TMath::Power( TMath::Power(@0,2) - TMath::Power(@1,2) , 2 ) + TMath::Power(@1,2)*TMath::Power(@2,2) )", RooArgSet(*mgg, meanBW, *sigmaBW));
    RooFFTConvPdf*  ConvolutedRes_CB;    
    ConvolutedRes_CB = new RooFFTConvPdf(TString::Format("mggSig_cat%d",c),TString::Format("mggSig_cat%d",c), *mgg,SigModelBW, ResAddPdf);
    w->import(*ConvolutedRes_CB);
    
    // RooHistFunc* rooFunc_norm = getRooHistFunc(c,w->var("MH"), model );
    RooHistFunc* rooFunc_norm = getNorm2D(c,w->var("MH"), w->var("MH")->getVal(), width,model);
    w->import(*rooFunc_norm);
    std::cout<<"SIG NORM ----->"<<rooFunc_norm->getVal(*w->var("MH"))<<std::endl;

    // RooHistFunc* rooFunc_norm2D_0 = getNorm2D(c,w->var("MH"),w->var("MH")->getVal(), 0);
    // RooHistFunc* rooFunc_norm2D_2 = getNorm2D(c,w->var("MH"),w->var("MH")->getVal(),2);
    // RooHistFunc* rooFunc_norm2D_5 = getNorm2D(c,w->var("MH"),w->var("MH")->getVal(),5);
    // RooHistFunc* rooFunc_norm2D_10 = getNorm2D(c,w->var("MH"),w->var("MH")->getVal(), 10);
    // std::cout<<"SIG NORM ----->"<<rooFunc_norm2D_0->getVal(*w->var("MH"))<<"      "<<rooFunc_norm2D_2->getVal(*w->var("MH"))<<"  "<<rooFunc_norm2D_5->getVal(*w->var("MH"))<<"      "<<rooFunc_norm2D_10->getVal(*w->var("MH"))<<std::endl;;*/
    // w->Print("V");

    if(width <2. && mass <= 150){ //if i want to plot the fit
      sigToFit[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c));

      RooFitResult* fitresults_CB = (RooFitResult* ) ConvolutedRes_CB.fitTo(*sigToFit[c], RooFit::Save(kTRUE));
      fitresults_CB->Print("V");
      
      RooPlot* plotOnlyResPdf = w->var("mgg")->frame(Range(minMassFit-20,maxMassFit+20),Bins(100));
      sigToFit[c]->plotOn(plotOnlyResPdf, LineColor(kRed), LineStyle(kDashed));
      double max = plotOnlyResPdf->GetMaximum();
      
      // Plot to verify everything is ok
      RooPlot* plotmggAll = w->var("mgg")->frame(Range(minMassFit-20,maxMassFit+20),Bins(100));
      sigToFit[c]->plotOn(plotmggAll);
      SigModelBW.plotOn(plotmggAll, LineColor(kGreen), LineStyle(kDashed));
      //  ResAddPdf_draw.plotOn(plotmggAll, LineColor(kRed), LineStyle(kDashed));
      ConvolutedRes_CB.plotOn(plotmggAll, LineColor(kBlue));

      TCanvas* c1 = new TCanvas("c1","mgg",0,0,800,800);
      c1->cd(1);
      plotPhotonsMassAll->Draw();  
      plotPhotonsMassAll->GetYaxis()->SetRangeUser(0.01, max*1.2);
      plotPhotonsMassAll->GetXaxis()->SetTitleFont(42);
      plotPhotonsMassAll->GetXaxis()->SetTitleSize(0.05);
      plotPhotonsMassAll->GetYaxis()->SetTitleFont(42);
      plotPhotonsMassAll->GetYaxis()->SetTitleSize(0.05);
      TLatex *lat  = new TLatex(0.55,0.9,TString::Format("Cat: %d", c));  
      lat->SetTextSize(0.038);
      lat->SetTextAlign(11);
      lat->SetTextFont(42); 
      lat->SetNDC();

      TLatex* latex = new TLatex(0.21, 0.76, "#splitline{m_{X}=150 GeV}{#splitline{}{Class 0}}");
      latex->SetTextSize(0.038);
      latex->SetTextAlign(11);
      latex->SetTextFont(42); 
      latex->SetNDC();

      TLegend *legmc = new TLegend(0.55, 0.69, 0.88, 0.95);
      legmc->AddEntry(plotPhotonsMassAll->getObject(0),"Simulation","PL");
      legmc->AddEntry(plotPhotonsMassAll->getObject(1),"BW","L");
      legmc->AddEntry(plotPhotonsMassAll->getObject(2),"BW #otimes Resolution","L");
      legmc->SetTextSize(0.0286044);
      legmc->SetTextFont(42);
      legmc->SetBorderSize(0);
      legmc->SetFillStyle(0);
      legmc->Draw();
      latex->Draw("same");
      // label_cms->Draw("same");
      // label_sqrt->Draw("same");
      int iPos=11 ;
      // CMS_lumi( c1,true,iPos );

      int massI(mass);
      c1->SetLogy();
     
      plotPhotonsMassAll->GetXaxis()->SetTitle("m_{#gamma#gamma}[GeV]"); 
      plotPhotonsMassAll->GetYaxis()->SetTitle("Events/1 GeV");
     
      c1->SetLogy(0);
      c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_"+model+".png").c_str(),massI, c));
      c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_"+model+".pdf").c_str(),massI, c));
      if(c==0){
	c1->SaveAs("~/www/plotsPAS/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_"+model+".png").c_str(),massI, c));
	c1->SaveAs("~/www/plotsPAS/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_"+model+".pdf").c_str(),massI, c));
	c1->SaveAs("~/www/plotsPAS/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_"+model+".C").c_str(),massI, c));
	c1->SaveAs("~/www/plotsPAS/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_"+model+".root").c_str(),massI, c));
      }
      c1->SetLogy();
      plotPhotonsMassAll->GetYaxis()->SetRangeUser(0.01,max*10. );
      plotPhotonsMassAll->GetXaxis()->SetRangeUser(210, 290);
      c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_LOG_"+model+".png").c_str(),massI,c));
      c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_LOG_"+model+".pdf").c_str(),massI,c));
      if(c==0){
	c1->SaveAs("~/www/plotsPAS/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_LOG_"+model+".png").c_str(),massI,c));
	c1->SaveAs("~/www/plotsPAS/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_LOG_"+model+".pdf").c_str(),massI,c));
	c1->SaveAs("~/www/plotsPAS/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_LOG_"+model+".C").c_str(),massI,c));
	c1->SaveAs("~/www/plotsPAS/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_LOG_"+model+".root").c_str(),massI,c));
      }
    }

    //plot signal model at different widths
    bool plotW = false;
    if(plotW && c==0){
      RooRealVar var_01("var_w01", "var_w01", 0.1);
      var_01.setConstant();	
      RooFormulaVar sigmaBW_01("w01", "w01","@0", var_01); 
      //	RooBreitWigner SiBW_01("sigBW_01","sigBW_01" , *mgg, meanBW, sigmaBW_01);
      RooGenericPdf SiBW_01("sigBW_01","1 / ( TMath::Power( TMath::Power(@0,2) - TMath::Power(@1,2) , 2 ) + TMath::Power(@1,2)*TMath::Power(@2,2) )", RooArgSet(*mgg, meanBW, sigmaBW_01));    
      RooFFTConvPdf  ConvolutedRes_01("conv01", "conv01", *mgg,SiBW_01, ResAddPdf);
      
      RooRealVar var_3("var_w3", "var_w3",3);
      var_3.setConstant();	
      RooFormulaVar sigmaBW_3("w3", "w3","@0",  var_3);     
      //	RooBreitWigner SiBW_3("sigBW_3","sigBW_3" , *mgg, meanBW, sigmaBW_3);
      RooGenericPdf SiBW_3("sigBW_3","1 / ( TMath::Power( TMath::Power(@0,2) - TMath::Power(@1,2) , 2 ) + TMath::Power(@1,2)*TMath::Power(@2,2) )", RooArgSet(*mgg, meanBW, sigmaBW_3));    
      RooFFTConvPdf  ConvolutedRes_3("conv3", "conv3", *mgg,SiBW_3, ResAddPdf);
      
      RooRealVar var_6("var_w6", "var_w6", 6);
      var_6.setConstant();	
      RooFormulaVar sigmaBW_6("w6", "w6","@0", var_6);     
      //RooBreitWigner SiBW_6("sigBW_6","sigBW_6" , *mgg, meanBW, sigmaBW_6);
      RooGenericPdf SiBW_6("sigBW_6","1 / ( TMath::Power( TMath::Power(@0,2) - TMath::Power(@1,2) , 2 ) + TMath::Power(@1,2)*TMath::Power(@2,2) )", RooArgSet(*mgg, meanBW, sigmaBW_3));    
      RooFFTConvPdf  ConvolutedRes_6("conv6", "conv6", *mgg,SiBW_6, ResAddPdf);
      
      RooRealVar var_10("var_w10", "var_w10", 10);
      var_10.setConstant();	
      RooFormulaVar sigmaBW_10("w10", "w10","@0", var_10);     
      //	RooBreitWigner SiBW_10("sigBW_10","sigBW_10" , *mgg, meanBW, sigmaBW_10);
      RooGenericPdf SiBW_10("sigBW_10","1 / ( TMath::Power( TMath::Power(@0,2) - TMath::Power(@1,2) , 2 ) + TMath::Power(@1,2)*TMath::Power(@2,2) )", RooArgSet(*mgg, meanBW, sigmaBW_10));    
      RooFFTConvPdf  ConvolutedRes_10("conv10", "conv10", *mgg,SiBW_10, ResAddPdf);
      
      RooRealVar var_15("var_w15", "var_w15",15);
      var_15.setConstant();	
      RooFormulaVar sigmaBW_15("w15", "w15","@0",  var_15);     
      //	RooBreitWigner SiBW_15("sigBW_15","sigBW_15" , *mgg, meanBW, sigmaBW_15);
      RooGenericPdf SiBW_15("sigBW_15","1 / ( TMath::Power( TMath::Power(@0,2) - TMath::Power(@1,2) , 2 ) + TMath::Power(@1,2)*TMath::Power(@2,2) )", RooArgSet(*mgg, meanBW, sigmaBW_15));    
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
      std::cout<<meanBW->getVal()<<"   --------"<<std::endl;
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
      
      c1->SaveAs("plots/SignalModels_differentWidths.png");
      c1->SaveAs("plots/SignalModels_differentWidths.pdf");
      c1->SaveAs("~/www/plotsNota/SignalModels_differentWidths.pdf");
      c1->SaveAs("~/www/plotsNota/SignalModels_differentWidths.png");
    }

    // IMPORTANT: fix all pdf parameters to constant
    w->defineSet(TString::Format("ConvolutedPdfParam_cat%d",c),RooArgSet( *w->var(TString::Format("ReducedMass_sig_sigma_cat%d",c)), 
									  *w->var(TString::Format("ReducedMass_sig_alphaCBpos_cat%d",c)),
									  *w->var(TString::Format("ReducedMass_sig_alphaCBneg_cat%d",c)),
									  *w->var(TString::Format("ReducedMass_sig_Npos_cat%d",c)),
									  *w->var(TString::Format("ReducedMass_sig_Nneg_cat%d",c)),	   
									  *w->var(TString::Format("ReducedMass_sig_frac_cat%d",c)),  
									  *w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)),
									  //	  *w->var(TString::Format("v1_cat%d",c)),
									  *w->var(TString::Format("sigmaBW_var_cat%d",c))));
    
    SetConstantParams(w->set(TString::Format("ConvolutedPdfParam_cat%d",c)));
  }
}

void SigModelFitCBC(RooWorkspace* w, Float_t mass) {

  Int_t ncat = NCAT;
  RooDataSet* sigToFit[NCAT];
  
  
  Float_t MASS(mass);  
  Float_t minMassFit(mass*0.8);
  Float_t maxMassFit(mass*1.2); 
  if(mass==150.) minMassFit = MINmass; 
  std::cout<<"----------------------------------------------------------------------------------------"<<std::endl;

  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
 
  // Fit Signal 
  for (int c=0; c<ncat; ++c) {

    cout << "---------- Category = " << c << endl;
    sigToFit[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c));
   
    RooRealVar* mgg = w->var("mgg"); 

    //cb
    RooFormulaVar CBpos_mean_draw(TString::Format("CBpos_mean_draw_cat%d",c),"","@0+250",*w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)));
    RooFormulaVar CBpos_sigma(TString::Format("CBpos_sigma_cat%d",c),"","@0*250",*w->var(TString::Format("ReducedMass_sig_sigma_cat%d",c) ));
    RooFormulaVar CBpos_alphaCB(TString::Format("CBpos_alphaCB_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_alphaCBpos_cat%d",c)) );
    RooFormulaVar CBneg_alphaCB(TString::Format("CBneg_alphaCB_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_alphaCBneg_cat%d",c)) );
    RooFormulaVar CBpos_n(TString::Format("CBpos_n_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_Npos_cat%d",c)) );
    RooFormulaVar CBneg_n(TString::Format("CBneg_n_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_Nneg_cat%d",c)) );
    RooFormulaVar CBpos_frac(TString::Format("CBpos_frac_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_frac_cat%d",c)) );
    
    RooCBShape ResCBpos_draw(TString::Format("ResCBpos_draw_cat%d",c),TString::Format("ResCBpos_draw_cat%d",c) , *mgg, CBpos_mean_draw, CBpos_sigma, CBpos_alphaCB, CBpos_n) ;
    RooCBShape ResCBneg_draw(TString::Format("ResCBneg_draw_cat%d",c),TString::Format("ResCBneg_draw_cat%d",c) , *mgg, CBpos_mean_draw, CBpos_sigma, CBneg_alphaCB, CBneg_n) ;

    //add CB neg + Gauss
    RooFormulaVar Gauss_frac(TString::Format("Gauss_frac_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_fracGauss_cat%d",c)));    
    RooFormulaVar Gauss_sigma(TString::Format("Gauss_sigma_cat%d",c),"","@0*250",*w->var(TString::Format("ReducedMass_sig_sigmaGauss_cat%d",c)));
    RooGaussian ResGauss_draw(TString::Format("ResGauss_draw_cat%d",c),TString::Format("ResGauss_draw_cat%d",c),*mgg, CBpos_mean_draw, Gauss_sigma );
    RooAddPdf ResAddGaussPdf_draw(TString::Format("mggSig_cat%d",c),TString::Format("mggSig_cat%d",c), RooArgList(ResCBneg_draw, ResGauss_draw), Gauss_frac);
    
    //CBC
    RooFormulaVar CBC_mean(TString::Format("CBC_mean_cat%d",c),"","@0",*w->var(TString::Format("mgg_sig_mean_cat%d",c)) );
    RooFormulaVar CBC_sigma(TString::Format("CBC_sigma_cat%d",c),"","@0",*w->var(TString::Format("mgg_sig_sigma_cat%d",c)) );
    RooFormulaVar CBC_alphaC(TString::Format("CBC_alphaC_cat%d",c),"","@0",*w->var(TString::Format("mgg_sig_alphaC_cat%d",c)) );
    RooFormulaVar CBC_alphaCB(TString::Format("CBC_alphaCB_cat%d",c),"","@0",*w->var(TString::Format("mgg_sig_alphaCB_cat%d",c)) );
    RooFormulaVar CBC_n(TString::Format("CBC_n_cat%d",c),"","@0",*w->var(TString::Format("mgg_sig_n_cat%d",c)) );

    
    RooCBCrujffPdf ResCBCPdf_draw(TString::Format("mggSig_cat%d",c),TString::Format("mggSig_cat%d",c) , *mgg, CBC_mean, CBC_sigma, CBC_alphaC, CBC_alphaCB, CBC_n) ; 
    w->import(ResCBCPdf_draw);


    double width=0.1;
    if(width < 2.){ //if i want to plot the fit

      RooFitResult* fitresults_Gauss = (RooFitResult* ) ResCBCPdf_draw.fitTo(*sigToFit[c],Range(minMassFit,maxMassFit),SumW2Error(kTRUE), RooFit::Save(kTRUE));
      std::cout<<TString::Format("******************************** Signal Fit results Gauss mass %f cat %d***********************************", mass, c)<<std::endl;
      fitresults_Gauss->Print("V");
      
      RooPlot* plotOnlyResPdf = w->var("mgg")->frame(Range(minMassFit-20,maxMassFit+20),Bins(100));
      sigToFit[c]->plotOn(plotOnlyResPdf, LineColor(kRed), LineStyle(kDashed));
      double max = plotOnlyResPdf->GetMaximum();
      
      // Plot to verify everything is ok
      RooPlot* plotmggAll = w->var("mgg")->frame(Range(minMassFit-20,maxMassFit+20),Bins(100));
      sigToFit[c]->plotOn(plotmggAll);

      ResCBCPdf_draw.plotOn(plotmggAll, LineColor(kRed), LineStyle(kDashed));

      TCanvas* c1 = new TCanvas("c1","mgg",0,0,800,800);
      c1->cd(1);
      
      plotmggAll->Draw();  
      
      TLegend *legmc = new TLegend(0.5491457,0.75,0.801457,0.9340659, TString::Format("Category %d",c), "brNDC");
      legmc->AddEntry(plotmggAll->getObject(0),"Simulation","LPE");
      //legmc->AddEntry(plotmggAll->getObject(2),"CB + Gauss","L");
      //legmc->AddEntry(plotmggAll->getObject(3),"Gauss","L");
      //legmc->AddEntry(plotmggAll->getObject(4),"CB","L");
      legmc->SetTextSize(0.0206044);
      legmc->SetTextFont(42);
      legmc->SetBorderSize(0);
      legmc->SetFillStyle(0);
      legmc->Draw();
      
      label_cms->Draw("same");
      label_sqrt->Draw("same");
      
      int massI(mass);
      c1->SetLogy();
      plotmggAll->GetYaxis()->SetRangeUser(0.0000000001,max*10. );
      c1->SaveAs("plots/prelimSignalCBGaussCconvBW"+TString::Format("_M%d_cat%d_LOG.png",massI,c));
      c1->SaveAs("plots/prelimSignalCBGaussCconvBW"+TString::Format("_M%d_cat%d_LOG.root",massI,c));
    }
    
    // IMPORTANT: fix all pdf parameters to constant
    w->defineSet(TString::Format("CBCParam_cat%d",c),RooArgSet(  *w->var(TString::Format("mgg_sig_alphaCB_cat%d",c)),
								 *w->var(TString::Format("mgg_sig_n_cat%d",c)),	   
								 *w->var(TString::Format("mgg_sig_alphaC_cat%d",c)),  
								 *w->var(TString::Format("mgg_sig_mean_cat%d",c)),  
								 *w->var(TString::Format("mgg_sig_sigma_cat%d",c))));  
    
    
    SetConstantParams(w->set(TString::Format("CBCParam_cat%d",c)));
  }
}

// Background model
RooFitResult* BkgModelFitExpPARFunc(RooWorkspace* w, Bool_t dobands, Float_t mass,Int_t c,  bool blind) {
  
  Int_t ncat = NCAT;
  
  RooDataSet* data;
 
  RooFitResult* fitresult;

  RooPlot* plotmggBkg;

  Float_t minMassFit, maxMassFit;
  minMassFit = MINmass;
  maxMassFit = MAXmass;
  
  // Fit data with background pdf for data limit
  RooRealVar* mgg = w->var("mgg");  
  mgg->setUnit("GeV");
 
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  
  data = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
  mgg->setRange("bkg range", MINmass, MAXmass);
    
  // fit con expol 
  RooFormulaVar *p1mod= new RooFormulaVar(TString::Format("par1ExpPAR_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_ExpPAR1_cat%d",c)));
  RooFormulaVar *p2mod= new RooFormulaVar(TString::Format("par2ExpPAR_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_ExpPAR2_cat%d",c)));
  RooAbsPdf* mggBkg = new RooGenericPdf(TString::Format("mggBkg_cat%d",c), "exp(-@1*@0)*pow(@0, @2)", RooArgList(*mgg, *p1mod, *p2mod));
    
  fitresult = mggBkg->fitTo(*data, Range(minMassFit,maxMassFit),RooFit::FitOptions("MHTER"), SumW2Error(kTRUE), Save(kTRUE));   
   
  std::cout<<TString::Format("******************************** Background Fit results mass %f cat %d ***********************************", mass, c)<<std::endl;
  fitresult->Print("V");
  w->import(*mggBkg);
 

  //************************************************
  // Plot mgg background fit results per categories 
  TCanvas* ctmp = new TCanvas("ctmp","mgg Background ",0,0,800,800);
  Int_t nBinsMass(60);
  plotmggBkg = mgg->frame(minMassFit, maxMassFit,nBinsMass);
  
  data->plotOn(plotmggBkg,RooFit::Invisible());    
  
  mggBkg->plotOn(plotmggBkg,LineColor(kBlue),Range(minMassFit,maxMassFit),NormRange("bkg range")); 
  
  double chi2 = plotmggBkg->chiSquare(2);
  Int_t ndof = nBinsMass-2;
  std::cout<<"------> "<< ndof<<std::endl;
  double prob = TMath::Prob(chi2*ndof, ndof);
  std::cout<<prob<<std::endl;
  
  blind=false;
  if( blind ) {
    RooDataSet* data_down = (RooDataSet*) data.reduce(*w->var("mgg"),"mgg < 327.666 ");     // chiara
    RooDataSet* data_up = (RooDataSet*) data.reduce(*w->var("mgg"),"mgg >650");             // chiara
    data_up->plotOn(plotmggBkg);    
    data_down->plotOn(plotmggBkg); 
  } else {
    data->plotOn(plotmggBkg,XErrorSize(0));    
  } 

  TH1F* h_data = new TH1F("h_data","h_data", 60,minMassFit, maxMassFit);
  TH1F* h_pdf = new TH1F("h_pdf","h_pdf", 60,minMassFit, maxMassFit);
  h_data = (TH1F*) data->createHistogram("mgg", 60, 0, 0);
  h_pdf =  (TH1F*) mggBkg->createHistogram("mgg",60);
  h_pdf ->Scale(h_data->Integral()/h_pdf->Integral());

  //-------pad 1-------//
  //TPad * pad1 = new TPad("pad1", "pad1",0.01256281,0.1304945,0.5741206,1);  
  TPad * pad1 = new TPad("pad1", "pad1",0., 0., 1., 1.);  
  pad1->SetLogy();
  pad1->Draw();
  pad1->cd();
  pad1->Range(154.1111,-5.730293,650.7778,4.501062);
  pad1->SetLeftMargin(0.1789709);
  pad1->SetRightMargin(0.01565995);
  pad1->SetTopMargin(0.04897314);
  pad1->SetBottomMargin(0.1691167);
  
  plotmggBkg->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
  plotmggBkg->GetYaxis()->SetTitle("Events/6.7 GeV");     // chiara
  plotmggBkg->SetAxisRange(0.001,plotmggBkg->GetMaximum()*1.5,"Y");
  plotmggBkg->Draw();  
  
  TLegend *legdata = new TLegend(0.5334677,0.680339,0.8245968,0.8958475, TString::Format("Class %d",c), "brNDC");
  if(mass<450.)legdata = new TLegend(0.2334677,0.300339,0.5645968,0.4958475, TString::Format("Class %d",c), "brNDC");
  legdata->AddEntry(plotmggBkg->getObject(2),"Data","LP");
  legdata->AddEntry(plotmggBkg->getObject(1),"Fit Model","L");
  
  dobands=false;
  //********************************************************************************/
  /*    if (dobands) {
	RooAbsPdf *cpdf; cpdf = mggBkg;
	TGraphAsymmErrors onesigma;
	TGraphAsymmErrors twosigma;
	
	RooRealVar *nlim = new RooRealVar(TString::Format("nlim%d",c),"",0.0,0.0,10.0);
	nlim->removeRange();
	
	RooCurve *nomcurve = dynamic_cast<RooCurve*>(plotmggBkg->getObject(1));
	
	double el1;
	double eh1;
	double el2;
	double eh2;
	
	int j = 0;
	for (int i=1; i<(plotmggBkg->GetXaxis()->GetNbins()+1); ++i) {
	double lowedge = plotmggBkg->GetXaxis()->GetBinLowEdge(i);
	double upedge  = plotmggBkg->GetXaxis()->GetBinUpEdge(i);
	double center  = plotmggBkg->GetXaxis()->GetBinCenter(i);
	
	double nombkg = nomcurve->interpolate(center);
	
	
	nlim->setVal(nombkg);
	mgg->setRange("errRange",lowedge,upedge);
	RooAbsPdf *epdf = 0;
	epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange");
	
	RooAbsReal *nll = epdf->createNLL(*(data),Extended());
	RooMinimizer minim(*nll);
	minim.setStrategy(0);
	double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
	double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
	
	minim.migrad();
	minim.minos(*nlim);
	printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());


	el1 = nlim->getErrorLo();
	eh1= nlim->getErrorHi();
	//	std::cout<<"-----------------------------------------------------------------> "<< minim.minos(*nlim)<<std::endl;
	minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2)); // the 0.5 is because qmu is -2*NLL
	minim.migrad();
	minim.minos(*nlim);
	el2 = nlim->getErrorLo();
	eh2= nlim->getErrorHi();


	delete nll;
	delete epdf;
	if( el1 != 0. && eh1 != 0. && el2 != 0. && eh2 != 0. &&  el1 != 1. && eh1 != 1.  && el2 != 1.  && eh2 != 1. ) {
	onesigma.SetPoint(j,center,nombkg);
	twosigma.SetPoint(j,center,nombkg);
	onesigma.SetPointError(j,0.,0.,-el1,eh1);
	twosigma.SetPointError(j,0.,0.,-el2,eh2);
	j++;
	}

	}

      mgg->setRange("errRange",minMassFit,maxMassFit);
      
      twosigma.SetLineColor(kGreen);
      twosigma.SetFillColor(kGreen);
      twosigma.SetMarkerColor(kGreen);
      twosigma.Draw("C3 SAME");
  
      onesigma.SetLineColor(kYellow);
      onesigma.SetFillColor(kYellow);
      onesigma.SetMarkerColor(kYellow);
      onesigma.Draw("C3 SAME");
   
      legdata->AddEntry(&onesigma, "#pm 1 #sigma", "F" );
      legdata->AddEntry(&twosigma, "#pm 2 #sigma","F" );
      plotmggBkg->Draw("SAME"); 
     
      }
  */

  legdata->SetTextSize(0.035);
  legdata->SetTextFont(42);
  legdata->SetBorderSize(0);
  legdata->SetFillStyle(0);
  legdata->Draw("same");

  TPaveText* label_cms = get_labelCMS(0, "2012", false);
  TPaveText* label_sqrt = get_labelSqrt(0);
  //    label_cms->Draw("same");
  //label_sqrt->Draw("same");
  int iPos=11 ;
  // CMS_lumi( pad1,false,iPos );
  
  // write down the chi2 of the fit on the
  TPaveText* label_chi2= new TPaveText(0.5744355,0.6050847,0.80583871,0.65822034,"brNDC");
  if(mass<450. || c >1)label_chi2 = new TPaveText(0.4479505,0.1861729,0.6095114,0.2933837,"brNDC");
  label_chi2->SetFillColor(kWhite);
  label_chi2->SetTextSize(0.035);
  label_chi2->SetTextFont(42);
  label_chi2->SetTextAlign(31); // align right
  label_chi2->AddText(TString::Format("Fit chi square/dof = %.3f", chi2));
  label_chi2->AddText(TString::Format("Chi square Prob = %.3f", prob));
  label_chi2->Draw("same");
  
  ctmp->cd();

  //-------pad 2------//
  TPad * pad2 = new TPad("pad2", "pad2",0.01507538,0.01510989,0.5854271,0.260989);
  pad2->SetGrid();
  pad2->Range(157.9171,-7.736842,659.5746,3.568421);
  pad2->SetLeftMargin(0.1696035);
  pad2->SetRightMargin(0.03303965);
  pad2->SetTopMargin(0.05027933);
  pad2->SetBottomMargin(0.4189944);
  pad2->cd();

  h_data->Sumw2();
  h_pdf->Sumw2();
  TH1F* pull = new TH1F("pull", "pull", 20., -3., 3.);
  TH1F* h_data1 = h_data->Clone();
  h_data->Add(h_pdf, -1);
  double ndata;
  double nfit;
  double errdata;
  double errfit;
  
  for(int i = 1; i< h_data->GetNbinsX()+1;i++){
    ndata=h_data1->GetBinContent(i);
    errdata=h_data1->GetBinError(i);
    nfit=h_pdf->GetBinContent(i);
    errfit=h_pdf->GetBinError(i);
    h_data->SetBinContent(i,h_data->GetBinContent(i)/errfit);
    pull->Fill(h_data->GetBinContent(i));
  }
  
  Double_t mean = pull->GetMean();
  Double_t sigma = pull->GetRMS();
  Double_t meanErr =pull->GetMeanError();
  Double_t sigmaErr = pull->GetRMSError();
  
  std::cout<<"MEAN: "<<pull->GetMean()<<"  RMS: "<<pull->GetRMS()<<std::endl;
  h_data->GetYaxis()->SetRangeUser(-3., 3.);
  h_data->GetYaxis()->SetNdivisions(505);
  h_data->SetMarkerSize(0.4);
  h_data->GetXaxis()->SetTitle("m_{#gamma#gamma}");
  h_data->GetXaxis()->SetTitleSize(0.15);
  h_data->GetXaxis()->SetLabelSize(0.1);
  h_data->GetYaxis()->SetLabelSize(0.1);
  h_data->GetYaxis()->SetTitleSize(0.13);
  h_data->GetYaxis()->SetTitle("#frac{N_{obs}-N_{exp}}{#sqrt{N_{exp}}}");
  h_data->GetYaxis()->SetTitleOffset(0.45);
  h_data->GetXaxis()->SetTitleOffset(0.8);
  for(int i=0; i<h_data->GetNbinsX()+1;i++) h_data->SetBinError(i, 0.);
  h_data->Draw("P");
  TH1F* h_dataCopy = h_data->Clone();
  for (int i = 0; i< h_dataCopy->GetNbinsX()+1;i++) if (h_dataCopy->GetBinContent(i)==0) h_dataCopy->SetBinContent(i, 1);
  
  h_dataCopy->Add(h_dataCopy, -1);
  for (int i = 0; i< h_dataCopy->GetNbinsX()+1;i++) h_dataCopy->SetBinError(i, 1);
  h_dataCopy->SetLineColor(kAzure-2);
  h_dataCopy->SetFillColor(kAzure-2);
  h_dataCopy->SetFillStyle(3002);
  h_dataCopy->SetMarkerSize(0.);
  h_dataCopy->Draw("HISTE23same");
  h_data->Draw("PSAME");
  
  //-------pad 3------//
  ctmp->cd();
  TPad * pad3 = new TPad("pad3", "pad3",0.5703518,0.01785714,0.7738693,0.2554945);
  pad3->SetGrid();
  pad3->cd();
  pad3->Range(-1.052224,-7.652275,12.89804,3.257984);
  pad3->SetLeftMargin(0.07542683);
  pad3->SetRightMargin(0.1718992);
  pad3->SetTopMargin(0.02364605);
  pad3->SetBottomMargin(0.4264129);
  
  pad3->SetFrameFillStyle(0);
  pad3->SetFrameBorderMode(0);
  pad3->SetFrameFillStyle(0);
  pad3->SetFrameBorderMode(0);
  
  pull->SetFillColor(kBlack);
  pull->Draw("hbar");
  
  ctmp->cd();

  //TPaveText* label2 = new TPaveText(0.5745968,0.1779661,0.7681452,0.3919492, "brNDC" );
  //label2->SetFillColor(kWhite);
  //label2->SetBorderSize(0.);
  //label2->SetTextSize(0.0213922);
  //label2->SetTextAlign(11);
  //label2->SetTextFont(42);
  //label2->AddText(TString::Format("#splitline{Mean: %.3f #pm %.3f }{RMS: %.3f #pm %.3f}", mean,meanErr,sigma, sigmaErr));
  // label2->Draw("same");
  //,   pad3->Draw();
  int massI(mass);
  ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPPAR_M%d.png",c,massI));
  ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPPAR_M%d.pdf",c,massI));
  
  ctmp->SetLogy();
  if(mass<650)plotmggBkg->SetAxisRange(0.5,plotmggBkg->GetMaximum()*10,"Y");
  else plotmggBkg->SetAxisRange(0.5,100,"Y");
  ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPPAR_LOG_M%d.png",c,massI));
  ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPPAR_LOG_M%d.pdf",c,massI));
  if(mass==350){
    ctmp->SaveAs("~/www/plotsNota/prelimBkg"+TString::Format("_cat%d_EXPPAR_LOG_M%d.png",c,massI));
    ctmp->SaveAs("~/www/plotsNota/prelimBkg"+TString::Format("_cat%d_EXPPAR_LOG_M%d.pdf",c,massI));
    ctmp->SaveAs("~/www/plotsPAS/prelimBkg"+TString::Format("_cat%d_EXPPAR_LOG_M%d.png",c,massI));
    ctmp->SaveAs("~/www/plotsPAS/prelimBkg"+TString::Format("_cat%d_EXPPAR_LOG_M%d.pdf",c,massI));
    ctmp->SaveAs("~/www/plotsPAS/prelimBkg"+TString::Format("_cat%d_EXPPAR_LOG_M%d.C",c,massI));
    ctmp->SaveAs("~/www/plotsPAS/prelimBkg"+TString::Format("_cat%d_EXPPAR_LOG_M%d.root",c,massI));
  }
  
  RooFitResult* r;
  
  return r;
}

RooFitResult* BkgModelFitDiJetFunc(RooWorkspace* w, Bool_t dobands, Float_t mass, bool blind) {

  Int_t ncat = NCAT;
  std::cout<<"isBlind: "<<blind<<std::endl;

  // retrieve pdfs and datasets from workspace to fit with pdf models
  RooDataSet* data[NCAT];
 
  RooFitResult* fitresult[NCAT];;

  RooPlot* plotmggBkg[NCAT];

  RooAbsPdf* mggSig[NCAT];
  
  Float_t minMassFit, maxMassFit;
  minMassFit = MINmass;
  maxMassFit = MAXmass;
   
  // Fit data with background pdf for data limit
  RooRealVar* mgg = w->var("mgg");  
  mgg->setUnit("GeV");
 
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  
  for (int c = 0; c < ncat; ++c) {
    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
  
    // fit a la dijets
    RooFormulaVar *p0mod = new RooFormulaVar(TString::Format("par0DiJet_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_norm_cat%d",c)));
    RooFormulaVar *p1mod = new RooFormulaVar(TString::Format("par1DiJet_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_slope1_cat%d",c)));
    RooFormulaVar *p2mod = new RooFormulaVar(TString::Format("par2DiJet_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_slope2_cat%d",c)));
    RooFormulaVar *p3mod = new RooFormulaVar(TString::Format("par3DiJet_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_slope3_cat%d",c)));
    mgg->setRange("bkg range", MINmass,MAXmass);
    RooFormulaVar *x     = new RooFormulaVar(TString::Format("xDiJet_cat%d",c),"","@0/8000.",*w->var("mgg"));
    RooAbsPdf* mggBkgTmp0 = new RooGenericPdf(TString::Format("mggBkg_cat%d",c), "pow(1-@0, @2)/pow(@0, @1+@3*log(@0))", RooArgList(*x, *p1mod, *p2mod,*p3mod));
   
    fitresult[c] = mggBkgTmp0->fitTo(*data[c], Range(minMassFit,maxMassFit),RooFit::FitOptions("MHTER"), SumW2Error(kTRUE), Save(kTRUE));
    w->import(*mggBkgTmp0);
   
    std::cout<<TString::Format("******************************** Background Fit results mass %f cat %d***********************************", mass, c)<<std::endl;
    fitresult[c]->Print("V");
    
    
    //************************************************
    // Plot mgg background fit results per categories 
    TCanvas* ctmp = new TCanvas("ctmp","mgg Background Categories",0,0,500,500);
    Int_t nBinsMass(60);
    plotmggBkg[c] = mgg->frame(minMassFit, maxMassFit,nBinsMass);
  
    data[c]->plotOn(plotmggBkg[c],RooFit::Invisible());    
   
    mggBkgTmp0->plotOn(plotmggBkg[c],LineColor(kBlue),Range(minMassFit,maxMassFit),NormRange("bkg range")); 
    double chi2 = plotmggBkg[c]->chiSquare(3);
    Int_t ndof = nBinsMass-3;
    std::cout<<"------> "<< ndof<<std::endl;
    double prob = TMath::Prob(chi2*ndof, ndof);
    std::cout<<prob<<std::endl;
    blind = true;

    if( blind ) {
      RooDataSet* data_down = (RooDataSet*) data[c].reduce(*w->var("mgg"),"mgg < 178.");   // chiara
      RooDataSet* data_up = (RooDataSet*) data[c].reduce(*w->var("mgg"),"mgg >402");       // chiara
      data_up->plotOn(plotmggBkg[c]);    
      data_down->plotOn(plotmggBkg[c]); 
    } else {
      data[c]->plotOn(plotmggBkg[c]);    
    } 
    
    plotmggBkg[c]->GetXaxis()->SetTitle("m_{#gamma #gamma}[GeV]");
    plotmggBkg[c]->GetYaxis()->SetTitle("Events/6.7 GeV");  // chiara
    plotmggBkg[c]->SetAxisRange(0.001,plotmggBkg[c]->GetMaximum()*1.5,"Y");
    plotmggBkg[c]->Draw();  
    
    TLegend *legdata = new TLegend(0.3790323,0.7775424,0.6290323,0.9279661, TString::Format("Class %d",c), "brNDC");
    legdata->AddEntry(plotmggBkg[c]->getObject(2),"Data","LPE");
    legdata->AddEntry(plotmggBkg[c]->getObject(1),"Parametric Model: DiJet","L");
    legdata->SetTextSize(0.035);
    legdata->SetTextFont(42);
    legdata->SetBorderSize(0);
    legdata->SetFillStyle(0);
    legdata->Draw("same");

    TPaveText* label_cms = get_labelCMS(0, "2012", false);
    TPaveText* label_sqrt = get_labelSqrt(0);
    label_cms->Draw("same");
    label_sqrt->Draw("same");

    // write down the chi2 of the fit on the
    TPaveText* label_chi2 = new TPaveText(0.5524194,0.6419492,0.796371,0.7690678, "brNDC");
    label_chi2->SetFillColor(kWhite);
    label_chi2->SetTextSize(0.035);
    label_chi2->SetTextFont(42);
    label_chi2->SetTextAlign(31); // align right
    label_chi2->AddText(TString::Format("Fit chi square/dof = %.3f", chi2));
    label_chi2->AddText(TString::Format("Chi square Prob = %.3f", prob));
    label_chi2->Draw("same");

    //********************************************************************************
    if (dobands) {

      RooAbsPdf *cpdf; cpdf = mggBkgTmp0;
      TGraphAsymmErrors *onesigma = new TGraphAsymmErrors();
      TGraphAsymmErrors *twosigma = new TGraphAsymmErrors();
      
      RooRealVar *nlim = new RooRealVar(TString::Format("nlim%d",c),"",0.0,0.0,10.0);
      nlim->removeRange();
      
      RooCurve *nomcurve = dynamic_cast<RooCurve*>(plotmggBkg[c]->getObject(1));
      
      for (int i=1; i<(plotmggBkg[c]->GetXaxis()->GetNbins()+1); ++i) {
	double lowedge = plotmggBkg[c]->GetXaxis()->GetBinLowEdge(i);
	double upedge  = plotmggBkg[c]->GetXaxis()->GetBinUpEdge(i);
	double center  = plotmggBkg[c]->GetXaxis()->GetBinCenter(i);
	
	double nombkg = nomcurve->interpolate(center);
	nlim->setVal(nombkg);
	mgg->setRange("errRange",lowedge,upedge);
	RooAbsPdf *epdf = 0;
	epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange");
	
	RooAbsReal *nll = epdf->createNLL(*(data[c]),Extended());
	RooMinimizer minim(*nll);
	minim.setStrategy(0);
	double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
	double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
	
	minim.migrad();
	minim.minos(*nlim);
	// printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());
	
	onesigma->SetPoint(i-1,center,nombkg);
	onesigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	
	minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2)); // the 0.5 is because qmu is -2*NLL
	// eventually if cl = 0.95 this is the usual 1.92!      
	
	minim.migrad();
	minim.minos(*nlim);
	
	twosigma->SetPoint(i-1,center,nombkg);
	twosigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	
	delete nll;
	delete epdf;
      }

      mgg->setRange("errRange",minMassFit,maxMassFit);
      
      twosigma->SetLineColor(kGreen);
      twosigma->SetFillColor(kGreen);
      twosigma->SetMarkerColor(kGreen);
      twosigma->Draw("L3 SAME");
      
      onesigma->SetLineColor(kYellow);
      onesigma->SetFillColor(kYellow);
      onesigma->SetMarkerColor(kYellow);
      onesigma->Draw("L3 SAME");
      
      plotmggBkg[c]->Draw("SAME"); 
    }

    int massI(mass);
    ctmp->SaveAs("plots/prelimBkg"+TString::Format("_cat%d_DIJET_M%d.png",c,massI));
    ctmp->SaveAs("plots/prelimBkg"+TString::Format("_cat%d_DIJET_M%d.pdf",c,massI));
    ctmp->SaveAs("plots/prelimBkg"+TString::Format("_cat%d_DIJET_M%d.root",c,massI));

    ctmp->SetLogy();
    plotmggBkg[c]->SetAxisRange(1.3,plotmggBkg[c]->GetMaximum()*1.5,"Y");
    ctmp->SaveAs("plots/prelimBkg"+TString::Format("_cat%d_DIJET_LOG_M%d.png",c,massI));
    ctmp->SaveAs("plots/prelimBkg"+TString::Format("_cat%d_DIJET_LOG_M%d.pdf",c,massI));
    ctmp->SaveAs("plots/prelimBkg"+TString::Format("_cat%d_DIJET_LOG_M%d.root",c,massI));

  }
  return fitresult;
}


/*
RooFitResult* BkgModelFitDiJetEXPFunc(RooWorkspace* w, Bool_t dobands, Float_t mass,Int_t c, bool blind) {

  Int_t ncat = NCAT;
  std::cout<<"isBlind: "<<blind<<std::endl;
  // retrieve pdfs and datasets from workspace to fit with pdf models
  RooDataSet* data;
 
  RooFitResult* fitresult;;
  
  RooPlot* plotmggBkg;

  // dobands and dosignal
  RooDataSet* signal;

  RooAbsPdf* mggSig;
  
  Float_t minMassFit, maxMassFit;

    minMassFit = MINmass;
    maxMassFit = MAXmass;

  
  // Fit data with background pdf for data limit
  RooRealVar* mgg = w->var("mgg");  
  mgg->setUnit("GeV");
 
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  
 
    data = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    
  
    // fit a la dijets
    
    RooFormulaVar *p1mod = new RooFormulaVar(TString::Format("par1DiJetEXP_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_slope1_3_cat%d",c)));
    RooFormulaVar *p2mod = new RooFormulaVar(TString::Format("par2DiJetEXP_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_slope2_3_cat%d",c)));
    RooFormulaVar *p3mod = new RooFormulaVar(TString::Format("par3DiJetEXP_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_slope3_3_cat%d",c)));
    RooFormulaVar *exp1 = new RooFormulaVar(TString::Format("expDiJetEXP_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_exp1DiJetEXP_cat%d",c)));
    RooFormulaVar *pFrac1 = new RooFormulaVar(TString::Format("fracDiJetEXP_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_fracDiJetEXP_cat%d",c)));
   
   
    RooGenericPdf* mggBkgTmp0DiJet = new RooGenericPdf(TString::Format("mggBkg_DIJETE_truth_cat%d",c), "pow(1-@0/8000., @2)/pow(@0/8000., @1+@3*log(@0/8000.))", RooArgList(*w->var("mgg"), *p1mod, *p2mod,*p3mod));
  
    RooExponential* mggBkgTmp0Exp = new RooExponential(TString::Format("mggBkg_EXP_truth_cat%d",c),"", *w->var("mgg"),  *exp1);
    
   
    RooAddPdf* mggBkgTmpAdd = new RooAddPdf(TString::Format("mggBkg_cat%d",c),TString::Format("mggBkg_cat%d",c) , RooArgList(*mggBkgTmp0DiJet, *mggBkgTmp0Exp), RooArgList(*pFrac1));
    
    fitresult = mggBkgTmpAdd->fitTo(*data,RooFit::FitOptions("MHTR"), Save(kTRUE));//RooFit::FitOptions("MHTER"), Range(minMassFit,maxMassFit),    
    w->import(*mggBkgTmpAdd);
  
    std::cout<<TString::Format("******************************** Background DiJetEXP Fit results mass %f cat %d ***********************************", mass, c)<<std::endl;
    fitresult->Print("V");
   

    //************************************************
    // Plot mgg background fit results per categories 
    TCanvas* ctmp = new TCanvas("ctmp","mgg Background Categories",0,0,500,500);
    Int_t nBinsMass(60);
    plotmggBkg = mgg->frame(minMassFit, maxMassFit,nBinsMass);
  
    data->plotOn(plotmggBkg,RooFit::Invisible());    
   
    mggBkgTmpAdd->plotOn(plotmggBkg,LineColor(kBlue)); 
    mggBkgTmpAdd->plotOn(plotmggBkg,Components(TString::Format("mggBkg_EXP_truth_cat%d",c)),LineColor(kViolet),LineStyle(kDashed)); 
    mggBkgTmpAdd->plotOn(plotmggBkg,Components(TString::Format("mggBkg_DIJETE_truth_cat%d",c)),LineColor(kOrange),LineStyle(kDashed));   

    double chi2 = plotmggBkg->chiSquare(3);
    Int_t ndof = nBinsMass-5;
    std::cout<<"------> "<< ndof<<std::endl;
    double prob = TMath::Prob(chi2*ndof, ndof);
    std::cout<<prob<<std::endl;
   
    if( blind ) {
 
      RooDataSet* data_down = (RooDataSet*) data.reduce(*w->var("mgg"),"mgg < 173.5");
      RooDataSet* data_up = (RooDataSet*) data.reduce(*w->var("mgg")," mgg >850");
      TH1F* h_up= new TH1F("h_up", "h_up",nBinsMass, 130, 1000);
      h_up->Sumw2();
      data_up->fillHistogram(h_up, RooArgList(*mgg));
      TH1F* h_down= new TH1F("h_down", "h_down",nBinsMass, 130, 1000);
      h_down->Sumw2();
      data_down->fillHistogram(h_down, RooArgList(*mgg));
   
    } else {
      data->plotOn(plotmggBkg);    
      } 
       
   
    plotmggBkg->GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
    plotmggBkg->SetAxisRange(0.1,10000,"Y");
    plotmggBkg->Draw();  
   if( blind ) {
       h_up->Draw("sameP");
       h_down->Draw("sameP");
     }
  
    TLegend *legdata = new TLegend(0.3790323,0.7775424,0.6290323,0.9279661, TString::Format("Category %d",c), "brNDC");
    legdata->AddEntry(plotmggBkg->getObject(4),"Data","LPE");
    legdata->AddEntry(plotmggBkg->getObject(1),"Parametric Model: DiJetEXP","L");  
    legdata->AddEntry(plotmggBkg->getObject(2),"Parametric Model: EXP","L");
    legdata->AddEntry(plotmggBkg->getObject(3),"Parametric Model: DiJet","L");
  
    legdata->SetTextSize(0.035);
    legdata->SetTextFont(42);
    // legdata->SetTextAlign(31);
    legdata->SetBorderSize(0);
    legdata->SetFillStyle(0);
    legdata->Draw("same");

    TPaveText* label_cms = get_labelCMS(0, "2012", false);
    TPaveText* label_sqrt = get_labelSqrt(0);
    label_cms->Draw("same");
    label_sqrt->Draw("same");

    //write down the chi2 of the fit on the
 
    TPaveText* label_chi2 = new TPaveText(0.5524194,0.6419492,0.796371,0.7690678, "brNDC");
    label_chi2->SetFillColor(kWhite);
    label_chi2->SetTextSize(0.035);
    label_chi2->SetTextFont(42);
    label_chi2->SetTextAlign(31); // align right
    label_chi2->AddText(TString::Format("Fit chi square/dof = %.3f", chi2));
    label_chi2->AddText(TString::Format("Chi square Prob = %.3f", prob));
    label_chi2->Draw("same");

    

    // ********************************************************************************
    if (dobands) {

      RooAbsPdf *cpdf; cpdf = mggBkgTmp0;
      TGraphAsymmErrors *onesigma = new TGraphAsymmErrors();
      TGraphAsymmErrors *twosigma = new TGraphAsymmErrors();
      
      RooRealVar *nlim = new RooRealVar(TString::Format("nlim%d",c),"",0.0,0.0,10.0);
      nlim->removeRange();
      
      RooCurve *nomcurve = dynamic_cast<RooCurve*>(plotmggBkg[c]->getObject(1));
      
      for (int i=1; i<(plotmggBkg->GetXaxis()->GetNbins()+1); ++i) {
	double lowedge = plotmggBkg->GetXaxis()->GetBinLowEdge(i);
	double upedge  = plotmggBkg->GetXaxis()->GetBinUpEdge(i);
	double center  = plotmggBkg->GetXaxis()->GetBinCenter(i);
	
	double nombkg = nomcurve->interpolate(center);
	nlim->setVal(nombkg);
	mgg->setRange("errRange",lowedge,upedge);
	RooAbsPdf *epdf = 0;
	epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange");
	
	RooAbsReal *nll = epdf->createNLL(*(data),Extended());
	RooMinimizer minim(*nll);
	minim.setStrategy(0);
	double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
	double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
	
	minim.migrad();
	minim.minos(*nlim);
	// printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());
	
	onesigma->SetPoint(i-1,center,nombkg);
	onesigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	
	minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2)); // the 0.5 is because qmu is -2*NLL
	// eventually if cl = 0.95 this is the usual 1.92!      
	
	minim.migrad();
	minim.minos(*nlim);
	
	twosigma->SetPoint(i-1,center,nombkg);
	twosigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	
	delete nll;
	delete epdf;
      }

      mgg->setRange("errRange",minMassFit,maxMassFit);
      
      twosigma->SetLineColor(kGreen);
      twosigma->SetFillColor(kGreen);
      twosigma->SetMarkerColor(kGreen);
      twosigma->Draw("L3 SAME");
      
      onesigma->SetLineColor(kYellow);
      onesigma->SetFillColor(kYellow);
      onesigma->SetMarkerColor(kYellow);
      onesigma->Draw("L3 SAME");
      
      plotmggBkg->Draw("SAME"); 
    }

    int massI(mass);
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXP_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXP_M%d.pdf",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXP_M%d.root",c,massI));

    ctmp->SetLogy();
    //  plotmggBkg->SetAxisRange(1.3,plotmggBkg[c]->GetMaximum()*1.5,"Y");
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXP_LOG_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXP_LOG_M%d.pdf",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXP_LOG_M%d.root",c,massI));

 

  return fitresult;
}
*/

/*
RooFitResult* BkgModelFitDiJetEXPOLFunc(RooWorkspace* w, Bool_t dobands, Float_t mass, Int_t c, bool blind) {

  Int_t ncat = NCAT;
 
  // retrieve pdfs and datasets from workspace to fit with pdf models
  RooDataSet* data;
 
  RooFitResult* fitresult;;

  RooPlot* plotmggBkg;

  // dobands and dosignal
  RooDataSet* signal;

  RooAbsPdf* mggSig;
  
  Float_t minMassFit, maxMassFit;

    minMassFit = MINmass;
    maxMassFit = MAXmass;

  
  // Fit data with background pdf for data limit
  RooRealVar* mgg = w->var("mgg");  
  mgg->setUnit("GeV");
 
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  
  
    data = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    
  
    // fit a la dijets
    
    RooFormulaVar *p1mod = new RooFormulaVar(TString::Format("par1DiJetEXPOL_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_slope1_4_cat%d",c)));
    RooFormulaVar *p2mod = new RooFormulaVar(TString::Format("par2DiJetEXPOL_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_slope2_4_cat%d",c)));
    RooFormulaVar *p3mod = new RooFormulaVar(TString::Format("par3DiJetEXPOL_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_slope3_4_cat%d",c)));
    RooFormulaVar *expol1 = new RooFormulaVar(TString::Format("expol1DiJetEXPOL_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_expol1_4_cat%d",c)));
    RooFormulaVar *expol2 = new RooFormulaVar(TString::Format("expol2DiJetEXPOL_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_expol2_4_cat%d",c)));
    RooFormulaVar *pFrac1 = new RooFormulaVar(TString::Format("fracDiJetEXPOL_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_fracDiJetEXPOL_cat%d",c)));
   
   
    RooGenericPdf* mggBkgTmp0DiJet = new RooGenericPdf(TString::Format("mggBkg_DiJetEx_cat%d",c), "pow(1-@0/8000., @2)/pow(@0/8000., @1+@3*log(@0/8000.))", RooArgList(*w->var("mgg"), *p1mod, *p2mod,*p3mod));
    RooGenericPdf* mggBkgTmp0Expol = new RooGenericPdf(TString::Format("mggBkg_ExpolDiJ_cat%d",c), "exp(-@0/(@1+@2*@0))", RooArgList(*w->var("mgg"), *expol1, *expol2));
    
   
    RooAddPdf* mggBkgTmpAdd = new RooAddPdf(TString::Format("mggBkg_cat%d",c),TString::Format("mggBkg_cat%d",c) , RooArgList(*mggBkgTmp0DiJet, *mggBkgTmp0Expol), RooArgList(*pFrac1));
    
    fitresult = mggBkgTmpAdd->fitTo(*data,RooFit::FitOptions("MHTR"), Save(kTRUE));//RooFit::FitOptions("MHTER"), Range(minMassFit,maxMassFit),    
    w->import(*mggBkgTmpAdd);
  
    std::cout<<TString::Format("******************************** Background DiJetEXPOL Fit results mass %f cat %d ***********************************", mass, c)<<std::endl;
    fitresult->Print("V");
   

    //************************************************
    // Plot mgg background fit results per categories 
    TCanvas* ctmp = new TCanvas("ctmp","mgg Background Categories",0,0,500,500);
    Int_t nBinsMass(60);
    plotmggBkg = mgg->frame(minMassFit, maxMassFit,nBinsMass);
  
    data->plotOn(plotmggBkg,RooFit::Invisible());    
   
    mggBkgTmpAdd->plotOn(plotmggBkg,LineColor(kBlue)); 
    mggBkgTmpAdd->plotOn(plotmggBkg,Components(TString::Format("mggBkg_ExpolDiJ_cat%d",c)),LineColor(kViolet),LineStyle(kDashed)); 
    mggBkgTmpAdd->plotOn(plotmggBkg,Components(TString::Format("mggBkg_DiJetEx_cat%d",c)),LineColor(kOrange),LineStyle(kDashed));   

    double chi2 = plotmggBkg->chiSquare(3);
    Int_t ndof = nBinsMass-6;
    std::cout<<"------> "<< ndof<<std::endl;
    double prob = TMath::Prob(chi2*ndof, ndof);
    std::cout<<prob<<std::endl;
    blind = false;
    if( blind ) {
     
      RooDataSet* data_down = (RooDataSet*) data.reduce(*w->var("mgg"),"mgg < 178.");
      RooDataSet* data_up = (RooDataSet*) data.reduce(*w->var("mgg"),"mgg >402");

      data_up->plotOn(plotmggBkg);    
      data_down->plotOn(plotmggBkg); 


   
    } else {
      data->plotOn(plotmggBkg);    
      } 
       
   
    plotmggBkg->GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
    plotmggBkg->SetAxisRange(0.001,plotmggBkg->GetMaximum()*1.5,"Y");
   
    plotmggBkg->Draw();  

    TLegend *legdata = new TLegend(0.3790323,0.7775424,0.6290323,0.9279661, TString::Format("Category %d",c), "brNDC");
    legdata->AddEntry(plotmggBkg->getObject(4),"Data","LPE");
    legdata->AddEntry(plotmggBkg->getObject(1),"Parametric Model: DiJetEXPOL","L");  
    legdata->AddEntry(plotmggBkg->getObject(2),"Parametric Model: Expol","L");
    legdata->AddEntry(plotmggBkg->getObject(3),"Parametric Model: DiJet","L");
  
    legdata->SetTextSize(0.035);
    legdata->SetTextFont(42);
    // legdata->SetTextAlign(31);
    legdata->SetBorderSize(0);
    legdata->SetFillStyle(0);
    legdata->Draw("same");

    TPaveText* label_cms = get_labelCMS(0, "2012", false);
    TPaveText* label_sqrt = get_labelSqrt(0);
    label_cms->Draw("same");
    label_sqrt->Draw("same");

    //write down the chi2 of the fit on the
 
    TPaveText* label_chi2 = new TPaveText(0.5524194,0.6419492,0.796371,0.7690678, "brNDC");
    label_chi2->SetFillColor(kWhite);
    label_chi2->SetTextSize(0.035);
    label_chi2->SetTextFont(42);
    label_chi2->SetTextAlign(31); // align right
    label_chi2->AddText(TString::Format("Fit chi square/dof = %.3f", chi2));
    label_chi2->AddText(TString::Format("Chi square Prob = %.3f", prob));
    label_chi2->Draw("same");

    

    //*******************************************************************************
    if (dobands) {

      RooAbsPdf *cpdf; cpdf = mggBkgTmp0;
      TGraphAsymmErrors *onesigma = new TGraphAsymmErrors();
      TGraphAsymmErrors *twosigma = new TGraphAsymmErrors();
      
      RooRealVar *nlim = new RooRealVar(TString::Format("nlim%d",c),"",0.0,0.0,10.0);
      nlim->removeRange();
      
      RooCurve *nomcurve = dynamic_cast<RooCurve*>(plotmggBkg->getObject(1));
      
      for (int i=1; i<(plotmggBkg->GetXaxis()->GetNbins()+1); ++i) {
	double lowedge = plotmggBkg->GetXaxis()->GetBinLowEdge(i);
	double upedge  = plotmggBkg->GetXaxis()->GetBinUpEdge(i);
	double center  = plotmggBkg->GetXaxis()->GetBinCenter(i);
	
	double nombkg = nomcurve->interpolate(center);
	nlim->setVal(nombkg);
	mgg->setRange("errRange",lowedge,upedge);
	RooAbsPdf *epdf = 0;
	epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange");
	
	RooAbsReal *nll = epdf->createNLL(*(data),Extended());
	RooMinimizer minim(*nll);
	minim.setStrategy(0);
	double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
	double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
	
	minim.migrad();
	minim.minos(*nlim);
	// printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());
	
	onesigma->SetPoint(i-1,center,nombkg);
	onesigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	
	minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2)); // the 0.5 is because qmu is -2*NLL
	// eventually if cl = 0.95 this is the usual 1.92!      
	
	minim.migrad();
	minim.minos(*nlim);
	
	twosigma->SetPoint(i-1,center,nombkg);
	twosigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	
	delete nll;
	delete epdf;
      }

      mgg->setRange("errRange",minMassFit,maxMassFit);
      
      twosigma->SetLineColor(kGreen);
      twosigma->SetFillColor(kGreen);
      twosigma->SetMarkerColor(kGreen);
      twosigma->Draw("L3 SAME");
      
      onesigma->SetLineColor(kYellow);
      onesigma->SetFillColor(kYellow);
      onesigma->SetMarkerColor(kYellow);
      onesigma->Draw("L3 SAME");
      
      plotmggBkg->Draw("SAME"); 
    }

    int massI(mass);
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXPOL_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXPOL_M%d.pdf",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXPOL_M%d.root",c,massI));

    ctmp->SetLogy();
    plotmggBkg->SetAxisRange(1.3,plotmggBkg->GetMaximum()*1.5,"Y");
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXPOL_LOG_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXPOL_LOG_M%d.pdf",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXPOL_LOG_M%d.root",c,massI));

  



  return fitresult;
}

*/

/*
RooFitResult* BkgModelFitExpolFunc(RooWorkspace* w, Bool_t dobands, Float_t mass, bool blind) {

  Int_t ncat = NCAT;
  std::cout<<"isBlind: "<<blind<<std::endl;
  // retrieve pdfs and datasets from workspace to fit with pdf models
  RooDataSet* data[NCAT];
 
  RooFitResult* fitresult[NCAT];;
  RooPlot* plotmggBkg[NCAT];

  // dobands and dosignal
  RooDataSet* signal[NCAT];

  RooAbsPdf* mggSig[NCAT];
  
  Float_t minMassFit, maxMassFit;
 
    minMassFit = MINmass;
    maxMassFit = MAXmass;
  
  // Fit data with background pdf for data limit
  RooRealVar* mgg = w->var("mgg");  
  mgg->setUnit("GeV");
 
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  
  for (int c = 0; c < ncat; ++c) {
    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    
  
    // fit con expo pol
   
    RooFormulaVar *p1mod = new RooFormulaVar(TString::Format("par1Expol_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_expol1_cat%d",c)));
    RooFormulaVar *p2mod = new RooFormulaVar(TString::Format("par2Expol_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_expol2_cat%d",c)));
    RooFormulaVar *p3mod = new RooFormulaVar(TString::Format("par3Expol_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_expol3_cat%d",c)));
 
    mgg->setRange("bkg range", MINmass,MAXmass);
    RooFormulaVar *x     = new RooFormulaVar(TString::Format("xExpol_cat%d",c),"","@0",*w->var("mgg"));

   
    RooAbsPdf* mggBkgTmp0 = new RooGenericPdf(TString::Format("mggBkg_cat%d",c), "exp(-@0*@0/(@1+@2*@0+@3*@0*@0))", RooArgList(*x, *p1mod, *p2mod,*p3mod));
   

    fitresult[c] = mggBkgTmp0->fitTo(*data[c], Range(minMassFit,maxMassFit),RooFit::FitOptions("MHTER"), SumW2Error(kTRUE), Save(kTRUE));
    w->import(*mggBkgTmp0);
   
    std::cout<<TString::Format("******************************** Background Fit results mass %f cat %d***********************************", mass, c)<<std::endl;
    fitresult[c]->Print("V");
   

    //************************************************
    // Plot mgg background fit results per categories 
    TCanvas* ctmp = new TCanvas("ctmp","mgg Background Categories",0,0,500,500);
    Int_t nBinsMass(60);
    plotmggBkg[c] = mgg->frame(minMassFit, maxMassFit,nBinsMass);
  
    data[c]->plotOn(plotmggBkg[c],RooFit::Invisible());    
   
    mggBkgTmp0->plotOn(plotmggBkg[c],LineColor(kBlue),Range(minMassFit,maxMassFit),NormRange("bkg range")); 
    double chi2 = plotmggBkg[c]->chiSquare(3);
    Int_t ndof = nBinsMass-2;
    std::cout<<"------> "<< ndof<<std::endl;
    double prob = TMath::Prob(chi2*ndof, ndof);
    std::cout<<prob<<std::endl;
    blind =false;
    if( blind ) {
     
      RooDataSet* data_down = (RooDataSet*) data[c].reduce(*w->var("mgg"),"mgg < 178.");
      RooDataSet* data_up = (RooDataSet*) data[c].reduce(*w->var("mgg"),"mgg >402");

      data_up->plotOn(plotmggBkg[c]);    
      data_down->plotOn(plotmggBkg[c]); 


   
    } else {
      data[c]->plotOn(plotmggBkg[c]);    
      } 
       
  
    plotmggBkg[c]->GetXaxis()->SetTitle("m_{#gamma #gamma}[GeV]");
    plotmggBkg[c]->SetAxisRange(0.001,plotmggBkg[c]->GetMaximum()*1.5,"Y");
    plotmggBkg[c]->Draw();  

    TLegend *legdata = new TLegend(0.3790323,0.7775424,0.6290323,0.9279661, TString::Format("Category %d",c), "brNDC");
    legdata->AddEntry(plotmggBkg[c]->getObject(2),"Data","LPE");
    legdata->AddEntry(plotmggBkg[c]->getObject(1),"Parametric Model: Expol","L");
    legdata->SetTextSize(0.035);
    legdata->SetTextFont(42);
    // legdata->SetTextAlign(31);
    legdata->SetBorderSize(0);
    legdata->SetFillStyle(0);
    legdata->Draw("same");

    TPaveText* label_cms = get_labelCMS(0, "2012", false);
    TPaveText* label_sqrt = get_labelSqrt(0);
    label_cms->Draw("same");
    label_sqrt->Draw("same");

    //write down the chi2 of the fit on the
      
    TPaveText* label_chi2 = new TPaveText(0.5524194,0.6419492,0.796371,0.7690678, "brNDC");
    label_chi2->SetFillColor(kWhite);
    label_chi2->SetTextSize(0.035);
    label_chi2->SetTextFont(42);
    label_chi2->SetTextAlign(31); // align right
    label_chi2->AddText(TString::Format("Fit chi square/dof = %.3f", chi2));
    label_chi2->AddText(TString::Format("Chi square Prob = %.3f", prob));
    label_chi2->Draw("same");

    
    dobands = false;
    //*******************************************************************************
    if (dobands) {

      RooAbsPdf *cpdf; cpdf = mggBkgTmp0;
      TGraphAsymmErrors *onesigma = new TGraphAsymmErrors();
      TGraphAsymmErrors *twosigma = new TGraphAsymmErrors();
      
      RooRealVar *nlim = new RooRealVar(TString::Format("nlim%d",c),"",0.0,0.0,10.0);
      nlim->removeRange();
      
      RooCurve *nomcurve = dynamic_cast<RooCurve*>(plotmggBkg[c]->getObject(1));
      
      for (int i=1; i<(plotmggBkg[c]->GetXaxis()->GetNbins()+1); ++i) {
	double lowedge = plotmggBkg[c]->GetXaxis()->GetBinLowEdge(i);
	double upedge  = plotmggBkg[c]->GetXaxis()->GetBinUpEdge(i);
	double center  = plotmggBkg[c]->GetXaxis()->GetBinCenter(i);
	
	double nombkg = nomcurve->interpolate(center);
	nlim->setVal(nombkg);
	mgg->setRange("errRange",lowedge,upedge);
	RooAbsPdf *epdf = 0;
	epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange");
	
	RooAbsReal *nll = epdf->createNLL(*(data[c]),Extended());
	RooMinimizer minim(*nll);
	minim.setStrategy(0);
	double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
	double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
	
	minim.migrad();
	minim.minos(*nlim);
	// printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());
	
	onesigma->SetPoint(i-1,center,nombkg);
	onesigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	
	minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2)); // the 0.5 is because qmu is -2*NLL
	// eventually if cl = 0.95 this is the usual 1.92!      
	
	minim.migrad();
	minim.minos(*nlim);
	
	twosigma->SetPoint(i-1,center,nombkg);
	twosigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	
	delete nll;
	delete epdf;
      }

      mgg->setRange("errRange",minMassFit,maxMassFit);
      
      twosigma->SetLineColor(kGreen);
      twosigma->SetFillColor(kGreen);
      twosigma->SetMarkerColor(kGreen);
      twosigma->Draw("L3 SAME");
      
      onesigma->SetLineColor(kYellow);
      onesigma->SetFillColor(kYellow);
      onesigma->SetMarkerColor(kYellow);
      onesigma->Draw("L3 SAME");
      
      plotmggBkg[c]->Draw("SAME"); 
    }
    int massI(mass);
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPOL_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPOL_M%d.pdf",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPOL_M%d.root",c,massI));

    ctmp->SetLogy();
    plotmggBkg[c]->SetAxisRange(1.3,plotmggBkg[c]->GetMaximum()*1.5,"Y");
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPOL_LOG_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPOL_LOG_M%d.pdf",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPOL_LOG_M%d.root",c,massI));

  }



  return fitresult;
}
*/

void SetConstantParams(const RooArgSet* params) {
  
  cout << endl; cout << "Entering SetConstantParams" << endl;
  TIterator* iter(params->createIterator());
  for (TObject *a = iter->Next(); a != 0; a = iter->Next()) {
    RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
    if (rrv) { rrv->setConstant(true); std::cout << " " << rrv->GetName(); }
  }  
}


// ---------- chiara: qui ----------------

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
  wAll->import(*w->var("mggTrue"));
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

// Write background pdfs and datasets into the workspace 
void MakeBkgWS(RooWorkspace* w, const char* fileBaseName, double mass) {

  TString wsDir = "workspaces/"+filePOSTfix;
  Int_t ncat = NCAT;  

  RooWorkspace *wAll = new RooWorkspace("w_all","w_all");
  
  //********************************//
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






Double_t effSigma(TH1 *hist) {
  
  TAxis *xaxis = hist->GetXaxis();
  Int_t nb = xaxis->GetNbins();
  if(nb < 10) {
    std::cout << "effsigma: Not a valid histo. nbins = " << nb << std::endl;
    return 0.;
  }
  
  Double_t bwid = xaxis->GetBinWidth(1);
  if(bwid == 0) {
    std::cout << "effsigma: Not a valid histo. bwid = " << bwid << std::endl;
    return 0.;
  }
  Double_t xmax = xaxis->GetXmax();
  Double_t xmin = xaxis->GetXmin();
  Double_t ave = hist->GetMean();
  Double_t rms = hist->GetRMS();

  Double_t total=0.;
  for(Int_t i=0; i<nb+2; i++) {
    total+=hist->GetBinContent(i);
  }
  if(total < 100.) {
    std::cout << "effsigma: Too few entries " << total << std::endl;
    return 0.;
  }
  Int_t ierr=0;
  Int_t ismin=999;

  Double_t rlim=0.683*total;
  Int_t nrms=rms/(bwid);    // Set scan size to +/- rms
  if(nrms > nb/10) nrms=nb/10; // Could be tuned...

  Double_t widmin=9999999.;
  for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
    Int_t ibm=(ave-xmin)/bwid+1+iscan;
    Double_t x=(ibm-0.5)*bwid+xmin;
    Double_t xj=x;
    Double_t xk=x;
    Int_t jbm=ibm;
    Int_t kbm=ibm;
    Double_t bin=hist->GetBinContent(ibm);
    total=bin;
    for(Int_t j=1;j<nb;j++){
      if(jbm < nb) {
        jbm++;
        xj+=bwid;
        bin=hist->GetBinContent(jbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
      if(kbm > 0) {
        kbm--;
        xk-=bwid;
        bin=hist->GetBinContent(kbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
    }
    Double_t dxf=(total-rlim)*bwid/bin;
    Double_t wid=(xj-xk+bwid-dxf)*0.5;
    if(wid < widmin) {
      widmin=wid;
      ismin=iscan;
    }
  }
  if(ismin == nrms || ismin == -nrms) ierr=3;
  if(ierr != 0) std::cout << "effsigma: Error of type " << ierr << std::endl;

  return widmin;
}

// preparing datacards
void MakeDataCard_1Channel(RooWorkspace* w, const char* fileBaseName, const char* fileBkgName, Float_t width,int iChan, std::string model) {

  TString cardDir = "/afs/cern.ch/work/s/soffi/CMSSW611-Analysis/src/h2gglobe/ChiaraFitLimits/datacardWithAllSyst/"+filePOSTfix;
  Int_t ncat = NCAT;
  TString wsDir   = "/afs/cern.ch/work/s/soffi/CMSSW611-Analysis/src/h2gglobe/ChiaraFitLimits/workspaces/"+filePOSTfix;

  // **********************
  // Retrieve the datasets
  cout << "Start retrieving dataset" << endl;
  
  RooDataSet* data[9];
  RooDataSet* signal[9];
  for (int c=0; c<ncat; ++c) {
    data[c]   = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    signal[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c));
  }

  RooRealVar*  lumi = w->var("lumi");
  /*
  // *****************************
  // Print Expected event yields
  cout << "======== Expected Events Number =====================" << endl;  
  cout << ".........Measured Data for L = " << lumi->getVal() << " pb-1 ............................" << endl;  
  cout << "#Events data: " <<  w->data("Data")->sumEntries()  << endl;
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events data cat%d:   ",c) << data[c]->sumEntries()  << endl;
  }
  cout << ".........Expected Signal for L = " << lumi->getVal() << " pb-1 ............................" << endl;  
  cout << "#Events Signal:      " << w->data("SigWeight")->sumEntries()  << endl;
  Float_t siglikeErr[6];
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events Signal cat%d: ",c) << signal[c]->sumEntries() << endl;
    siglikeErr[c]=0.6*signal[c]->sumEntries();
  }
  cout << "====================================================" << endl;  

 */
  // *****************************
  // Printdata Data Card int file
  TString filename(cardDir+TString(fileBaseName)+"_"+"8TeV"+Form(("_m%.2f_w%.2f_channel%d_"+model+".txt").c_str(),w->var("MH")->getVal(),width,iChan));
  ofstream outFile(filename);

  outFile << "#CMS-HGG HighMass DataCard for Unbinned Limit Setting, " << lumi->getVal() <<  " pb-1 " << endl;
  outFile << "#Run with: combine -d datacardName.txt -U -m *mass* -H ProfileLikelihood -M MarkovChainMC --rMin=0 --rMax=20.0  -b 3000 -i 50000 --optimizeSim=1 --tries 30" << endl;
  outFile << "# Lumi =  " << lumi->getVal() << " pb-1" << endl;
  outFile << "imax *" << endl;
  outFile << "jmax *" << endl;
  outFile << "kmax *" << endl;
  outFile << "---------------" << endl;

  outFile << "shapes data_obs * " << wsDir+TString(fileBkgName)+TString::Format("_m%.2f.root",w->var("MH")->getVal()) << " w_all:data_obs_$CHANNEL" << endl;
  outFile << "shapes sig * "      << wsDir+TString(fileBaseName)+"_8TeV"+TString::Format(("_m%.2f_w%.2f.inputsig_"+model+".root").c_str(),w->var("MH")->getVal(),width) << " w_all:mggSig_$CHANNEL" << endl;
  outFile << "shapes bkg * "      << wsDir+TString(fileBkgName)+TString::Format("_m%.2f.root",w->var("MH")->getVal()) << " w_all:mggBkg_$CHANNEL" << endl;

  outFile << "---------------" << endl;
  outFile << Form("bin          cat%d", iChan) << endl;
  outFile <<  "observation   "  <<  Form("%.10lg",data[iChan]->sumEntries()) << endl;
  outFile << "------------------------------" << endl;
  outFile << "bin                 " << Form("cat%d      cat%d", iChan, iChan) << endl;
  outFile << "process                 sig      bkg" << endl;
  outFile << "process                   0        1" << endl;
  // if(signalScaler==1.)
  // signalScaler=1./signal[2]->sumEntries()*20;
  outFile << "rate                   " 
    //	  << signal[iChan]->sumEntries()*signalScaler << " " << data[iChan]->sumEntries() << endl;
    // << 1 << " " << data[iChan]->sumEntries() << endl;
	  <<19620. << "  "<<data[iChan]->sumEntries() << endl;
  outFile << "--------------------------------" << endl;
  outFile << "# signal scaled by " << signalScaler << endl;

 
  outFile << "lumi_8TeV     lnN     1.026000  - " << endl;
  outFile << "eff_trig     lnN     1.010000  - " << endl;
  outFile << "global_syst     lnN     1.050000  - " << endl;
  if(width==0.1) outFile << "bw_syst     lnN     1.00100  - " << endl;
  if(width==2) outFile << "bw_syst     lnN     1.0200  - " << endl;
  if(width==5) outFile << "bw_syst     lnN     1.0500  - " << endl;
  if(width==7) outFile << "bw_syst     lnN     1.0700  - " << endl;
  if(width==10) outFile << "bw_syst     lnN     1.1000  - " << endl;
  if(iChan==0){    
    outFile << "id_eff_eb     lnN     1.005000  - " << endl;    
    outFile << "id_eff_ee     lnN     1.000000  - " << endl;    
    outFile << "r9Eff   lnN   1.0145/0.9915   - " << endl;
    outFile << "vtxEff   lnN   0.996/1.008   - " << endl; 
  }else if(iChan==1){    
    outFile << "id_eff_eb     lnN     1.005000  - " << endl;    
    outFile << "id_eff_ee     lnN     1.000000  - " << endl;    
    outFile << "r9Eff   lnN   0.985/1.0085   - " << endl;
    outFile << "vtxEff   lnN   0.998/1.005   - " << endl; 
  }else if(iChan==2){    
    outFile << "id_eff_eb     lnN     1.00000 - " << endl;    
    outFile << "id_eff_ee     lnN     1.02600  - " << endl;    
    outFile << "r9Eff   lnN   1.0193/0.964   - " << endl;
    outFile << "vtxEff   lnN   0.996/1.007   - " << endl; 
  }else if(iChan==3){    
    outFile << "id_eff_eb     lnN     1.000000  - " << endl;    
    outFile << "id_eff_ee     lnN     1.026000  - " << endl;    
    outFile << "r9Eff   lnN   0.981/1.0357   - " << endl;
    outFile << "vtxEff   lnN   0.998/1.003   - " << endl; 
  }
    
    outFile << Form("mShift_cat%d    param   0 1 ",iChan) << endl;
    outFile << Form("mSmear_cat%d    param   0 1 ",iChan) << endl;



  // outFile << "CMS_VV_eff_g         lnN  0.8/1.20      - # Signal Efficiency" << endl;
  // outFile << "# Parametric shape uncertainties, entered by hand." << endl;
  // outFile << Form("CMS_hgg_sig_m0_absShift    param   1   0.0125   # displacement of the mean w.r.t. nominal in EB*EX category, good R9",iChan) << endl;
  // outFile << Form("CMS_hgg_sig_sigmaScale     param   1   0.1   # multiplicative correction to sigmas in EB*EX category, good R9",iChan) << endl;
  // outFile << Form("rooHistFunc_cat%d_norm       ",iChan) << endl;
  // outFile << Form("CMS_hgg_bkg_8TeV_slope2_cat%d         flatParam  # Mean and absolute uncertainty on background slope",iChan) << endl;
  // outFile << Form("CMS_hgg_bkg_8TeV_slope3_cat%d         flatParam  # Mean and absolute uncertainty on background slope",iChan) << endl;
  // if (iChan != 2 )  outFile << Form("CMS_hgg_bkg_8TeV_slope1_cat%d         flatParam  # Mean and absolute uncertainty on background slope",iChan) << endl;
  outFile.close();

  cout << "Write data card in: " << filename << " file" << endl;

  return;
}





Double_t effSigma(TH1 *hist) {

  TAxis *xaxis = hist->GetXaxis();
  Int_t nb = xaxis->GetNbins();
  if(nb < 10) {
    std::cout << "effsigma: Not a valid histo. nbins = " << nb << std::endl;
    return 0.;
  }

  Double_t bwid = xaxis->GetBinWidth(1);
  if(bwid == 0) {
    std::cout << "effsigma: Not a valid histo. bwid = " << bwid << std::endl;
    return 0.;
  }
  Double_t xmax = xaxis->GetXmax();
  Double_t xmin = xaxis->GetXmin();
  Double_t ave = hist->GetMean();
  Double_t rms = hist->GetRMS();

  Double_t total=0.;
  for(Int_t i=0; i<nb+2; i++) {
    total+=hist->GetBinContent(i);
  }
  if(total < 100.) {
    std::cout << "effsigma: Too few entries " << total << std::endl;
    return 0.;
  }
  Int_t ierr=0;
  Int_t ismin=999;

  Double_t rlim=0.683*total;
  Int_t nrms=rms/(bwid);    // Set scan size to +/- rms
  if(nrms > nb/10) nrms=nb/10; // Could be tuned...

  Double_t widmin=9999999.;
  for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
    Int_t ibm=(ave-xmin)/bwid+1+iscan;
    Double_t x=(ibm-0.5)*bwid+xmin;
    Double_t xj=x;
    Double_t xk=x;
    Int_t jbm=ibm;
    Int_t kbm=ibm;
    Double_t bin=hist->GetBinContent(ibm);
    total=bin;
    for(Int_t j=1;j<nb;j++){
      if(jbm < nb) {
        jbm++;
        xj+=bwid;
        bin=hist->GetBinContent(jbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
      if(kbm > 0) {
        kbm--;
        xk-=bwid;
        bin=hist->GetBinContent(kbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
    }
    Double_t dxf=(total-rlim)*bwid/bin;
    Double_t wid=(xj-xk+bwid-dxf)*0.5;
    if(wid < widmin) {
      widmin=wid;
      ismin=iscan;
    }
  }
  if(ismin == nrms || ismin == -nrms) ierr=3;
  if(ierr != 0) std::cout << "effsigma: Error of type " << ierr << std::endl;

  return widmin;
}

RooHistFunc* getRooHistFunc(int cat, RooRealVar* var, std::string model){


  double mass[7] = {150., 200., 250., 300., 400., 600, 800};

  double c0[7];
  double c1[7];
  double c2[7];
  double c3[7];
  double all[7];
  //changed 10/04
  /*  double c0g[5] = {0.00300125 , 0.00377226 , 0.00446634 , 0.0051454 , 0.00592569};
  double c1g[5] = {0.00416221 , 0.00503098 , 0.00494583 , 0.00483801 , 0.00465511 };
  double c2g[5] = {0.0015748 , 0.00188656 , 0.00215004 , 0.00246588 , 0.00260045 };
  double c3g[5] = {0.00209722 , 0.00255168 , 0.0025782 , 0.00251383 , 0.00218195 };
  double allg[5] = {0.0108564 , 0.0132374 , 0.0141401 , 0.014937 , 0.0153788 };  
  double c0h[5] = {0.00302 , 0.003772 , 0.0044645 , 0.005133 , 0.00595};
  double c1h[5] = {0.00414 , 0.004974 , 0.0049136 , 0.0048243 , 0.00463 };
  double c2h[5] = {0.00158 , 0.001885 , 0.0021557 , 0.0024568 , 0.0026 };
  double c3h[5] = {0.002092 , 0.002533 , 0.00256584 , 0.0024928 , 0.002173 };
  double allh[5] = {0.010853 , 0.01316 , 0.0140993 , 0.014881 , 0.0153686};*/

  double c0g[7] = {0.00281259889626 ,0.00360405499841 ,0.0042512557194 ,0.00480861314285 ,0.00555793387016 ,0.00711359319501 ,0.00823700633031};
  double c1g[7] = { 0.00384737104156 ,0.00409835498596 ,0.00396547370455 ,0.003904278683 ,0.00368952356883 ,0.00312694674604 ,0.00264070618382   };
  double c2g[7] = {0.00147564281183 ,0.00190380886793 ,0.00200106257462 ,0.00237975631151 ,0.00252602251139 ,0.00252633118101 ,0.00243976040669 };
  double c3g[7] = {0.00179971742008 ,0.001998562755 ,0.00209254398016 ,0.00195718051131 ,0.00166752834056 ,0.00124744484509 ,0.000929175591328};
  double allg[7] = { 0.00993533016973 ,0.0116047816073 ,0.0123103359787 ,0.0130498286487 ,0.0134410082909 ,0.0140143159671 ,0.0142466485122};  

  double c0v[7] = {0.00295597,   0.00351098,   0.00404037,   0.00454414,   0.00547479,0.00702859,   0.00817238};
  double c1v[7] = {0.00402382,   0.00410415,   0.00415866,   0.00418736,   0.00416732,0.00381745,   0.00305454   };
  double c2v[7] = {0.00144884,   0.00165231,   0.00183088,   0.00198455,   0.0022172,0.00238368 ,  0.00215175 };
  double c3v[7] = {0.00205624,   0.00204125,   0.00201627,   0.00198131,   0.00188142,0.0015618 ,  0.00108238 };
  double allv[7] = {0.0104849,   0.0113087,   0.0120462,   0.0126974   0.0137407,0.0147915,   0.0144611   };
  
  for(int i = 0; i<7; i++){
    if (model=="GGH"){
      c0[i] = c0g[i];
      c1[i] = c1g[i];
      c2[i] = c2g[i];
      c3[i] = c3g[i];
      all[i] = allg[i];

    }else if(model=="VBF"){
      c0[i] = c0v[i];
      c1[i] = c1v[i];
      c2[i] = c2v[i];
      c3[i] = c3v[i];
      all[i] = allv[i];


    }
  }
  
  TH1F* h_all = new TH1F("h_all", "h_all", 15, 130., 875.);
  //m150-400
  for(int i=0;i<5;i++){
    std::cout<<"i: "<<i<<" mass: "<<mass[i]<<std::endl;
    if(cat==0) h_all->SetBinContent(h_all->FindBin(mass[i]),c0[i]);
    if(cat==1) h_all->SetBinContent(h_all->FindBin(mass[i]),c1[i]);
    if(cat==2) h_all->SetBinContent(h_all->FindBin(mass[i]),c2[i]);
    if(cat==3) h_all->SetBinContent(h_all->FindBin(mass[i]),c3[i]);
    if(cat==4) h_all->SetBinContent(h_all->FindBin(mass[i]),all[i]);
   
  }
  std::cout<<"--------------------------->"<<cat<<std::endl;
  //m350
  if(cat==0) h_all->SetBinContent(h_all->FindBin(350),(c0[3]+c0[4])/2);
  if(cat==1) h_all->SetBinContent(h_all->FindBin(350),(c1[3]+c1[4])/2);
  if(cat==2) h_all->SetBinContent(h_all->FindBin(350),(c2[3]+c2[4])/2);
  if(cat==3) h_all->SetBinContent(h_all->FindBin(350),(c3[3]+c3[4])/2);
  if(cat==4) h_all->SetBinContent(h_all->FindBin(350),(all[3]+all[4])/2);

  //m450-550 
  double r1mass[3]={450, 500, 550};
  for(int i = 0; i<3;i++){
  if(cat==0) h_all->SetBinContent(h_all->FindBin(r1mass[i]),(c0[4]+c0[5])/2);
  if(cat==1) h_all->SetBinContent(h_all->FindBin(r1mass[i]),(c1[4]+c1[5])/2);
  if(cat==2) h_all->SetBinContent(h_all->FindBin(r1mass[i]),(c2[4]+c2[5])/2);
  if(cat==3) h_all->SetBinContent(h_all->FindBin(r1mass[i]),(c3[4]+c3[5])/2);
  if(cat==4) h_all->SetBinContent(h_all->FindBin(r1mass[i]),(all[4]+all[5])/2);
  }

  //m600 e 850
  double highmass[2]={ 600., 850.};
  for(int i = 0; i<2;i++){
  if(cat==0) h_all->SetBinContent(h_all->FindBin(highmass[i]),c0[i+5]);
  if(cat==1) h_all->SetBinContent(h_all->FindBin(highmass[i]),c1[i+5]);
  if(cat==2) h_all->SetBinContent(h_all->FindBin(highmass[i]),c2[i+5]);
  if(cat==3) h_all->SetBinContent(h_all->FindBin(highmass[i]),c3[i+5]);
  if(cat==4) h_all->SetBinContent(h_all->FindBin(highmass[i]),all[i+5]);
  }
  

  //m650-800 
  double r2mass[4]={650, 700, 750, 800};
  for(int i = 0; i<4;i++){
  if(cat==0) h_all->SetBinContent(h_all->FindBin(r2mass[i]),(c0[5]+c0[6])/2);
  if(cat==1) h_all->SetBinContent(h_all->FindBin(r2mass[i]),(c1[5]+c1[6])/2);
  if(cat==2) h_all->SetBinContent(h_all->FindBin(r2mass[i]),(c2[5]+c2[6])/2);
  if(cat==3) h_all->SetBinContent(h_all->FindBin(r2mass[i]),(c3[5]+c3[6])/2);
  if(cat==4) h_all->SetBinContent(h_all->FindBin(r2mass[i]),(all[5]+all[6])/2);
  }

  for(int i=0;i<16;i++){
    std::cout<<cat<<std::endl;
    // std::cout<<mass[i]<<"  "<<h_all->FindBin(mass[i])<<"   "<<c0[i]<<std::endl;
    std::cout<<i<<"  "<<h_all->GetBinCenter(i)<<"  "<<h_all->GetBinContent(i)<<std::endl;

 }
  TCanvas* c = new TCanvas("c", "c", 1);
  c->cd();
  //  TPaveText* label_cms = get_labelCMS(0, false);
  //TPaveText* label_sqrt = get_labelSqrt(0);  
  h_all->GetYaxis()->SetRangeUser(0.0001, 0.1);
  h_all->GetXaxis()->SetTitle("m_{X} [GeV]");
  h_all->GetYaxis()->SetTitle("Signal Yield/1pb");
  h_all->Draw();

  RooDataHist* rooData_all = new RooDataHist("rooData_all", "roData_all",*var,h_all);
  RooHistFunc* rooFunc_all = new RooHistFunc(TString::Format("mggSig_cat%d_norm",cat),TString::Format("mggSig_cat%d_norm",cat), *var,*rooData_all, 3);
  /*  RooPlot* plot = var->frame();
  rooData_all->plotOn(plot);
  rooFunc_all->plotOn(plot);
  plot->Draw("same");
  c->SetLogy();

  c->SaveAs("plots/signalYield.png");
  c->SaveAs("plots/signalYield.pdf");*/
  
  return rooFunc_all;

}


RooHistFunc* getNorm2D(int cat,RooRealVar* var,double mass, double width, std::string model ){
  
  double norm;
  TFile* f;
  if(model=="GRAVITON") f= new TFile("effXacc2D_GRAVITON.root", "READ");
  else f= new TFile("effXacc2D.root", "READ");
  f->Print();
  TH2D* h2 = (TH2D*) f->Get(TString::Format("h2_cat%d",cat));
  int bin = h2->GetYaxis()->FindBin(width);
  std::cout<<bin<<std::endl;
 
  TH1D* h1 = (TH1D*)h2->ProjectionX("py", bin, bin);

  double eff;
  eff = h1->GetBinContent(h1->FindBin(mass));
  norm = 1.87004999e-02*eff/100.;
  for(int i = 0; i<h1->GetNbinsX()+1; i++){
    std::cout<< i<<" M: "<<h1->GetBinCenter(i)<<"   "<<h1->GetBinContent(i)<<std::endl;
    h1->SetBinContent(i,1.87004999e-02*h1->GetBinContent(i)/100);
    //std::cout<<h1->GetBinCenter(i)<<"   "<<h1->GetBinContent(i)<<std::endl;
  }
  std::cout<<norm<<std::endl;
  RooDataHist* rooData_all = new RooDataHist("rooData_all", "roData_all",*var,h1);
  RooHistFunc* rooFunc_all = new RooHistFunc(TString::Format("mggSig_cat%d_norm",cat),TString::Format("mggSig_cat%d_norm",cat), *var,*rooData_all, 3);
  // std::cout<<rooFunc_all->getVal(mass)<<std::endl;
  /* TCanvas* c = new TCanvas("c", "c", 1);
  c->cd();
  //  TPaveText* label_cms = get_labelCMS(0, false);
  //TPaveText* label_sqrt = get_labelSqrt(0);  
  // h1->GetYaxis()->SetRangeUser(0.0001, 0.1);
  h1->GetXaxis()->SetTitle("m_{X} [GeV]");
  h1->GetYaxis()->SetTitle("Signal Yield/1pb");
  h1->Draw();
  RooPlot* plot = var->frame();
  rooData_all->plotOn(plot);
  rooFunc_all->plotOn(plot);
  plot->Draw("same");
  c->SetLogy();

  c->SaveAs(TString::Format("~/www/signalYield_cat%d.png",cat));
  c->SaveAs(TString::Format("~/www/signalYield_cat%d.pdf",cat));*/
  return rooFunc_all;

}

// chiara: questi vanno aggiornati
// min of the mass fit range - values from the bias study
RooHistFunc* getRooHistFuncFitMIN(int cat, RooRealVar* var){

  double mass[3] = {1500., 3000., 5000.};
  double c0[3]   = { 500., 2000., 4000.};
  double c1[3]   = { 500., 2000., 4000.};
  double c2[3]   = { 500., 2000., 4000.};
  double c3[3]   = { 500., 2000., 4000.};

  TH1F* h_all = new TH1F("h_all", "h_all", 3, 500, 6000);
  h_all->GetXaxis()->SetTitle("m_{X} [GeV]");
  h_all->GetYaxis()->SetTitle("Fit Range: Minimum");

  for (int i=0; i<3; i++) {
    //if (cat==0) std::cout << cat << " " << mass[i] << " " << h_all->FindBin(mass[i]) << " " << c0[i] << std::endl;
    //if (cat==1) std::cout << cat << " " << mass[i] << " " << h_all->FindBin(mass[i]) << " " << c1[i] << std::endl;
    //if (cat==2) std::cout << cat << " " << mass[i] << " " << h_all->FindBin(mass[i]) << " " << c2[i] << std::endl;
    //if (cat==3) std::cout << cat << " " << mass[i] << " " << h_all->FindBin(mass[i]) << " " << c3[i] << std::endl;
    //
    if (cat==0) h_all->SetBinContent(h_all->FindBin(mass[i]),c0[i]);
    if (cat==1) h_all->SetBinContent(h_all->FindBin(mass[i]),c1[i]);
    if (cat==2) h_all->SetBinContent(h_all->FindBin(mass[i]),c2[i]);
    if (cat==3) h_all->SetBinContent(h_all->FindBin(mass[i]),c3[i]);
  }

  TCanvas* c = new TCanvas("c", "c", 1);
  c->cd();
  TPaveText* label_cms  = get_labelCMS(0,"2015", true);
  TPaveText* label_sqrt = get_labelSqrt(0);  
  h_all->Draw("");
 
  RooDataHist* rooData_all = new RooDataHist("rooData_all","roData_all",*var,h_all);
  RooHistFunc* rooFunc_all = new RooHistFunc(TString::Format("mggBkg_cat%d_min",cat),TString::Format("mggBkg_cat%d_min",cat), *var,*rooData_all, 3);
  RooPlot* plot = var->frame();
  rooData_all->plotOn(plot);
  rooFunc_all->plotOn(plot);
  plot->Draw("same");
  c->SaveAs("plots/Bkg_fitMinimum.png");
  
  return rooFunc_all;
}

// chiara: questi vanno aggiornati
// min of the mass fit range - extracted from a fit to the values from the bias study
TF1* getFuncFitMIN(int cat, RooRealVar* var){

  double mass[3] = {1500., 3000., 5000.};
  double c0[3]   = { 500., 2000., 4000.};
  double c1[3]   = { 500., 2000., 4000.};
  double c2[3]   = { 500., 2000., 4000.};
  double c3[3]   = { 500., 2000., 4000.};
 
  TF1* f = new TF1("f", "[0]+[1]*x", 500, 6000);
  f->SetParameter(0,1000 );
  f->SetParameter(1, -0.5);
  f->SetParameter(2, -1000.);
  
  TH1F* h_all = new TH1F("h_all", "h_all", 3, 500, 6000);
  h_all->GetXaxis()->SetTitle("m_{X} [GeV]");
  h_all->GetYaxis()->SetTitle("Fit Range: Minimum");

  for(int i=0;i<3;i++){
    if (cat==0) h_all->SetBinContent(h_all->FindBin(mass[i]),c0[i]);
    if (cat==1) h_all->SetBinContent(h_all->FindBin(mass[i]),c1[i]);
    if (cat==2) h_all->SetBinContent(h_all->FindBin(mass[i]),c2[i]);
    if (cat==3) h_all->SetBinContent(h_all->FindBin(mass[i]),c3[i]);
  }

  TCanvas* c = new TCanvas("c", "c", 1);
  c->cd();
  TPaveText* label_cms  = get_labelCMS(0,"2015", false);
  TPaveText* label_sqrt = get_labelSqrt(0);  
  h_all->Draw("");
  h_all->Fit("f");
  f->Draw("Lsame");
  c->SaveAs("plots/Bkg_fitMinimum_withFit.png");

  // for(int i=0;i<3;i++){
  //  if(cat==0) std::cout << "cat = " << cat <<" => mass = " << mass[i] << ", point = " << c0[i] << ", fit = " << f->Eval(mass[i]) << std::endl;
  //  if(cat==1) std::cout << "cat = " << cat <<" => mass = " << mass[i] << ", point = " << c1[i] << ", fit = " << f->Eval(mass[i]) << std::endl;
  //  if(cat==2) std::cout << "cat = " << cat <<" => mass = " << mass[i] << ", point = " << c2[i] << ", fit = " << f->Eval(mass[i]) << std::endl;
  //  if(cat==3) std::cout << "cat = " << cat <<" => mass = " << mass[i] << ", point = " << c3[i] << ", fit = " << f->Eval(mass[i]) << std::endl;
  // }
  
  return f;
}

// chiara: questi vanno aggiornati
// max of the mass fit range - values from the bias study
RooHistFunc* getRooHistFuncFitMAX(int cat, RooRealVar* var){

  double mass[3] = {1500., 3000., 5000.};
  double c0[3]   = {2500., 4000., 6000.};
  double c1[3]   = {2500., 4000., 6000.};
  double c2[3]   = {2500., 4000., 6000.};
  double c3[3]   = {2500., 4000., 6000.};

  TH1F* h_all = new TH1F("h_all", "h_all", 3, 500, 6000);
  h_all->GetXaxis()->SetTitle("m_{X} [GeV]");
  h_all->GetYaxis()->SetTitle("Fit Range: Maximum");

  for(int i=0;i<3;i++){
    if(cat==0) h_all->SetBinContent(h_all->FindBin(mass[i]),c0[i]);
    if(cat==1) h_all->SetBinContent(h_all->FindBin(mass[i]),c1[i]);
    if(cat==2) h_all->SetBinContent(h_all->FindBin(mass[i]),c2[i]);
    if(cat==3) h_all->SetBinContent(h_all->FindBin(mass[i]),c3[i]);
    // 
    // if (cat==0) std::cout << cat << " " << mass[i] << " " << h_all->FindBin(mass[i]) << " " << c0[i] << std::endl;
    // if (cat==1) std::cout << cat << " " << mass[i] << " " << h_all->FindBin(mass[i]) << " " << c1[i] << std::endl;
    // if (cat==2) std::cout << cat << " " << mass[i] << " " << h_all->FindBin(mass[i]) << " " << c2[i] << std::endl;
    // if (cat==3) std::cout << cat << " " << mass[i] << " " << h_all->FindBin(mass[i]) << " " << c3[i] << std::endl;
  }

  TCanvas* c = new TCanvas("c", "c", 1);
  c->cd();
  TPaveText* label_cms = get_labelCMS(0,"2015", false);
  TPaveText* label_sqrt = get_labelSqrt(0);  
  h_all->Draw();
  
  RooDataHist* rooData_all = new RooDataHist("rooData_all", "roData_all",*var,h_all);
  RooHistFunc* rooFunc_all = new RooHistFunc(TString::Format("mggBkg_cat%d_max",cat),TString::Format("mggBkg_cat%d_max",cat), *var,*rooData_all, 3);

  RooPlot* plot = var->frame();
  rooData_all->plotOn(plot);
  rooFunc_all->plotOn(plot);
  plot->Draw("same");
  c->SaveAs("plots/Bkg_fitMaximum.png");

  return rooFunc_all;
}

// chiara: questi vanno aggiornati
// max of the mass fit range - extracted from a fit to the values from the bias study
TF1* getFuncFitMAX(int cat, RooRealVar* var){

  double mass[3] = {1500., 3000., 5000.};
  double c0[3]   = {2500., 4000., 6000.};
  double c1[3]   = {2500., 4000., 6000.};
  double c2[3]   = {2500., 4000., 6000.};
  double c3[3]   = {2500., 4000., 6000.};  

  TF1* f = new TF1("f", "[0]+[2]*pow(x,[1])", 500, 6000);
  f->SetParameter(0,1000 );
  // f->SetParameter(1, -0.5);
  f->SetParameter(1, 0.5);
  f->SetParameter(2, -1000.);
  
  TH1F* h_all = new TH1F("h_all", "h_all", 3, 500, 6000);
  h_all->GetXaxis()->SetTitle("m_{X} [GeV]");
  h_all->GetYaxis()->SetTitle("Fit Range: Maximum");
  
  for(int i=0;i<3;i++){
    if(cat==0) h_all->SetBinContent(h_all->FindBin(mass[i]),c0[i]);
    if(cat==1) h_all->SetBinContent(h_all->FindBin(mass[i]),c1[i]);
    if(cat==2) h_all->SetBinContent(h_all->FindBin(mass[i]),c2[i]);
    if(cat==3) h_all->SetBinContent(h_all->FindBin(mass[i]),c3[i]);
  }

  TCanvas* c = new TCanvas("c", "c", 1);
  c->cd();
  TPaveText* label_cms  = get_labelCMS(0,"2015", false);
  TPaveText* label_sqrt = get_labelSqrt(0);  

  h_all->Draw();
  h_all->Fit("f");
  f->Draw("Lsame");
  c->SaveAs("plots/Bkg_fitMaximum_withFit.png");

  // for(int i=0;i<3;i++){
  // if(cat==0) std::cout << "cat = " << cat <<" => mass = " << mass[i] << ", point = " << c0[i] << ", fit = " << f->Eval(mass[i]) << std::endl;
  // if(cat==1) std::cout << "cat = " << cat <<" => mass = " << mass[i] << ", point = " << c1[i] << ", fit = " << f->Eval(mass[i]) << std::endl;
  // if(cat==2) std::cout << "cat = " << cat <<" => mass = " << mass[i] << ", point = " << c2[i] << ", fit = " << f->Eval(mass[i]) << std::endl;
  // if(cat==3) std::cout << "cat = " << cat <<" => mass = " << mass[i] << ", point = " << c3[i] << ", fit = " << f->Eval(mass[i]) << std::endl;
  // }
 
  return f;
}

// ok
// plot of the fit mass range - max and min from the fit to the values extracted from the bias study 
makeRooHistPlot(RooHistFunc* rooFitMin,RooHistFunc* rooFitMax ,TF1* fmin, TF1* fmax, RooRealVar* var){

  TCanvas* cfit = new TCanvas("cfit", "cfit", 1);
  cfit->cd();
  
  RooPlot* p = var->frame();
  rooFitMin->plotOn(p, LineColor(kAzure+9), LineStyle(kDashed));
  rooFitMax->plotOn(p, LineColor(kViolet+9));
  fmax->SetLineColor(kViolet+9);
  fmax->SetLineWidth(2);
  fmin->SetLineWidth(2);
  fmin->SetLineColor(kAzure+9);
  fmin->SetLineStyle(kDashed);
  fmax->GetYaxis()->SetRangeUser(0, 6000);
  fmax->Draw("L");
  fmin->Draw("Lsame");
  fmax->GetYaxis()->SetTitle("m_{#gamma#gamma} fit range in GeV");
  fmax->GetXaxis()->SetTitle("m_{X} [GeV]");

  TLegend* leg = new TLegend(0.52, 0.16, 0.85, 0.35, "", "brNDC");
  leg->SetTextFont(42);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(fmax, "Fit Range Maximum", "L");
  leg->AddEntry(fmin, "Fit Range Minimum", "L");
  
  TPaveText* label_cms = get_labelCMS(0, "2015", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
  // label_cms->Draw("same");
  // label_sqrt->Draw("same");

  // int iPos=11 ;
  // CMS_lumi( cfit,true,iPos );
  cfit->SaveAs("plots/bkg_fitRange.png");
}

// ok
// CMS stuffs
TPaveText* get_labelCMS( int legendQuadrant, std::string year, bool sim) {

  if( legendQuadrant!=0 && legendQuadrant!=1 && legendQuadrant!=2 && legendQuadrant!=3 ) {
    std::cout << "WARNING! Legend quadrant '" << legendQuadrant << "' not yet implemented for CMS label. Using 2." << std::endl;
    legendQuadrant = 2;
  }

  float x1, y1, x2, y2;
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

// ok
TPaveText* get_labelSqrt( int legendQuadrant ) {

  if( legendQuadrant!=0 && legendQuadrant!=1 && legendQuadrant!=2 ) {
    std::cout << "WARNING! Legend quadrant '" << legendQuadrant << "' not yet implemented for Sqrt label. Using 2." << std::endl;
    legendQuadrant = 2;
  }

  float x1, y1, x2, y2;
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
