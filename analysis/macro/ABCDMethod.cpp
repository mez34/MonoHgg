#include "ABCDMethod.hh"

ABCDMethod::ABCDMethod( SamplePairVec Samples, const Double_t inLumi, const TString outdir){

  lumi = inLumi;
  fInDir = outdir;
  fOutDir = outdir+"ABCD";

  mgg_minAB1 = 100.;
  mgg_minCD  = 110.;
  mgg_maxCD  = 130.;
  mgg_maxAB2 = 180.; 
  met_minB   = 150.;
  met_minD   = 250.;
  met_maxD   = 400.;

  // make output root file
  MakeOutDirectory(Form("%s",fOutDir.Data()));
  fOutFile = new TFile(Form("%s/analysis.root",fOutDir.Data()),"RECREATE");
  CheckValidFile(fOutFile, Form("%s/analysis.root",fOutDir.Data())); 

  // make vectors with names of samples
  for (SamplePairVecIter iter = Samples.begin(); iter != Samples.end(); ++iter){
    if ( (*iter).second == 1 ) fBkgNames.push_back((*iter).first);
    else if ( (*iter).second == 0 ) fSigNames.push_back((*iter).first);
    else fDataNames.push_back((*iter).first); 
  }
  fNBkg  = fBkgNames.size();
  fNSig  = fSigNames.size();
  fNData = fDataNames.size();

  // initialize histo names 
  ABCDMethod::InitVariables();
  fNTH1D = fTH1DNames.size();
  fNTH2D = fTH2DNames.size();

  ABCDMethod::InitHists();


}

ABCDMethod::~ABCDMethod(){
  std::cout << "Finished ABCD & Deleting" << std::endl;


  // delete input histos
  for (UInt_t th1d = 0; th1d < fNTH1D; th1d++){
    for (UInt_t data = 0; data < fNData; data++) delete fInDataTH1DHists[th1d][data];
    for (UInt_t mc = 0; mc < fNBkg; mc++) delete fInBkgTH1DHists[th1d][mc];
    for (UInt_t mc = 0; mc < fNSig; mc++) delete fInSigTH1DHists[th1d][mc];
  }
  for (UInt_t th2d = 0; th2d < fNTH2D; th2d++){
    for (UInt_t data = 0; data < fNData; data++) delete fInDataTH2DHists[th2d][data];
    for (UInt_t mc = 0; mc < fNBkg; mc++) delete fInBkgTH2DHists[th2d][mc];
    for (UInt_t mc = 0; mc < fNSig; mc++) delete fInSigTH2DHists[th2d][mc];
  }
 
  delete fOutFile;
}


void ABCDMethod::DoAnalysis(){

  // add all data files together
  // data : copy first histogram & add all others too it 
  for (UInt_t data = 0; data < fNData; data++){
    if (data == 0) fOutDataTH2DHists[0] = (TH2D*)fInDataTH2DHists[0][data]->Clone(); 
    else fOutDataTH2DHists[0]->Add(fInDataTH2DHists[0][data]);
  }

  // scale bkg and then make one copy of histos where bkg added together
  for (UInt_t mc = 0; mc < fNBkg; mc++){
    fInBkgTH2DHists[0][mc]->Scale(lumi);
    if (mc == 0) fOutBkgTH2DHists[0] = (TH2D*)fInBkgTH2DHists[0][mc]->Clone();
    else fOutBkgTH2DHists[0]->Add(fInBkgTH2DHists[0][mc]);
  } 

  // just scale the signal by lumi, don't add together 
  for (UInt_t mc = 0; mc < fNSig; mc++){
    fInSigTH2DHists[0][mc]->Scale(lumi);
  }
 
  Int_t fNCat = 6; // for 6 categories A1,B1,A2,B2,D,C
  Data_Int.resize(fNCat);
  Data_IntErr.resize(fNCat);
  Bkg_Int.resize(fNCat);
  Bkg_IntErr.resize(fNCat);
  Sig_Int.resize(fNCat);
  Sig_IntErr.resize(fNCat);

  DblVec min_x;
  DblVec max_x;
  DblVec min_y;
  DblVec max_y;
  min_x.resize(fNCat);
  max_x.resize(fNCat);
  min_y.resize(fNCat);
  max_y.resize(fNCat); 
 
  min_x[0]=mgg_minAB1; // cat0 = A1
  min_x[1]=mgg_minAB1; // cat1 = B1
  min_x[2]=mgg_maxCD;  // cat2 = A2
  min_x[3]=mgg_maxCD;  // cat3 = B2
  min_x[4]=mgg_minCD;  // cat4 = D 
  min_x[5]=mgg_minCD;  // cat5 = C 

  max_x[0]=mgg_minCD;  // cat0 = A1
  max_x[1]=mgg_minCD;  // cat1 = B1
  max_x[2]=mgg_maxAB2; // cat2 = A2
  max_x[3]=mgg_maxAB2; // cat3 = B2
  max_x[4]=mgg_maxCD;  // cat4 = D 
  max_x[5]=mgg_maxCD;  // cat5 = C 

  min_y[0]=met_minD;   // cat0 = A1
  min_y[1]=met_minB;   // cat1 = B1
  min_y[2]=met_minD;   // cat2 = A2
  min_y[3]=met_minB;   // cat3 = B2
  min_y[4]=met_minD;   // cat4 = D 
  min_y[5]=met_minB;   // cat5 = C 

  max_y[0]=met_maxD;   // cat0 = A1
  max_y[1]=met_minD;   // cat1 = B1
  max_y[2]=met_maxD;   // cat2 = A2
  max_y[3]=met_minD;   // cat3 = B2
  max_y[4]=met_maxD;   // cat4 = D 
  max_y[5]=met_minD;   // cat5 = C 

 
  for (UInt_t cat = 0; cat < fNCat; cat++){ // loop over each category
    Data_Int[cat].resize(1); 		// only one group for data since it is lumped together
    Data_IntErr[cat].resize(1);		
    Bkg_Int[cat].resize(fNBkg+1);	// do all Bkg separately and then one where all combined
    Bkg_IntErr[cat].resize(fNBkg+1);
    Sig_Int[cat].resize(fNSig);		// do all Sig separately
    Sig_IntErr[cat].resize(fNSig);

    Data_Int[cat][0] = ABCDMethod::ComputeIntAndErr( fOutBkgTH2DHists[0], Data_IntErr[cat][0],  min_x[cat], max_x[cat], min_y[cat], max_y[cat], cat); 
    for (UInt_t mc = 0; mc < fNBkg; mc++){ 
      Bkg_Int[cat][mc] = ABCDMethod::ComputeIntAndErr( fInBkgTH2DHists[0][mc], Bkg_IntErr[cat][mc],  min_x[cat], max_x[cat], min_y[cat], max_y[cat], cat); 
    }
    // after finished with bkg samples separately, look at the combined sample
    Bkg_Int[cat][fNBkg+1] = ABCDMethod::ComputeIntAndErr( fOutBkgTH2DHists[0], Bkg_IntErr[cat][fNBkg+1],  min_x[cat], max_x[cat], min_y[cat], max_y[cat], cat);
    for (UInt_t mc = 0; mc < fNSig; mc++){ 
      Sig_Int[cat][mc] = ABCDMethod::ComputeIntAndErr( fInSigTH2DHists[0][mc], Sig_IntErr[cat][mc],  min_x[cat], max_x[cat], min_y[cat], max_y[cat], cat); 
    } 

    //std::cout << "Data " << Data_Int[cat][0] << " " << Data_IntErr[cat][0] << std::endl;
    //for (UInt_t mc = 0; mc < fNBkg+1; mc++){ std::cout << "Bkg " << Bkg_Int[cat][mc] << " " << Bkg_IntErr[cat][mc] << std::endl; }
    //for (UInt_t mc = 0; mc < fNSig; mc++){   std::cout << "Sig " << Sig_Int[cat][mc] << " " << Sig_IntErr[cat][mc] << std::endl; }

  }// end cat loop    

    
    ABCDMethod::FillTable( "Data",0, Data_Int[0][0], Data_IntErr[0][0]);

}

void ABCDMethod::FillTable( const TString fSampleName, const UInt_t reg, const UInt_t Integral, const UInt_t Error){


}



Double_t ABCDMethod::ComputeIntAndErr(TH2D *& h, Double_t & error, const Double_t minX, const Double_t maxX, const Double_t minY, const Double_t maxY, const UInt_t isReg ){

  Double_t integral = 0.;

  if(h == (TH2D*) NULL) std::cout << "NULL TH2D" << std::endl;

  //std::cout << "minx = " << minX << " maxX = " << maxX << " minY = " << minY << " maxY = " << maxY << std::endl; 

  Int_t binXmin;
  Int_t binXmax;
  Int_t binYmin;
  Int_t binYmax;    
  if (isReg == 4 || isReg == 5){ // if signal find the exact bins
    binXmin = h->GetXaxis()->FindBin(minX);
    binXmax = h->GetXaxis()->FindBin(maxX);
    binYmin = h->GetYaxis()->FindBin(minY);
    if (isReg == 4) binYmax = h->GetYaxis()->FindBin(maxY);
    if (isReg == 5) binYmax = h->GetYaxis()->FindBin(maxY)-1;
  }
  else if (isReg == 0 || isReg == 1){ // if to left of signal region 
    binXmin = h->GetXaxis()->FindBin(minX);
    binXmax = h->GetXaxis()->FindBin(maxX)-1;
    binYmin = h->GetYaxis()->FindBin(minY);
    if (isReg == 0) binYmax = h->GetYaxis()->FindBin(maxY);
    if (isReg == 1) binYmax = h->GetYaxis()->FindBin(maxY)-1;
  }
  else if (isReg == 2 || isReg == 3){ // if to right of signal region 
    binXmin = h->GetXaxis()->FindBin(minX)+1;
    binXmax = h->GetXaxis()->FindBin(maxX);
    binYmin = h->GetYaxis()->FindBin(minY);
    if (isReg == 2) binYmax = h->GetYaxis()->FindBin(maxY);
    if (isReg == 3) binYmax = h->GetYaxis()->FindBin(maxY)-1;
  }
 
  //std::cout << "binXmin " << binXmin << std::endl;
  //std::cout << "binXmax " << binXmax << std::endl;
  //std::cout << "binYmin " << binYmin << std::endl;
  //std::cout << "binYmax " << binYmax << std::endl;

  integral = h->IntegralAndError(binXmin,binXmax,binYmin,binYmax,error);
  //std::cout << "integral = " << integral << " error = " << error << std::endl;
  return integral;
}

void ABCDMethod::InitHists(){
  // open input files into TFileVec for data
  fDataFiles.resize(fNData);
  for (UInt_t data = 0; data < fNData; data++) {
    TString datafile = Form("%s%s/plots_%s.root",fInDir.Data(),fDataNames[data].Data(),fDataNames[data].Data());
    fDataFiles[data] = TFile::Open(datafile.Data());
    CheckValidFile(fDataFiles[data],datafile);
  }

  // open input files into TFileVec for bkg
  fBkgFiles.resize(fNBkg);
  for (UInt_t mc = 0; mc < fNBkg; mc++) {
    TString bkgfile = Form("%s%s/plots_%s.root",fInDir.Data(),fBkgNames[mc].Data(),fBkgNames[mc].Data());
    fBkgFiles[mc] = TFile::Open(bkgfile.Data());
    CheckValidFile(fBkgFiles[mc],bkgfile);
  }

  // open input files into TFileVec for bkg
  fSigFiles.resize(fNSig);
  for (UInt_t mc = 0; mc < fNSig; mc++) {
    TString sigfile = Form("%s%s/plots_%s.root",fInDir.Data(),fSigNames[mc].Data(),fSigNames[mc].Data());
    fSigFiles[mc] = TFile::Open(sigfile.Data());
    CheckValidFile(fSigFiles[mc],sigfile);
  }

  fInDataTH1DHists.resize(fNTH1D);
  fInBkgTH1DHists.resize(fNTH1D);
  fInSigTH1DHists.resize(fNTH1D);

  for (UInt_t th1d = 0; th1d < fNTH1D; th1d++){ // loop over 1d hists
    fInDataTH1DHists[th1d].resize(fNData); 
    for (UInt_t data = 0; data < fNData; data++) { // init data double hists
      fInDataTH1DHists[th1d][data] = (TH1D*)fDataFiles[data]->Get(Form("%s",fTH1DNames[th1d].Data()));
      CheckValidTH1D(fInDataTH1DHists[th1d][data],fTH1DNames[th1d],fDataFiles[data]->GetName());
    }
    fInBkgTH1DHists[th1d].resize(fNBkg); 
    for (UInt_t mc = 0; mc < fNBkg; mc++) { // init bkg double hists
      fInBkgTH1DHists[th1d][mc] = (TH1D*)fBkgFiles[mc]->Get(Form("%s",fTH1DNames[th1d].Data()));
      CheckValidTH1D(fInBkgTH1DHists[th1d][mc],fTH1DNames[th1d],fBkgFiles[mc]->GetName());
    }
    fInSigTH1DHists[th1d].resize(fNSig); 
    for (UInt_t mc = 0; mc < fNSig; mc++) { // init sig double hists
      fInSigTH1DHists[th1d][mc] = (TH1D*)fSigFiles[mc]->Get(Form("%s",fTH1DNames[th1d].Data()));
      CheckValidTH1D(fInSigTH1DHists[th1d][mc],fTH1DNames[th1d],fSigFiles[mc]->GetName());
    }
  }

  fInDataTH2DHists.resize(fNTH2D);
  fInBkgTH2DHists.resize(fNTH2D);
  fInSigTH2DHists.resize(fNTH2D);

  for (UInt_t th2d = 0; th2d < fNTH2D; th2d++){ // loop over 1d hists
    fInDataTH2DHists[th2d].resize(fNData); 
    for (UInt_t data = 0; data < fNData; data++) { // init data double hists
      fInDataTH2DHists[th2d][data] = (TH2D*)fDataFiles[data]->Get(Form("%s",fTH2DNames[th2d].Data()));
      CheckValidTH2D(fInDataTH2DHists[th2d][data],fTH2DNames[th2d],fDataFiles[data]->GetName());
    }
    fInBkgTH2DHists[th2d].resize(fNBkg); 
    for (UInt_t mc = 0; mc < fNBkg; mc++) { // init bkg double hists
      fInBkgTH2DHists[th2d][mc] = (TH2D*)fBkgFiles[mc]->Get(Form("%s",fTH2DNames[th2d].Data()));
      CheckValidTH2D(fInBkgTH2DHists[th2d][mc],fTH2DNames[th2d],fBkgFiles[mc]->GetName());
    }
    fInSigTH2DHists[th2d].resize(fNSig); 
    for (UInt_t mc = 0; mc < fNSig; mc++) { // init sig double hists
      fInSigTH2DHists[th2d][mc] = (TH2D*)fSigFiles[mc]->Get(Form("%s",fTH2DNames[th2d].Data()));
      CheckValidTH2D(fInSigTH2DHists[th2d][mc],fTH2DNames[th2d],fSigFiles[mc]->GetName());
    }
  }

  fOutDataTH2DHists.resize(fNTH2D);
  fOutBkgTH2DHists.resize(fNTH2D);

}

void ABCDMethod::InitVariables(){
  // 1D histograms of interest
  fTH1DNames.push_back("mgg");
  fTH1DNames.push_back("t1pfmet");

  // 2D histograms of interest
  fTH2DNames.push_back("t1pfmet_mgg");
}

