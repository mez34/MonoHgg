#include "Combiner.hh"

Combiner::Combiner( SamplePairVec Samples, const Double_t inLumi, const ColorMap colorMap, const TString outdir){

  lumi	= inLumi;
  fOutDir = outdir;
  TString fOut = "comb";

  //MakeOutDirectory(Form("%s%s",fOutDir.Data(),fOut.Data()));
  fOutFile = new TFile(Form("%s%s/combplots.root",fOutDir.Data(),fOut.Data()));
 

  for (SamplePairVecIter iter = Samples.begin(); iter != Samples.end(); ++iter){
    if ( (*iter).second == 1 ) {fBkgNames.push_back((*iter).first);} // background
    if ( (*iter).second == 0 ) {fSigNames.push_back((*iter).first);} // signal
    else {fDataNames.push_back((*iter).first);}			     // data
  }
  
  fNData = fDataNames.size();
  fNBkg  = fBkgNames.size();
  fNSig  = fSigNames.size();
  Combiner::InitTH1DNames();
  fNTH1D = 2*fTH1DNames.size(); // 2x for regular and n-1 plots

  // define colorMap and title
  fColorMap = colorMap;
  fSampleTitleMap["qcd"] 	= "QCD";
  fSampleTitleMap["gjets"]	= "G + Jets";
  fSampleTitleMap["gluglu"]	= "GluGlu #rightarrow H #rightarrow GG";
  fSampleTitleMap["dmhgg1"]	= "DM + H #rightarrow GG, M1000GeV";
  fSampleTitleMap["dmhgg2"]	= "DM + H #rightarrow GG, M100GeV";
  fSampleTitleMap["dmhgg3"]	= "DM + H #rightarrow GG, M10GeV";
  fSampleTitleMap["dmhgg4"]	= "DM + H #rightarrow GG, M1GeV";

  // open input files into TFileVec for data
  fDataFiles.resize(fNData);
  for (UInt_t data = 0; data < fNData; data++) {
    TString datafile = Form("%s/%s/plots_%s.root",fOutDir.Data(),fDataNames[data].Data(),fDataNames[data].Data());
    fDataFiles[data] = TFile::Open(datafile.Data());
  }

  // open input files into TFileVec for bkg
  fBkgFiles.resize(fNBkg);
  for (UInt_t mc = 0; mc < fNBkg; mc++) {
    TString bkgfile = Form("%s/%s/plots_%s.root",fOutDir.Data(),fBkgNames[mc].Data(),fBkgNames[mc].Data());
    fBkgFiles[mc] = TFile::Open(bkgfile.Data());
  }

  // open input files into TFileVec for bkg
  fSigFiles.resize(fNSig);
  for (UInt_t mc = 0; mc < fNSig; mc++) {
    TString sigfile = Form("%s/%s/plots_%s.root",fOutDir.Data(),fSigNames[mc].Data(),fSigNames[mc].Data());
    fSigFiles[mc] = TFile::Open(sigfile.Data());
  }

  fInBkgTH1DHists.resize(fNTH1D);
/*  for (UInt_t th1d = 0; th1d < fNTH1D; th1d++){ // loop over double hists
    for (UInt_t mc = 0; mc < fNBkg; mc++) { // init mc double hists
      //fInBkgTH1DHists[th1d][mc] = new TH1D(Form("%s",fTH1DNames[th1d].Data()), Form("%s",fTH1DNames[th1d].Data()),100,0,100);
 
      fInBkgTH1DHists[th1d][mc] = (TH1D*)fBkgFiles[mc]->Get(Form("%s",fTH1DNames[th1d].Data()));
      //fInBkgTH1DHists[th1d][mc]->SetFillColor(fColorMap[fBkgNames[mc]]);
      //fInBkgTH1DHists[th1d][mc]->SetLineColor(kBlack);
    }
  }
*/
/*  fInDataTH1DHists.resize(fNTH1D);
  fInBkgTH1DHists.resize(fNTH1D);
  fInSigTH1DHists.resize(fNTH1D);
  
  for (UInt_t th1d = 0; th1d < fNTH1D; th1d++){ // loop over double hists
    // data first
    fInDataTH1DHists[th1d].resize(fNData); 
    for (UInt_t data = 0; data < fNData; data++) { // init data double hists
      fInDataTH1DHists[th1d][data] = new TH1D(Form("%s",fTH1DNames[th1d].Data()));
      //fInDataTH1DHists[th1d][data] = (TH1D*)fDataFiles[data]->Get(Form("%s",fTH1DNames[th1d].Data()));
    }

    // bkg second
    fInBkgTH1DHists[th1d].resize(fNBkg); 
    for (UInt_t mc = 0; mc < fNBkg; mc++) { // init mc double hists
      fInBkgTH1DHists[th1d][mc] = new (TH1D*)fBkgFiles[mc]->Get(Form("%s",fTH1DNames[th1d].Data()));
      fInBkgTH1DHists[th1d][mc]->SetFillColor(fColorMap[fMCNames[mc]]);
      fInBkgTH1DHists[th1d][mc]->SetLineColor(kBlack);
    }

    // sig second
    fInSigTH1DHists[th1d].resize(fNSig); 
    for (UInt_t mc = 0; mc < fNSig; mc++) { // init mc double hists
      fInSigTH1DHists[th1d][mc] = new (TH1D*)fSigFiles[mc]->Get(Form("%s",fTH1DNames[th1d].Data()));
      fInSigTH1DHists[th1d][mc]->SetFillColor(fColorMap[fMCNames[mc+3]]); //FIXME sig are after first 3 bkg samples
      fInSigTH1DHists[th1d][mc]->SetLineColor(kBlack);
    }
  }*/

/*  fOutDataTH1DHists.resize(fNTH1D);
  fOutBkgTH1DHists.resize(fNTH1D);
  fOutSigTH1DHists.resize(fNTH1D);
  for (UInt_t th1d = 0; th1d < fNTH1D; th1d++){
    fOutBkgTH1DStacks[th1d] = new THStack("","");
  }*/

}// end Combiner::Combiner

Combiner::~Combiner(){
  std::cout << "Finished & Deleting" << std::endl;

 // delete all pointers
  for (UInt_t th1d = 0; th1d < fNTH1D; th1d++){
    delete fOutDataTH1DHists[th1d];
    delete fOutBkgTH1DHists[th1d];
    delete fOutSigTH1DHists[th1d];
//    delete fOutBkgTH1DStacks[th1d];
    delete fTH1DLegends[th1d];
    delete fOutTH1DCanvases[th1d];
    
    /*for (UInt_t data = 0; data < fNData; data++) {
      delete fInDataTH1DHists[th1d][data];
    }
   */ 
   delete fInBkgTH1DHists[th1d];
    for (UInt_t mc = 0; mc < fNBkg; mc++) {

//      delete fInBkgTH1DHists[th1d][mc];
    }
   /*
    for (UInt_t mc = 0; mc < fNSig; mc++) {
      delete fInSigTH1DHists[th1d][mc];
    }*/
  }

  for (UInt_t data = 0; data < fNData; data++) {
    delete fDataFiles[data];
  }

  for (UInt_t mc = 0; mc < fNBkg; mc++) {
    delete fBkgFiles[mc];
  }

  for (UInt_t mc = 0; mc < fNSig; mc++) {
    delete fSigFiles[mc];
  }

  delete fOutFile;

}// end Combiner::~Combiner

void Combiner::DoComb(){
  Combiner::OverlayPlots();
  Combiner::StackPlots();
}// end Combiner::DoComb


void Combiner::OverlayPlots(){

}// end Combiner::OverlayPlots

void Combiner::StackPlots(){

}// end Combiner::StackPlots

void Combiner::InitTH1DNames(){
  // higgs & met variables
  fTH1DNames.push_back("mgg");
  fTH1DNames.push_back("ptgg");
  fTH1DNames.push_back("phiH");
  fTH1DNames.push_back("t1pfmetPhi");
  fTH1DNames.push_back("phiHMET");
  fTH1DNames.push_back("t1pfmet");
  fTH1DNames.push_back("nvtx");

  // photon variables
  fTH1DNames.push_back("pt1");
  fTH1DNames.push_back("pt2");
  fTH1DNames.push_back("eta1");
  fTH1DNames.push_back("eta2");
  fTH1DNames.push_back("phi1");
  fTH1DNames.push_back("phi2");

  fTH1DNames.push_back("r91");
  fTH1DNames.push_back("r92");

  // photon ID variables
  fTH1DNames.push_back("hoe1");
  fTH1DNames.push_back("hoe2");
  fTH1DNames.push_back("sieie1");
  fTH1DNames.push_back("sieie2");
  fTH1DNames.push_back("phoiso1");
  fTH1DNames.push_back("phoiso2");
  fTH1DNames.push_back("chiso1");
  fTH1DNames.push_back("chiso2");
  fTH1DNames.push_back("neuiso1");
  fTH1DNames.push_back("neuiso2");

}// end Combiner::InitTH1DNames
