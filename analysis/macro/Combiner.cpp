#include "Combiner.hh"

Combiner::Combiner( SamplePairVec Samples, const Double_t inLumi, const ColorMap colorMap, const TString outdir, const Bool_t doNmin1){

  if (doNmin1) addText = "_n-1";
  else addText="";

  lumi	= inLumi;
  fOutDir = outdir;
  TString fOut = "comb";

  MakeOutDirectory(Form("%s%s",fOutDir.Data(),fOut.Data()));
  fOutFile = new TFile(Form("%s%s/combplots%s.root",fOutDir.Data(),fOut.Data(),addText.Data()),"RECREATE");
  CheckValidFile(fOutFile, Form("%s%s/combplots%s.root",fOutDir.Data(),fOut.Data(),addText.Data())); 

  for (SamplePairVecIter iter = Samples.begin(); iter != Samples.end(); ++iter){
    if ( (*iter).second == 1 ) {fBkgNames.push_back((*iter).first);}      // background
    else if ( (*iter).second == 0 ) {fSigNames.push_back((*iter).first);} // signal
    else {fDataNames.push_back((*iter).first);}			          // data
  }
  
  fNData = fDataNames.size();
  fNBkg  = fBkgNames.size();
  fNSig  = fSigNames.size();
  Combiner::InitTH1DNames();
  fNTH1D = fTH1DNames.size(); // 2x for regular and n-1 plots

  // define colorMap and title
  fColorMap = colorMap;
  
  fSampleTitleMap["QCD"] 		= "QCD";
  fSampleTitleMap["GJets"]		= "G + Jets";
  fSampleTitleMap["GluGluHToGG"]	= "GluGlu #rightarrow H #rightarrow GG";
  fSampleTitleMap["DMHtoGG_M1000"]	= "DM + H #rightarrow GG, M1000GeV";
  fSampleTitleMap["DMHtoGG_M100"]	= "DM + H #rightarrow GG, M100GeV";
  fSampleTitleMap["DMHtoGG_M10"]	= "DM + H #rightarrow GG, M10GeV";
  fSampleTitleMap["DMHtoGG_M1"]		= "DM + H #rightarrow GG, M1GeV";

  //for (std::map<TString,TString>::iterator iter = fSampleTitleMap.begin(); iter != fSampleTitleMap.end(); ++iter) {
  //  std::cout << (*iter).first << "  " << (*iter).second << std::endl;
  //}

  // open input files into TFileVec for data
  fDataFiles.resize(fNData);
  for (UInt_t data = 0; data < fNData; data++) {
    TString datafile = Form("%s%s/plots_%s.root",fOutDir.Data(),fDataNames[data].Data(),fDataNames[data].Data());
    fDataFiles[data] = TFile::Open(datafile.Data());
    CheckValidFile(fDataFiles[data],datafile);
  }

  // open input files into TFileVec for bkg
  fBkgFiles.resize(fNBkg);
  for (UInt_t mc = 0; mc < fNBkg; mc++) {
    TString bkgfile = Form("%s%s/plots_%s.root",fOutDir.Data(),fBkgNames[mc].Data(),fBkgNames[mc].Data());
    fBkgFiles[mc] = TFile::Open(bkgfile.Data());
    CheckValidFile(fBkgFiles[mc],bkgfile);
  }

  // open input files into TFileVec for bkg
  fSigFiles.resize(fNSig);
  for (UInt_t mc = 0; mc < fNSig; mc++) {
    TString sigfile = Form("%s%s/plots_%s.root",fOutDir.Data(),fSigNames[mc].Data(),fSigNames[mc].Data());
    fSigFiles[mc] = TFile::Open(sigfile.Data());
    CheckValidFile(fSigFiles[mc],sigfile);
  }

  fInDataTH1DHists.resize(fNTH1D);
  fInBkgTH1DHists.resize(fNTH1D);
  fInSigTH1DHists.resize(fNTH1D);


  for (UInt_t th1d = 0; th1d < fNTH1D; th1d++){ // loop over double hists
    fInDataTH1DHists[th1d].resize(fNData); 
    for (UInt_t data = 0; data < fNData; data++) { // init data double hists
      fInDataTH1DHists[th1d][data] = (TH1D*)fDataFiles[data]->Get(Form("%s_%s%s",fTH1DNames[th1d].Data(),fDataNames[data].Data(),addText.Data()));
      CheckValidTH1D(fInDataTH1DHists[th1d][data],fTH1DNames[th1d],fDataFiles[data]->GetName());
    }
    fInBkgTH1DHists[th1d].resize(fNBkg); 
    for (UInt_t mc = 0; mc < fNBkg; mc++) { // init bkg double hists
      fInBkgTH1DHists[th1d][mc] = (TH1D*)fBkgFiles[mc]->Get(Form("%s_%s%s",fTH1DNames[th1d].Data(),fBkgNames[mc].Data(),addText.Data()));
      CheckValidTH1D(fInBkgTH1DHists[th1d][mc],fTH1DNames[th1d],fBkgFiles[mc]->GetName());
      fInBkgTH1DHists[th1d][mc]->SetFillColor(fColorMap[fBkgNames[mc]]);
      fInBkgTH1DHists[th1d][mc]->SetLineColor(kBlack);
    }
    fInSigTH1DHists[th1d].resize(fNSig); 
    for (UInt_t mc = 0; mc < fNSig; mc++) { // init sig double hists
      fInSigTH1DHists[th1d][mc] = (TH1D*)fSigFiles[mc]->Get(Form("%s_%s%s",fTH1DNames[th1d].Data(),fSigNames[mc].Data(),addText.Data()));
      CheckValidTH1D(fInSigTH1DHists[th1d][mc],fTH1DNames[th1d],fSigFiles[mc]->GetName());
      fInSigTH1DHists[th1d][mc]->SetLineColor(fColorMap[fSigNames[mc]]);
    }
  }

  // output histos
  if (fNData > 0) fOutDataTH1DHists.resize(fNTH1D);
  fOutBkgTH1DStacks.resize(fNTH1D);
  for (UInt_t th1d = 0; th1d < fNTH1D; th1d++){
    fOutBkgTH1DStacks[th1d] = new THStack("","");
  }

  fTH1DLegends.resize(fNTH1D);
  for (UInt_t th1d = 0; th1d < fNTH1D; th1d++){
    fTH1DLegends[th1d] = new TLegend(0.70,0.7,0.9,0.89); // (x1,y1,x2,y2)
    fTH1DLegends[th1d]->SetBorderSize(4);
    fTH1DLegends[th1d]->SetLineColor(kBlack);
  }

  fOutTH1DCanvases.resize(fNTH1D);
  for (UInt_t th1d = 0; th1d < fNTH1D; th1d++){
    fOutTH1DCanvases[th1d] = new TCanvas(fTH1DNames[th1d].Data(),"");
    fOutTH1DCanvases[th1d]->cd();
  }
  fOutTH1DStackPads.resize(fNTH1D);
  for (UInt_t th1d = 0; th1d < fNTH1D; th1d++){
    fOutTH1DStackPads[th1d] = new TPad("","",0,0.3,1.0,0.99);
    fOutTH1DStackPads[th1d]->SetBottomMargin(0); // upper and lower pad are joined
  }

  // for ratio plots
  fOutBkgTH1DHists.resize(fNTH1D);
  if (fNData > 0) fOutRatioTH1DHists.resize(fNTH1D); 

}// end Combiner::Combiner

Combiner::~Combiner(){
  std::cout << "Finished & Deleting" << std::endl;

 // delete all pointers
  for (UInt_t th1d = 0; th1d < fNTH1D; th1d++){
    if (fNData > 0 ) delete fOutDataTH1DHists[th1d];
    delete fOutBkgTH1DHists[th1d];
    delete fOutBkgTH1DStacks[th1d];
    delete fTH1DLegends[th1d];
    delete fOutTH1DStackPads[th1d];
    delete fOutTH1DCanvases[th1d];
    
    for (UInt_t data = 0; data < fNData; data++) { 
      delete fInDataTH1DHists[th1d][data];
    }
    for (UInt_t mc = 0; mc < fNBkg; mc++) { 
      delete fInBkgTH1DHists[th1d][mc];
    }
    for (UInt_t mc = 0; mc < fNSig; mc++) {
      delete fInSigTH1DHists[th1d][mc];
    }
  }

  for (UInt_t data = 0; data < fNData; data++) { delete fDataFiles[data]; }
  for (UInt_t mc = 0; mc < fNBkg; mc++) { delete fBkgFiles[mc]; }
  for (UInt_t mc = 0; mc < fNSig; mc++) { delete fSigFiles[mc]; }

  delete fOutFile;

}// end Combiner::~Combiner

void Combiner::DoComb(){

  Combiner::OverlayPlots();
}// end Combiner::DoComb


void Combiner::OverlayPlots(){
  // copy th1d plots into output hists/stacks
  for (UInt_t th1d = 0; th1d < fNTH1D; th1d++){

    // data : copy first histogram & add all others too it 
    if (fNData > 0){
      for (UInt_t data = 0; data < fNData; data++){
        if (data == 0) fOutDataTH1DHists[th1d] = (TH1D*)fInDataTH1DHists[th1d][data]->Clone();
        else fOutDataTH1DHists[th1d]->Add(fInDataTH1DHists[th1d][data]);
      }
      //fTH1DLegends[th1d]->AddEntry(fOutDataTH1DHists[th1d],"Data","pl"); //add data entry to legend
    }// end if ndata>0

    // bkg : copy histos and add to stacks
    for (UInt_t mc = 0; mc < fNBkg; mc++){
      fInBkgTH1DHists[th1d][mc]->Scale(lumi);
      fOutBkgTH1DStacks[th1d]->Add(fInBkgTH1DHists[th1d][mc]);
      fTH1DLegends[th1d]->AddEntry(fInBkgTH1DHists[th1d][mc],fSampleTitleMap[fBkgNames[mc]],"lf");
      if (mc == 0){
        fOutBkgTH1DHists[th1d] = (TH1D*)fInBkgTH1DHists[th1d][mc]->Clone();
      }
      else{
        fOutBkgTH1DHists[th1d]->Add(fInBkgTH1DHists[th1d][mc]);
      }
    } 
  
    // sig: just add to legend
    for (UInt_t mc = 0; mc < fNSig; mc++){
      fTH1DLegends[th1d]->AddEntry(fInSigTH1DHists[th1d][mc],fSampleTitleMap[fSigNames[mc]],"l");
    }
  }// end loop over th1d histos
  Combiner::MakeOutputCanvas();

}// end Combiner::OverlayPlots


void Combiner::MakeOutputCanvas(){
  for (UInt_t th1d = 0; th1d < fNTH1D; th1d++){
    Bool_t isLogY = true;
    // do stack plots first
    Combiner::DrawCanvasStack(th1d,isLogY);
    isLogY = false;
    Combiner::DrawCanvasStack(th1d,isLogY);
    isLogY = true;
    // do overlay next 
    isLogY = true;
    Combiner::DrawCanvasOverlay(th1d,isLogY);
    isLogY = false;
    Combiner::DrawCanvasOverlay(th1d,isLogY);
  }
}// end Combiner::MakeOutputCanvas

void Combiner::DrawCanvasOverlay(const UInt_t th1d, const Bool_t isLogY){
  gStyle->SetOptStat(0);
  fOutTH1DCanvases[th1d]->cd();
  fOutTH1DStackPads[th1d]->Draw();
  fOutTH1DStackPads[th1d]->cd();

   
  for (UInt_t mc = 0; mc < fNSig; mc++){
    if (fInSigTH1DHists[th1d][mc]->Integral() > 0){
      fInSigTH1DHists[th1d][mc]->Scale(1.0/fInSigTH1DHists[th1d][mc]->Integral());
    }
  }
  for (UInt_t mc = 0; mc < fNBkg; mc++){
    if (fInBkgTH1DHists[th1d][mc]->Integral() > 0 ){
      fInBkgTH1DHists[th1d][mc]->Scale(1.0/fInBkgTH1DHists[th1d][mc]->Integral());
    }
    fInBkgTH1DHists[th1d][mc]->SetFillColor(0);
    fInBkgTH1DHists[th1d][mc]->SetLineColor(fColorMap[fBkgNames[mc]]);
  }

  Double_t max = -100;
  max = Combiner::GetMaximum(th1d, false);

  Double_t minval = 1E9;
  minval = Combiner::GetMinimum(th1d, false);

  // start by drawing the sig first
  if (isLogY) fInSigTH1DHists[th1d][0]->SetMaximum(max*10);
  else fInSigTH1DHists[th1d][0]->SetMaximum(max*1.1);
  fInSigTH1DHists[th1d][0]->SetTitle("");
  fInSigTH1DHists[th1d][0]->Draw("hist");

  for (UInt_t mc = 0; mc < fNBkg; mc++){
    fInBkgTH1DHists[th1d][mc]->Draw("HIST SAME");
  }
  for (UInt_t mc = 0; mc < fNSig; mc++){
    fInSigTH1DHists[th1d][mc]->Draw("HIST SAME");
  }
  fTH1DLegends[th1d]->Draw("SAME"); 

  TString suffix = "";
  if (isLogY) suffix="_log";

  fOutTH1DStackPads[th1d]->SetLogy(isLogY);
  fOutTH1DCanvases[th1d]->cd();

  CMSLumi(fOutTH1DCanvases[th1d],0,lumi);

  fOutTH1DCanvases[th1d]->SaveAs(Form("%scomb/%s_comb%s%s.png",fOutDir.Data(),fTH1DNames[th1d].Data(),addText.Data(),suffix.Data()));  
  fOutFile->cd();
  fOutTH1DCanvases[th1d]->Write(Form("%s%s_comb%s",fTH1DNames[th1d].Data(),suffix.Data(),addText.Data()));


}// end Combiner::DrawCanvasOverlay

void Combiner::DrawCanvasStack(const UInt_t th1d, const Bool_t isLogY){
  gStyle->SetOptStat(0);
  fOutTH1DCanvases[th1d]->cd();
  fOutTH1DStackPads[th1d]->Draw();
  fOutTH1DStackPads[th1d]->cd();


  for (UInt_t mc = 0; mc < fNSig; mc++){
    fInSigTH1DHists[th1d][mc]->Scale(lumi);
  }


  Double_t max = -100;
  max = Combiner::GetMaximum(th1d, true);

  Double_t minval = 1E20;
  minval = Combiner::GetMinimum(th1d, true);

  // start by drawing the sig first
  if (isLogY){
    fInSigTH1DHists[th1d][0]->SetMaximum(max*10);
  //  fInSigTH1DHists[th1d][0]->SetMinimum(minval*0.9);
  }
  else{
    fInSigTH1DHists[th1d][0]->SetMaximum(max*1.1);
  //  fInSigTH1DHists[th1d][0]->SetMinimum(minval*0.1);
  }
  fInSigTH1DHists[th1d][0]->SetTitle("");
  fInSigTH1DHists[th1d][0]->Draw("HIST");

  if (fNData > 0) fOutDataTH1DHists[th1d]->Draw("PE SAME");
  fOutBkgTH1DStacks[th1d]->Draw("HIST SAME");

  for (UInt_t mc = 0; mc < fNSig; mc++){
    fInSigTH1DHists[th1d][mc]->Draw("HIST SAME");
    if ( mc == fNSig-1 ){ // on last draw redraw axis because they are overwritten by stack
      fInSigTH1DHists[th1d][mc]->Draw("AXIS SAME");
    }
  }

  fTH1DLegends[th1d]->Draw("SAME"); 

  TString suffix = "";
  if (isLogY) suffix="_log";

  fOutTH1DStackPads[th1d]->SetLogy(isLogY);

  if (fNData > 0){
    Combiner::MakeRatioPlots();
  } 


  fOutTH1DCanvases[th1d]->cd();

  CMSLumi(fOutTH1DCanvases[th1d],0,lumi);

  fOutTH1DCanvases[th1d]->SaveAs(Form("%scomb/%s_stack%s%s.png",fOutDir.Data(),fTH1DNames[th1d].Data(),addText.Data(),suffix.Data()));  
  fOutFile->cd();
  fOutTH1DCanvases[th1d]->Write(Form("%s%s_stack%s",fTH1DNames[th1d].Data(),suffix.Data(),addText.Data()));

}// end Combiner::DrawCanvasStack


void Combiner::MakeRatioPlots(){

  for (UInt_t th1d = 0; th1d < fNTH1D; th1d++){ // double hists
    fOutRatioTH1DHists[th1d] = (TH1D*)fOutDataTH1DHists[th1d]->Clone();
    fOutRatioTH1DHists[th1d]->Divide(fOutBkgTH1DHists[th1d]);  
    fOutRatioTH1DHists[th1d]->SetLineColor(kBlack);
    fOutRatioTH1DHists[th1d]->SetMinimum(-0.1);  // Define Y ..
    fOutRatioTH1DHists[th1d]->SetMaximum(2.1);   // .. range
    fOutRatioTH1DHists[th1d]->SetStats(0);      // No statistics on lower plot
    fOutRatioTH1DHists[th1d]->GetYaxis()->SetTitle("Data/MC");
  }

}// end Combiner::MakeRatioPlots

Double_t Combiner::GetMaximum(const UInt_t th1d, const Bool_t stack) {
  Double_t max = -100;

  std::vector<Double_t> tmpmax;
  if (fNData > 0) tmpmax.push_back(fOutDataTH1DHists[th1d]->GetBinContent(fOutDataTH1DHists[th1d]->GetMaximumBin()));
  for (UInt_t mc = 0; mc < fNSig; mc++){
    tmpmax.push_back( fInSigTH1DHists[th1d][mc]->GetBinContent(fInSigTH1DHists[th1d][mc]->GetMaximumBin()));
  }
  if (stack) tmpmax.push_back(fOutBkgTH1DStacks[th1d]->GetMaximum());
  else{
    for (UInt_t mc = 0; mc < fNBkg; mc++){
      tmpmax.push_back(fInBkgTH1DHists[th1d][mc]->GetBinContent(fInBkgTH1DHists[th1d][mc]->GetMaximumBin()));
    }
  }

  //else{
  //  for (UInt_t mc = 0; mc < fNBkg; mc++){
  //    tmpmax.push_back( fInBkgTH1DHists[th1d][mc]->GetBinContent(fInBkgTH1DHists[th1d][mc]->GetMaximumBin()));
  //  }
  //}

  for (UInt_t i = 0; i < tmpmax.size(); i++){
    if ( tmpmax[i] > max ) max = tmpmax[i];
  }

  return max;
}


Double_t Combiner::GetMinimum(const UInt_t th1d, const Bool_t stack) {
  // need to loop through to check bin != 0
  Double_t datamin  = 1e9;
 
  
  std::vector<Double_t> tmpmin;
//  tmpmax.push_back(fOutDataTH1DHists[th1d]->GetBinContent(fOutDataTH1DHists[th1d]->GetMaximumBin()));
  for (UInt_t mc = 0; mc < fNSig; mc++){
    tmpmin.push_back( fInSigTH1DHists[th1d][mc]->GetBinContent(fInSigTH1DHists[th1d][mc]->GetMinimumBin()));
  }
  if (stack) tmpmin.push_back(fOutBkgTH1DStacks[th1d]->GetMinimum());
  else{
    for (UInt_t mc = 0; mc < fNBkg; mc++){
      tmpmin.push_back(fInBkgTH1DHists[th1d][mc]->GetBinContent(fInBkgTH1DHists[th1d][mc]->GetMinimumBin()));
    }
  }
 
  for (UInt_t i = 0; i < tmpmin.size(); i++){
    if ( tmpmin[i] < datamin ) datamin = tmpmin[i];
  }
//  for (Int_t bin = 1; bin <= fOutDataTH1DHists[th1d]->GetNbinsX(); bin++){
//    Float_t tmpmin = fOutDataTH1DHists[th1d]->GetBinContent(bin);
//    if ((tmpmin < datamin) && (tmpmin > 0)) {
//      datamin    = tmpmin;
//      newdatamin = true;
//    }
//  }
  return datamin;

}

void Combiner::InitTH1DNames(){
  // higgs & met variables
  fTH1DNames.push_back("mgg");
/*  fTH1DNames.push_back("ptgg");
  fTH1DNames.push_back("t1pfmetPhi");
  fTH1DNames.push_back("t1pfmet");
  fTH1DNames.push_back("nvtx");
  if (addText!="_n-1"){ fTH1DNames.push_back("phi_H"); }
  if (addText!="_n-1"){ fTH1DNames.push_back("phi_HMET"); }

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
*/
}// end Combiner::InitTH1DNames
