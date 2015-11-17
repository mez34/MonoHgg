#include "ABCDMethod.hh"

ABCDMethod::ABCDMethod( SamplePairVec Samples, const Double_t inLumi, const TString outdir, Bool_t Blind){

  // load RooFit
  gSystem->Load("libRooFit");
  //using namespace RooFit;

  lumi = inLumi;
  doBlind = Blind;
  fInDir = outdir;
  fOutDir = outdir+"ABCD";

  // values for the different bin edges
  mgg_minAB1 = 100.;
  mgg_minCD  = 120.;
  mgg_maxCD  = 130.;
  mgg_maxAB2 = 150.; 
  met_minB   = 0.;
  met_minD   = 250.;
  met_maxD   = 999.;

  // titles for output Latex table
  fSampleTitleMap["Data"]		= "Data";
  fSampleTitleMap["QCD"] 		= "QCD";
  fSampleTitleMap["GJets"]		= "$\\gamma$ + Jets";
  fSampleTitleMap["VH"]			= "V + H";
  fSampleTitleMap["DYJetsToLL"]		= "Drell-Yan";
  fSampleTitleMap["GluGluHToGG"]	= "$H \\rightarrow \\gamma \\gamma$ (ggH)";
  fSampleTitleMap["DiPhoton"]		= "$\\gamma\\gamma$";
  //fSampleTitleMap["DMHtoGG_M1"]	= "$\\bar{\\chi}\\chi HH, m_{\\chi}$ = 1 GeV (10 fb)";
  //fSampleTitleMap["DMHtoGG_M10"]	= "$\\bar{\\chi}\\chi HH, m_{\\chi}$ = 10 GeV (10 fb)";
  //fSampleTitleMap["DMHtoGG_M100"]	= "$\\bar{\\chi}\\chi HH, m_{\\chi}$ = 100 GeV (10 fb)";
  //fSampleTitleMap["DMHtoGG_M1000"]	= "$\\bar{\\chi}\\chi HH, m_{\\chi}$ = 1000 GeV (10 fb)";
  fSampleTitleMap["2HDM_mZP600"]	= "2HDM, $m_{Z'} = 600 GeV, m_{A0} = 300 GeV$";
  fSampleTitleMap["2HDM_mZP800"]	= "2HDM, $m_{Z'} = 800 GeV, m_{A0} = 300 GeV$";
  fSampleTitleMap["2HDM_mZP1000"]	= "2HDM, $m_{Z'} = 1000 GeV, m_{A0} = 300 GeV$";
  fSampleTitleMap["2HDM_mZP1200"]	= "2HDM, $m_{Z'} = 1200 GeV, m_{A0} = 300 GeV$";
  fSampleTitleMap["2HDM_mZP1400"]	= "2HDM, $m_{Z'} = 1400 GeV, m_{A0} = 300 GeV$";
  fSampleTitleMap["2HDM_mZP1700"]	= "2HDM, $m_{Z'} = 1700 GeV, m_{A0} = 300 GeV$";
  fSampleTitleMap["2HDM_mZP2000"]	= "2HDM, $m_{Z'} = 2000 GeV, m_{A0} = 300 GeV$";
  fSampleTitleMap["2HDM_mZP2500"]	= "2HDM, $m_{Z'} = 2500 GeV, m_{A0} = 300 GeV$";

  // make output txt file with output table
  fOutTableTxtFile.open(Form("%s/ResultsTableForLatex.tex",fOutDir.Data()));

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
    delete fOutDataTH2DHists[th2d];
    delete fOutBkgTH2DHists[th2d];
    delete fOutSelBkgTH2DHists[th2d];
  }

  for (UInt_t cat = 0; cat < 7; cat++){
    delete fRooData[cat][0];
    for (UInt_t mc = 0; mc < fNBkg+2; mc++){ delete fRooBkg[cat][mc]; }
    for (UInt_t mc = 0; mc < fNSig; mc++){ delete fRooSig[cat][mc]; }
  }

   
  delete fOutFile;
  // close output text files
  fOutTableTxtFile.close();
  
  std::cout << "Finished Deleting" << std::endl;
}


void ABCDMethod::DoAnalysis(){

  // add all data files together
  // data : copy first histogram & add all others too it 
  for (UInt_t data = 0; data < fNData; data++){
    //std::cout << "number entries in data in " << fInDataTH2DHists[0][data]->GetEntries() << std::endl;
    if (data == 0) fOutDataTH2DHists[0] = (TH2D*)fInDataTH2DHists[0][data]->Clone(); 
    else fOutDataTH2DHists[0]->Add(fInDataTH2DHists[0][data]);
  }
  //std::cout << "number entries in data" << fOutDataTH2DHists[0]->GetEntries() << std::endl;

  // scale bkg and then make one copy of histos where bkg added together

  fSampleTitleMap["QCD"] 		= "QCD";
  fSampleTitleMap["GJets"]		= "$\\gamma$ + Jets";
  fSampleTitleMap["VH"]			= "V + H";
  fSampleTitleMap["DYJetsToLL"]		= "Drell-Yan";
  fSampleTitleMap["GluGluHToGG"]	= "$H \\rightarrow \\gamma \\gamma$ (ggH)";

  // find indices for the different bkg samples
  for (UInt_t mc = 0; mc < fNBkg; mc++){
    if (fBkgNames[mc] == "VH")		i_vh  = mc;
    if (fBkgNames[mc] == "QCD")		i_qcd = mc;
    if (fBkgNames[mc] == "DiPhoton")	i_gg  = mc;
    if (fBkgNames[mc] == "GJets")	i_gj  = mc;
    if (fBkgNames[mc] == "DYJetsToLL")  i_dy  = mc;
    if (fBkgNames[mc] == "GluGluHToGG") i_hgg = mc;
  }

  // sum over nonresonant bkgs only
  fOutSelBkgTH2DHists[0] = (TH2D*)fInBkgTH2DHists[0][i_dy]->Clone();
  fOutSelBkgTH2DHists[0]->Add(fInBkgTH2DHists[0][i_gg]);
  fOutSelBkgTH2DHists[0]->Add(fInBkgTH2DHists[0][i_gj]); 
  fOutSelBkgTH2DHists[0]->Add(fInBkgTH2DHists[0][i_qcd]);    
 
  for (UInt_t mc = 0; mc < fNBkg; mc++){
    //fInBkgTH2DHists[0][mc]->Scale(3000./1260.);// in order to scale to 3fb-1
    //std::cout << "number entries in bkg in " << fInBkgTH2DHists[0][mc]->GetEntries() << std::endl;
    // use below if summing over all backgrounds
    if (mc == 0) fOutBkgTH2DHists[0] = (TH2D*)fInBkgTH2DHists[0][mc]->Clone();
    else fOutBkgTH2DHists[0]->Add(fInBkgTH2DHists[0][mc]);
  } 
  // just scale the signal by lumi, don't add together 
  for (UInt_t mc = 0; mc < fNSig; mc++){
    //fInSigTH2DHists[0][mc]->Scale(3000./1260.);// in order to scale to 3fb-1
    //std::cout << "number entries in sig in " << fInSigTH2DHists[0][mc]->GetEntries() << std::endl;
  }
 
  Int_t fNCat = 7; // for 6 categories A1,B1,A2,B2,D,C,All
  Data_Int.resize(fNCat);
  Data_IntErr.resize(fNCat);
  Bkg_Int.resize(fNCat);
  Bkg_IntErr.resize(fNCat);
  Sig_Int.resize(fNCat);
  Sig_IntErr.resize(fNCat);

  std::vector<UInt_t> min_x;
  std::vector<UInt_t> max_x;
  std::vector<UInt_t> min_y;
  std::vector<UInt_t> max_y;
  min_x.resize(fNCat);
  max_x.resize(fNCat);
  min_y.resize(fNCat);
  max_y.resize(fNCat); 
 
  // find bin boundaries for each category using BkgTH2DHist[0][0]
  // cat0 = A1 
  min_x[0]=fInBkgTH2DHists[0][0]->GetXaxis()->FindBin(mgg_minAB1); 
  max_x[0]=fInBkgTH2DHists[0][0]->GetXaxis()->FindBin(mgg_minCD)-1;
  min_y[0]=fInBkgTH2DHists[0][0]->GetYaxis()->FindBin(met_minD);   
  max_y[0]=fInBkgTH2DHists[0][0]->GetYaxis()->FindBin(met_maxD);   
  // cat1 = B1
  min_x[1]=fInBkgTH2DHists[0][0]->GetXaxis()->FindBin(mgg_minAB1); 
  max_x[1]=fInBkgTH2DHists[0][0]->GetXaxis()->FindBin(mgg_minCD)-1;
  min_y[1]=fInBkgTH2DHists[0][0]->GetYaxis()->FindBin(met_minB);   
  max_y[1]=fInBkgTH2DHists[0][0]->GetYaxis()->FindBin(met_minD)-1; 
  // cat2 = A2
  min_x[2]=fInBkgTH2DHists[0][0]->GetXaxis()->FindBin(mgg_maxCD)+1;
  max_x[2]=fInBkgTH2DHists[0][0]->GetXaxis()->FindBin(mgg_maxAB2); 
  min_y[2]=fInBkgTH2DHists[0][0]->GetYaxis()->FindBin(met_minD);   
  max_y[2]=fInBkgTH2DHists[0][0]->GetYaxis()->FindBin(met_maxD);   
  // cat3 = B2
  min_x[3]=fInBkgTH2DHists[0][0]->GetXaxis()->FindBin(mgg_maxCD)+1;
  max_x[3]=fInBkgTH2DHists[0][0]->GetXaxis()->FindBin(mgg_maxAB2); 
  min_y[3]=fInBkgTH2DHists[0][0]->GetYaxis()->FindBin(met_minB);   
  max_y[3]=fInBkgTH2DHists[0][0]->GetYaxis()->FindBin(met_minD)-1; 
  // cat4 = D
  min_x[4]=fInBkgTH2DHists[0][0]->GetXaxis()->FindBin(mgg_minCD);  
  max_x[4]=fInBkgTH2DHists[0][0]->GetXaxis()->FindBin(mgg_maxCD);  
  min_y[4]=fInBkgTH2DHists[0][0]->GetYaxis()->FindBin(met_minD);   
  max_y[4]=fInBkgTH2DHists[0][0]->GetYaxis()->FindBin(met_maxD);   
  // cat5 = C
  min_x[5]=fInBkgTH2DHists[0][0]->GetXaxis()->FindBin(mgg_minCD);  
  max_x[5]=fInBkgTH2DHists[0][0]->GetXaxis()->FindBin(mgg_maxCD);  
  min_y[5]=fInBkgTH2DHists[0][0]->GetYaxis()->FindBin(met_minB);   
  max_y[5]=fInBkgTH2DHists[0][0]->GetYaxis()->FindBin(met_minD);   
  // cat6 = ALL
  min_x[6]=fInBkgTH2DHists[0][0]->GetXaxis()->FindBin(mgg_minAB1); 
  max_x[6]=fInBkgTH2DHists[0][0]->GetXaxis()->FindBin(mgg_maxAB2); 
  min_y[6]=fInBkgTH2DHists[0][0]->GetYaxis()->FindBin(met_minB);   
  max_y[6]=fInBkgTH2DHists[0][0]->GetYaxis()->FindBin(met_maxD);   
 
  for (UInt_t cat = 0; cat < fNCat; cat++){
  //  std::cout << "bin Xmin = " << min_x[cat] << std::endl;
  //  std::cout << "bin Xmax = " << max_x[cat] << std::endl;
  //  std::cout << "bin Ymin = " << min_y[cat] << std::endl;
  //  std::cout << "bin Ymax = " << max_y[cat] << std::endl;
  }
 
  for (UInt_t cat = 0; cat < fNCat; cat++){ // loop over each category
    Data_Int[cat].resize(1); 		// only one group for data since it is lumped together
    Data_IntErr[cat].resize(1);		
    Bkg_Int[cat].resize(fNBkg+2);	// do all Bkg separately and then one where all combined, then non-res bkg
    Bkg_IntErr[cat].resize(fNBkg+2);
    Sig_Int[cat].resize(fNSig);		// do all Sig separately
    Sig_IntErr[cat].resize(fNSig);

    Data_Int[cat][0] = ABCDMethod::ComputeIntAndErr( fOutDataTH2DHists[0], Data_IntErr[cat][0],  min_x[cat], max_x[cat], min_y[cat], max_y[cat]); 

    for (UInt_t mc = 0; mc < fNBkg; mc++){ 
      Bkg_Int[cat][mc] = ABCDMethod::ComputeIntAndErr( fInBkgTH2DHists[0][mc], Bkg_IntErr[cat][mc],  min_x[cat], max_x[cat], min_y[cat], max_y[cat]); 
    }// finished with bkg samples separately, look at the combined sample
    // total bkg
    Bkg_Int[cat][fNBkg] = ABCDMethod::ComputeIntAndErr( fOutBkgTH2DHists[0], Bkg_IntErr[cat][fNBkg],  min_x[cat], max_x[cat], min_y[cat], max_y[cat]);
    // non-resonant bkg
    Bkg_Int[cat][fNBkg+1] = ABCDMethod::ComputeIntAndErr( fOutSelBkgTH2DHists[0], Bkg_IntErr[cat][fNBkg+1],  min_x[cat], max_x[cat], min_y[cat], max_y[cat]);

    for (UInt_t mc = 0; mc < fNSig; mc++){ 
      Sig_Int[cat][mc] = ABCDMethod::ComputeIntAndErr( fInSigTH2DHists[0][mc], Sig_IntErr[cat][mc],  min_x[cat], max_x[cat], min_y[cat], max_y[cat]); 
    } 

    std::cout << fBkgNames[i_gg] << " in " << cat << " int = " << Bkg_Int[cat][i_gg] << " pm " << Bkg_IntErr[cat][i_gg] << std::endl;
    //for (UInt_t mc = 0; mc < fNBkg+1; mc++){ std::cout << fBkgNames[mc] << " in " << cat << " reg "<< Bkg_Int[cat][mc] << " " << Bkg_IntErr[cat][mc] << std::endl; }
    //for (UInt_t mc = 0; mc < fNSig; mc++){   std::cout << fSigNames[mc] << " in " << cat << " reg "<< Sig_Int[cat][mc] << " " << Sig_IntErr[cat][mc] << std::endl; }

  }// end cat loop over A1,B1,A2,B2,C,D,All  

  ABCDMethod::GetFinalValuesForABCDReg(); // merge A1&A2->A and B1&B2->B
  for (UInt_t reg=0; reg <4; reg++){// print out A,B,D,C values
    std::cout << fBkgNames[i_gg] << " in " << reg << " int = " << fBkg_Int[reg][i_gg] << " pm " << fBkg_IntErr[reg][i_gg] << std::endl;
  //  for (UInt_t mc = 0; mc < fNBkg+1; mc++){ std::cout << fBkgNames[mc] << " " << reg << " has int " << fBkg_Int[reg][mc] << " pm " << fBkg_IntErr[reg][mc] << std::endl; }
  //  for (UInt_t mc = 0; mc < fNSig; mc++){   std::cout << fSigNames[mc] << " " << reg << " has int " << fSig_Int[reg][mc] << " pm " << fSig_IntErr[reg][mc] << std::endl; }
  }


  ABCDMethod::DoABCDCalculations(); // calculate corr & diff values

  fExpData.resize(1);
  fExpBkg.resize(fNBkg+2);
  fExpSig.resize(fNSig);
  fExpErrData.resize(1);
  fExpErrBkg.resize(fNBkg+2);
  fExpErrSig.resize(fNSig);

  fExpData[0]=ABCDMethod::FindExpectedValuesInD(fData_Int[0][0],fData_Int[1][0],fData_Int[3][0],fData_IntErr[0][0],fData_IntErr[1][0],fData_IntErr[3][0],fExpErrData[0]);
  //std::cout << "Data: Exp D = " << fExpData[0] << " Exp D err " << fExpErrData[0] << std::endl;
  for (UInt_t mc = 0; mc < fNBkg+2; mc++){
    fExpBkg[mc]=ABCDMethod::FindExpectedValuesInD(fBkg_Int[0][mc],fBkg_Int[1][mc],fBkg_Int[3][mc],fBkg_IntErr[0][mc],fBkg_IntErr[1][mc],fBkg_IntErr[3][mc],fExpErrBkg[mc]);
    if (mc == i_gg){
      std::cout << fBkgNames[mc] << ": nA = " << fBkg_Int[0][mc] << " nB = " << fBkg_Int[1][mc] << " nC = " << fBkg_Int[3][mc] << std::endl;
      std::cout << fBkgNames[mc] << ": Exp D = " << fExpBkg[mc] << " Exp D err " << fExpErrBkg[mc] << std::endl;
    }
  }
  for (UInt_t mc = 0; mc < fNSig; mc++){
    fExpSig[mc]=ABCDMethod::FindExpectedValuesInD(fSig_Int[0][mc],fSig_Int[1][mc],fSig_Int[3][mc],fSig_IntErr[0][mc],fSig_IntErr[1][mc],fSig_IntErr[3][mc],fExpErrSig[mc]);
    //std::cout << fSigNames[mc] <<": Exp D = " << fExpSig[mc] << " Exp D err " << fExpErrSig[mc] << std::endl;
  }

  ABCDMethod::FillTable();

  for (UInt_t mc = 0; mc < fNSig; mc++){
    ABCDMethod::WriteDataCard(fSigNames[mc].Data(),fRooSig[2][mc],fExpSig[mc],fBkg_Int,fRooBkg);
  }
}

void ABCDMethod::GetFinalValuesForABCDReg(){
 
  //combine A1&A2 to A region and B1&B2 to B region

  UInt_t fNReg = 5; // for 4 regions A,B,D,C, 5th is ALL
  fData_Int.resize(fNReg);
  fData_IntErr.resize(fNReg);
  fBkg_Int.resize(fNReg);
  fBkg_IntErr.resize(fNReg);
  fSig_Int.resize(fNReg);
  fSig_IntErr.resize(fNReg);
  for (UInt_t cat = 0; cat < fNReg; cat++){ // loop over each category
    fData_Int[cat].resize(1); 		// only one group for data since it is lumped together
    fData_IntErr[cat].resize(1);		
    fBkg_Int[cat].resize(fNBkg+2);	// do all Bkg separately and then one where all combined
    fBkg_IntErr[cat].resize(fNBkg+2);
    fSig_Int[cat].resize(fNSig);		// do all Sig separately
    fSig_IntErr[cat].resize(fNSig);

    if (cat == 0 || cat == 1){ // A1 or B1 region sum the integral and get error from sqrt(errReg1^2+errReg2^2)
      fData_Int[cat][0] = (Data_Int[cat][0]+Data_Int[cat+2][0]); //cat+2 is the corresponding A2 and B2 regions 
      fData_IntErr[cat][0] = (std::sqrt(Data_IntErr[cat][0]*Data_IntErr[cat][0]+ Data_IntErr[cat+2][0]*Data_IntErr[cat+2][0])); 
      for (UInt_t mc = 0; mc < fNSig; mc++){ 
        fSig_Int[cat][mc] = (Sig_Int[cat][mc]+Sig_Int[cat+2][mc]); 
        fSig_IntErr[cat][mc] = (std::sqrt(Sig_IntErr[cat][mc]*Sig_IntErr[cat][mc]+ Sig_IntErr[cat+2][mc]*Sig_IntErr[cat+2][mc])); 
      }
      for (UInt_t mc = 0; mc < fNBkg+2; mc++){ 
        fBkg_Int[cat][mc] = (Bkg_Int[cat][mc]+Bkg_Int[cat+2][mc]); 
        fBkg_IntErr[cat][mc] = (std::sqrt(Bkg_IntErr[cat][mc]*Bkg_IntErr[cat][mc]+ Bkg_IntErr[cat+2][mc]*Bkg_IntErr[cat+2][mc])); 
      }
    }
    else { // D,C, or ALL region, just take value from calculations above
      fData_Int[cat][0] = (Data_Int[cat+2][0]); //cat+2 is the corresponding C and D regions 
      fData_IntErr[cat][0] = (Data_IntErr[cat+2][0]); 
      for (UInt_t mc = 0; mc < fNSig; mc++){ 
        fSig_Int[cat][mc] = (Sig_Int[cat+2][mc]); 
        fSig_IntErr[cat][mc] = (Sig_IntErr[cat+2][mc]);
      }
      for (UInt_t mc = 0; mc < fNBkg+2; mc++){ 
        fBkg_Int[cat][mc] = (Bkg_Int[cat+2][mc]); 
        fBkg_IntErr[cat][mc] = (Bkg_IntErr[cat+2][mc]);
      }
    }  

  }// end cat loop over A,B,D,C
}

void ABCDMethod::DoABCDCalculations(){

  // calculate correlation for each sample
  fCorrData.push_back(fOutDataTH2DHists[0]->GetCorrelationFactor(1,2)); 
  for (UInt_t mc = 0; mc < fNBkg; mc++){ 
    fCorrBkg.push_back(fInBkgTH2DHists[0][mc]->GetCorrelationFactor(1,2));
  }
  fCorrBkg.push_back(fOutBkgTH2DHists[0]->GetCorrelationFactor(1,2)); // all bkg samples added together
  fCorrBkg.push_back(fOutSelBkgTH2DHists[0]->GetCorrelationFactor(1,2)); // all non-resonant bkg samples
  for (UInt_t mc = 0; mc < fNSig; mc++){ 
    fCorrSig.push_back(fInSigTH2DHists[0][mc]->GetCorrelationFactor(1,2));
  } 

  //calculate Diff = TMath::Abs((NC*NA/NB-ND)/(NC*NA/NB));
  fDiffData.resize(1);
  fDiffBkg.resize(fNBkg+2);
  fDiffSig.resize(fNSig);

  fDiffData[0] = ABCDMethod::FindDiff(fData_Int[0][0],fData_Int[1][0],fData_Int[3][0],fData_Int[2][0]);
  //std::cout << "Data Diff = " << fDiffData[0] << std::endl;
  for (UInt_t mc = 0; mc < fNBkg+2; mc++){
    fDiffBkg[mc]= ABCDMethod::FindDiff(fBkg_Int[0][mc],fBkg_Int[1][mc],fBkg_Int[3][mc],fBkg_Int[2][mc]); 
    //std::cout << fBkgNames[mc] << " Diff = " << fDiffBkg[mc] << std::endl;
  }
  for (UInt_t mc = 0; mc < fNSig; mc++){ 
    fDiffSig[mc]= ABCDMethod::FindDiff(fSig_Int[0][mc],fSig_Int[1][mc],fSig_Int[3][mc],fSig_Int[2][mc]); 
    //std::cout << fSigNames[mc] << " Diff = " << fDiffSig[mc] << std::endl;
  }

}

Double_t ABCDMethod::FindDiff(const Double_t NA, const Double_t NB, const Double_t NC, const Double_t ND){
  Double_t Diff = 0.;
  if ( NB-ND != 0 && NA > 0 && NB > 0 && NC > 0 ){
    Diff = TMath::Abs(((NC*NA/NB)-ND)/(NC*NA/NB));// This does: (Dexp-Dobs)/Dexp 
  }  
  return Diff;
}

Double_t ABCDMethod::FindExpectedValuesInD(const Double_t NA, const Double_t NB, const Double_t NC, const Double_t NAerr, const Double_t NBerr, const Double_t NCerr, Double_t & NDerr){ // find expected values in the D (signal) region
  Double_t ExpND = 0.;
  NDerr = 0.;
  if (NB > 0){ 
    ExpND = NC*NA/NB;                                                                                                                                
    NDerr = std::sqrt((NCerr*NCerr*NA*NA/(NB*NB))+(NAerr*NAerr*NC*NC/(NB*NB))+(NBerr*NBerr*NC*NC*NA*NA/(std::pow(NB,4))));
  }
  return ExpND;
}

void ABCDMethod::FillTable(){
  if (fOutTableTxtFile.is_open()){

     //setup Latex doc
     fOutTableTxtFile << "\\documentclass[a4paper,landscape]{article}" << std::endl;
     fOutTableTxtFile << "\\usepackage[paperheight=15.0in,paperwidth=6.0in,margin=1.0in,headheight=0.0in,footskip=0.5in,includehead,includefoot]{geometry}" << std::endl;
     fOutTableTxtFile << "\\begin{document}" << std::endl;
 
     // first table has the summary of total number of events for each sample
     fOutTableTxtFile << "\% Summary of Number of Events for each sample" << std::endl;
     fOutTableTxtFile << "\\begin{table}[bthp]" <<std::endl;
     fOutTableTxtFile << "\\begin{tabular}{cc}" <<std::endl;
     fOutTableTxtFile << "\\hline \\hline" <<std::endl;
     fOutTableTxtFile << Form("$\\sqrt{s}$ = 13 TeV; L = %3.1f $pb^{-1}$",lumi) <<" \\\\" <<std::endl;
     fOutTableTxtFile << "\\hline" <<std::endl;

     ABCDMethod::SetRooVariables();
     fOutTableTxtFile << "Data &  " << *(fRooData[4][0]->format(2,"EXPF")) <<  " \\\\" << std::endl;
     //std::cout << "Data &  " << *(fRooData[4][0]->format(2,"EXPF")) << std::endl;
     fOutTableTxtFile << "\\hline" << std::endl;

     TString name = "";
     for (UInt_t mc = 0; mc < fNBkg; mc++){
       fOutTableTxtFile << fSampleTitleMap[fBkgNames[mc]] << " &  " << *(fRooBkg[4][mc]->format(2,"EXPF")) <<  " \\\\" << std::endl;
       //std::cout << fBkgNames[mc].Data() <<  " &  " << *(fRooBkg[4][mc]->format(2,"EXPF")) << std::endl;
     }
     fOutTableTxtFile << "\\hline" << std::endl;
       fOutTableTxtFile << "Total Bkg &  " << *(fRooBkg[4][fNBkg]->format(2,"EXPF")) <<  " \\\\" << std::endl;
     fOutTableTxtFile << "\\hline" << std::endl;
     
     for (UInt_t mc = 0; mc < fNSig; mc++){
       fOutTableTxtFile << fSampleTitleMap[fSigNames[mc]] << " &  " << *(fRooSig[4][mc]->format(2,"EXPF")) <<  " \\\\" << std::endl; 
       //std::cout << fSigNames[mc] <<  " &  " << *(fRooSig[4][mc]->format(2,"EXPF")) << std::endl;
     }
 
     fOutTableTxtFile << "\\hline \\hline" <<std::endl;
     fOutTableTxtFile << "\\end{tabular}" <<std::endl;
     fOutTableTxtFile << "\\end{table}" <<std::endl;
     // end first table

     // second table has the summary of total number of events in signal region for each sample
     fOutTableTxtFile << "\% Summary of Events in Signal Region for each sample" << std::endl;
     fOutTableTxtFile << "\\begin{table}[bthp]" <<std::endl;
     fOutTableTxtFile << "\\begin{tabular}{cc}" <<std::endl;
     fOutTableTxtFile << "\\hline \\hline" <<std::endl;
     fOutTableTxtFile << Form("$\\sqrt{s}$ = 13 TeV; L = %3.1f $pb^{-1}$",lumi) <<" \\\\" <<std::endl;
     fOutTableTxtFile << Form("$m_{\\gamma \\gamma}$ in $[%3.f-%3.f]$ GeV and MET $>$ %3.f GeV",mgg_minCD,mgg_maxCD,met_minD)<<" \\\\" <<std::endl;
     fOutTableTxtFile << "\\hline" <<std::endl;

     ABCDMethod::SetRooVariables();
     fOutTableTxtFile << "Data &  " << *(fRooData[2][0]->format(2,"EXPF")) <<  " \\\\" << std::endl;
     //std::cout << "Data &  " << *(fRooData[2][0]->format(2,"EXPF")) << std::endl;
     fOutTableTxtFile << "\\hline" << std::endl;

     for (UInt_t mc = 0; mc < fNBkg; mc++){
       fOutTableTxtFile << fSampleTitleMap[fBkgNames[mc]] << " &  " << *(fRooBkg[2][mc]->format(2,"EXPF")) <<  " \\\\" << std::endl;
       //std::cout << fBkgNames[mc].Data() <<  " &  " << *(fRooBkg[2][mc]->format(2,"EXPF")) << std::endl;
     }
     fOutTableTxtFile << "\\hline" << std::endl;
       fOutTableTxtFile << "Total Bkg &  " << *(fRooBkg[2][fNBkg]->format(2,"EXPF")) <<  " \\\\" << std::endl;
     fOutTableTxtFile << "\\hline" << std::endl;
     
     for (UInt_t mc = 0; mc < fNSig; mc++){
       fOutTableTxtFile << fSampleTitleMap[fSigNames[mc]] << " &  " << *(fRooSig[2][mc]->format(2,"EXPF")) <<  " \\\\" << std::endl; 
       //std::cout << fSigNames[mc] <<  " &  " << *(fRooSig[2][mc]->format(2,"EXPF")) << std::endl;
     }
 
     fOutTableTxtFile << "\\hline \\hline" <<std::endl;
     fOutTableTxtFile << "\\end{tabular}" <<std::endl;
     fOutTableTxtFile << "\\end{table}" <<std::endl;
     // end second table
    
     
     // third table has values for A,B,C,ExpD,D,Diff,Corr 

     fOutTableTxtFile << "\% Table with values of A,B,C,ExpD,D,Diff,Corr for each sample" << std::endl;
     fOutTableTxtFile << "\\begin{table}[bthp]" <<std::endl;
     fOutTableTxtFile << "\\begin{tabular}{ccccccccc}" <<std::endl;
     fOutTableTxtFile << "\\hline \\hline" <<std::endl;
     fOutTableTxtFile << Form("$\\sqrt{s}$ = 13 TeV; L = %3.1f $pb^{-1}$",lumi) <<" \\\\" <<std::endl;
     fOutTableTxtFile << "\\hline" <<std::endl;
     fOutTableTxtFile << "Sample & Corr($m_{\\gamma\\gamma}$,MET) & $N_A$ & $N_B$ & $N_C$ & $N_D$ & Expected D & Diff  \\\\" << std::endl;
     fOutTableTxtFile << "\\hline" <<std::endl;

     fOutTableTxtFile << "Data &  " << 
	*(fRooData[5][0]->format(2,"XPF")) << " &  " <<
        *(fRooData[0][0]->format(2,"EXPF")) << " &  " << 
        *(fRooData[1][0]->format(2,"EXPF")) << " &  " << 
        *(fRooData[3][0]->format(2,"EXPF")) << " &  " <<  
        *(fRooData[2][0]->format(2,"EXPF")) << " &  " << 
	*(fRooData[7][0]->format(2,"EXPF")) << " &  $" <<
        Form("%f",fDiffData[0]) <<"$ \\\\" << std::endl;
     fOutTableTxtFile << "\\hline" << std::endl;

     for (UInt_t mc = 0; mc < fNBkg; mc++){
       fOutTableTxtFile << fSampleTitleMap[fBkgNames[mc]] << " &  " << 
         *(fRooBkg[5][mc]->format(2,"XPF")) << " &  " << 
         *(fRooBkg[0][mc]->format(2,"EXPF")) << " &  " << 
         *(fRooBkg[1][mc]->format(2,"EXPF")) << " &  " << 
         *(fRooBkg[3][mc]->format(2,"EXPF")) << " &  " <<  
         *(fRooBkg[2][mc]->format(2,"EXPF")) << " &  " <<
         *(fRooBkg[7][mc]->format(2,"EXPF")) << " &  $" <<
         Form("%f",fDiffBkg[mc]) <<"$ \\\\" << std::endl;
     }
     fOutTableTxtFile << "\\hline" << std::endl;
       fOutTableTxtFile << "NonRes Bkg" << " &  " << // only non-resonant bkg here 
         *(fRooBkg[5][fNBkg+1]->format(2,"XPF")) << " &  " << 
         *(fRooBkg[0][fNBkg+1]->format(2,"EXPF")) << " &  " << 
         *(fRooBkg[1][fNBkg+1]->format(2,"EXPF")) << " &  " << 
         *(fRooBkg[3][fNBkg+1]->format(2,"EXPF")) << " &  " <<  
         *(fRooBkg[2][fNBkg+1]->format(2,"EXPF")) << " &  " << 
         *(fRooBkg[7][fNBkg+1]->format(2,"EXPF")) << " &  $" << 
         Form("%f",fDiffBkg[fNBkg+1]) <<"$ \\\\" << std::endl;
     fOutTableTxtFile << "\\hline" << std::endl;
     for (UInt_t mc = 0; mc < fNSig; mc++){
       fOutTableTxtFile << fSampleTitleMap[fSigNames[mc]] << " &  " << 
         *(fRooSig[5][mc]->format(2,"XPF")) << " &  " << 
         *(fRooSig[0][mc]->format(2,"EXPF")) << " &  " << 
         *(fRooSig[1][mc]->format(2,"EXPF")) << " &  " << 
         *(fRooSig[3][mc]->format(2,"EXPF")) << " &  " <<  
         *(fRooSig[2][mc]->format(2,"EXPF")) << " &  " << 
         *(fRooSig[7][mc]->format(2,"EXPF")) << " &  $" << 
         Form("%f",fDiffSig[mc]) <<"$ \\\\" << std::endl;
     }
 
     fOutTableTxtFile << "\\hline \\hline" <<std::endl;
     fOutTableTxtFile << "\\end{tabular}" <<std::endl;
     fOutTableTxtFile << "\\end{table}" <<std::endl;
     // end third table

     // finish Latex doc
     fOutTableTxtFile << "\\end{document}" <<std::endl;
     std::cout << "Writing ResultsTable in " << Form("%s/ResultsTableForLatex.tex",fOutDir.Data()) << std::endl;

     //std::cout << "Data: A    =  " << *(fRooData[0][0]->format(2,"EXPF")) << std::endl;
     //std::cout << "Data: B    =  " << *(fRooData[1][0]->format(2,"EXPF")) << std::endl;
     //std::cout << "Data: C    =  " << *(fRooData[3][0]->format(2,"EXPF")) << std::endl;
     //std::cout << "Data: D    =  " << *(fRooData[2][0]->format(2,"EXPF")) << std::endl;
     //std::cout << "Data: ExpD =  " << *(fRooData[7][0]->format(2,"EXPF")) << std::endl;
     //std::cout << "Data: corr =  " << *(fRooData[5][0]->format(2,"EXPF")) << std::endl;
  
     //TString bkgname="";
     //for (UInt_t mc = 0; mc < fNBkg+2; mc++){
     //  if (mc == fNBkg) bkgname="Total Bkg";
     //  else if (mc==fNBkg+1) bkgname="Non-res Bkg";
     //  else bkgname = fBkgNames[mc]; 
     //  std::cout << bkgname << ": A    =  " << *(fRooBkg[0][mc]->format(2,"EXPF")) << std::endl;
     //  std::cout << bkgname << ": B    =  " << *(fRooBkg[1][mc]->format(2,"EXPF")) << std::endl;
     //  std::cout << bkgname << ": C    =  " << *(fRooBkg[3][mc]->format(2,"EXPF")) << std::endl;
     //  std::cout << bkgname << ": D    =  " << *(fRooBkg[2][mc]->format(2,"EXPF")) << std::endl;
     //  std::cout << bkgname << ": ExpD =  " << *(fRooBkg[5][mc]->format(2,"EXPF")) << std::endl;
     //  std::cout << bkgname << ": corr =  " << *(fRooBkg[7][mc]->format(2,"EXPF")) << std::endl;
     //}
 
     //for (UInt_t mc = 0; mc < fNSig; mc++){
     //  std::cout << fSigNames[mc] << ": A    =  " << *(fRooSig[0][mc]->format(2,"EXPF")) << std::endl;
     //  std::cout << fSigNames[mc] << ": B    =  " << *(fRooSig[1][mc]->format(2,"EXPF")) << std::endl;
     //  std::cout << fSigNames[mc] << ": C    =  " << *(fRooSig[3][mc]->format(2,"EXPF")) << std::endl;
     //  std::cout << fSigNames[mc] << ": D    =  " << *(fRooSig[2][mc]->format(2,"EXPF")) << std::endl;
     //  std::cout << fSigNames[mc] << ": ExpD =  " << *(fRooSig[7][mc]->format(2,"EXPF")) << std::endl;
     //  std::cout << fSigNames[mc] << ": corr =  " << *(fRooSig[5][mc]->format(2,"EXPF")) << std::endl;
     //}
  }
  else std::cout << "Unable to open ResultsTable Output File" <<std::endl;

}

void ABCDMethod::SetRooVariables(){

  UInt_t fNReg = 8;
  fRooData.resize(fNReg);
  fRooBkg.resize(fNReg);
  fRooSig.resize(fNReg);
  TString cat_name = "";
  TString name = "";
  Double_t multData,multSig,multBkg;

  for (UInt_t cat = 0; cat < fNReg; cat++){// loop over A,B,C,D,ExpD,Corr
    fRooData[cat].resize(1);
    fRooBkg[cat].resize(fNBkg+2);
    fRooSig[cat].resize(fNSig);

    if (cat==0) cat_name = "_A";
    if (cat==1) cat_name = "_B";
    if (cat==2) cat_name = "_D";
    if (cat==3) cat_name = "_C";
    if (cat==4) cat_name = "_All";
    if (cat==5) cat_name = "_Corr";
    if (cat==6) cat_name = "_mult";
    if (cat==7) cat_name = "_Exp";
 
    name = Form("Data%s",cat_name.Data());
    if (cat < 5){
      fRooData[cat][0] = new RooRealVar(name,name,fData_Int[cat][0]); 
      fRooData[cat][0]->setError(fData_IntErr[cat][0]);
    }
    if (cat == 5){
      fRooData[cat][0] = new RooRealVar(name,name,fCorrData[0]);
    } 
    if (cat == 6){
      multData = fData_Int[0][0]/fData_Int[1][0];
      fRooData[cat][0] = new RooRealVar(name,name,multData);
    } 
    if (cat == 7){
      fRooData[cat][0] = new RooRealVar(name,name,fExpData[0]);
      fRooData[cat][0]->setError(fExpErrData[0]);
    }

    for (UInt_t mc = 0; mc < fNBkg+2; mc++){
      if (mc<fNBkg) name = Form("%s%s",fBkgNames[mc].Data(),cat_name.Data());
      else if (mc == fNBkg) name = Form("TotBkg%s",cat_name.Data());
      else name = Form("SelBkg%s",cat_name.Data());
      if (cat < 5){
        fRooBkg[cat][mc] = new RooRealVar(name,name,fBkg_Int[cat][mc]);
        fRooBkg[cat][mc]->setError(fBkg_IntErr[cat][mc]);
      }
      if (cat == 5){
        fRooBkg[cat][mc] = new RooRealVar(name,name,fCorrBkg[mc]);
      } 
      if (cat == 6){
        multBkg = fBkg_Int[0][mc]/fBkg_Int[1][mc];
        fRooBkg[cat][mc] = new RooRealVar(name,name,multBkg);
      } 
      if (cat == 7){
        fRooBkg[cat][mc] = new RooRealVar(name,name,fExpBkg[mc]);
        fRooBkg[cat][mc]->setError(fExpErrBkg[mc]);
      }
    }// end loop over Bkg mc
   
    for (UInt_t mc = 0; mc < fNSig; mc++){
      name = Form("%s%s",fSigNames[mc].Data(),cat_name.Data());
      if (cat < 5){
        fRooSig[cat][mc] = new RooRealVar(name,name,fSig_Int[cat][mc]);
        fRooSig[cat][mc]->setError(fSig_IntErr[cat][mc]);
      }
      if (cat == 5){
        fRooSig[cat][mc] = new RooRealVar(name,name,fCorrSig[mc]);
      } 
      if (cat == 6){
        multSig = fSig_Int[0][mc]/fSig_Int[1][mc];
        fRooSig[cat][mc] = new RooRealVar(name,name,multSig);
      } 
      if (cat == 7){
        fRooSig[cat][mc] = new RooRealVar(name,name,fExpSig[mc]);
        fRooSig[cat][mc]->setError(fExpErrSig[mc]);
      }
    }// end loop over Sig mc
  }// end loop over categories

}// end SetRooVariables

void ABCDMethod::WriteDataCard( const TString fSigName, const RooRealVar* sigrate, const Double_t expsig, const DblVecVec bkgrates, const RooVecVec bkgrate){
  TString sig = *sigrate->format(2,"");
  TString vh  = *bkgrate[2][i_vh]->format(2,"");
  TString hgg = *bkgrate[2][i_hgg]->format(2,"");
  TString dy  = *bkgrate[7][i_dy]->format(2,"");
  TString gg  = *bkgrate[7][i_gg]->format(2,"");
  TString qcd = *bkgrate[7][i_qcd]->format(2,""); 
  TString gj  = *bkgrate[7][i_gj]->format(2,"");
  //std::cout << "sig = " << sig << " vh " << vh << " hgg " << hgg << " dy " << dy << " gg " << gg << " qcd " << qcd << " gj " << gj << std::endl; 

  TStrVec N_C;
  N_C.resize(fNBkg);
  for (UInt_t mc = 0; mc < fNBkg; mc++){
    N_C[mc] = *bkgrate[3][mc]->format(2,"");
    //std::cout << fBkgNames[mc] << ": C = " << N_C[mc] <<std::endl;
  }

  TString fac_vh  = *bkgrate[6][i_vh]->format(2,"");
  TString fac_hgg = *bkgrate[6][i_hgg]->format(2,"");
  TString fac_dy  = *bkgrate[6][i_dy]->format(2,"");
  TString fac_gg  = *bkgrate[6][i_gg]->format(2,"");
  TString fac_qcd = *bkgrate[6][i_qcd]->format(2,""); 
  TString fac_gj  = *bkgrate[6][i_gj]->format(2,"");


  DblVec N_A,N_B,mult; // N_C,mult;
  N_A.resize(fNBkg);
  N_B.resize(fNBkg);
  mult.resize(fNBkg);

  for (UInt_t mc=0; mc < fNBkg; mc++){
    N_A[mc] = bkgrates[0][mc];
    N_B[mc] = bkgrates[1][mc];
    mult[mc]= N_A[mc]/N_B[mc];
    //std::cout << fBkgNames[mc].Data() << " nA = " << N_A[mc] << " nB = " << N_B[mc] << " nA/nB = " << mult[mc] << std::endl;
  }
 
  std::cout << "Writing data card in: " << fOutDir.Data() << "/DataCard_" << fSigName.Data() <<".txt" << std::endl;
  fOutTxtFile.open(Form("%s/DataCard_%s.txt",fOutDir.Data(),fSigName.Data())); 
  // print out the Data Card file
  if (fOutTxtFile.is_open()){
    fOutTxtFile << Form("#MonoHgg DataCard for C&C Limit Setting, %f pb-1 ",lumi) << std::endl;
    fOutTxtFile << "#Run with:combine -M Asymptotic cardname.txt --run blind " << std::endl;
    fOutTxtFile << Form("# Lumi =  %f pb-1",lumi) << std::endl;
    fOutTxtFile << "imax 1" << std::endl;
    fOutTxtFile << "jmax *" << std::endl;
    fOutTxtFile << "kmax *" << std::endl;
    fOutTxtFile << "---------------" << std::endl;
 
    fOutTxtFile << "bin 1"<< std::endl;
    fOutTxtFile << "observation 0 "  << std::endl;
    fOutTxtFile << "------------------------------" << std::endl;
    fOutTxtFile << "bin		1		1		1		1		1		1		1"<< std::endl;
    fOutTxtFile << "process		DM		gg		dy		qcd		gj		hgg		vh" << std::endl;
    fOutTxtFile << "process		0		1		2		3		4		5 		6" << std::endl;
    fOutTxtFile << Form("rate	       %s	       %s	       %s	       %s	       %s	       %s	       %s",sig.Data(),gg.Data(),dy.Data(),qcd.Data(),gj.Data(),hgg.Data(),vh.Data()) << std::endl; 
    fOutTxtFile << "--------------------------------" << std::endl;
    fOutTxtFile << "#signal related" << std::endl; //just took these numbers from Livia's example (all estimates from 8TeV)
    fOutTxtFile << "lumi_13TeV lnN	1.1000		-		-		-		-		1.1000     	1.1000" << std::endl;
    fOutTxtFile << "eff_trig   lnN	1.010000	-		-		-		-		1.01000    	1.01000" << std::endl;
    fOutTxtFile << "id_eff_eb  lnN	1.02000		-		-		-		-		1.02000    	1.02000   " << std::endl;    
    fOutTxtFile << "vtxEff     lnN	0.996/1.008	-		-		-		-		0.996/1.008	0.996/1.008" << std::endl; 
    fOutTxtFile << "#background related" << std::endl;
    //fOutTxtFile << "abcd_estimate lnN	-	1.27000		1.27000	   1.27000	1.27000		-	-	-  " << std::endl;
    fOutTxtFile << Form("gg_norm  gmN 	%s	-	%s	   -		-		-	-	-  ",N_C[i_gg].Data(),fac_gg.Data()) << std::endl;
    fOutTxtFile << Form("dy_norm  gmN 	%s	-	-	   %s		-		-	-	-  ",N_C[i_dy].Data(),fac_dy.Data()) << std::endl;
    fOutTxtFile << Form("qcd_norm gmN 	%s	-	- 	   -		%s		-	-	-  ",N_C[i_qcd].Data(),fac_qcd.Data()) << std::endl;
    fOutTxtFile << Form("gj_norm  gmN 	%s	-	-	   -		-		%s	-	-  ",N_C[i_gj].Data(),fac_gj.Data()) << std::endl;

  }
  else std::cout << "Unable to open DataCard Output File" << std::endl;
 
  fOutTxtFile.close();
  std::cout << "Finished Writing DataCard" << std::endl;
}

//void ABCDMethod::FindBins(TH2D *& h, con

Double_t ABCDMethod::ComputeIntAndErr(TH2D *& h, Double_t & error, const UInt_t minX, const UInt_t maxX, const UInt_t minY, const UInt_t maxY ){
  Double_t integral = 0.;
  if(h == (TH2D*) NULL) std::cout << "NULL TH2D" << std::endl;
  //std::cout << isReg <<  " minx = " << minX << " maxX = " << maxX << " minY = " << minY << " maxY = " << maxY << std::endl; 
  integral = h->IntegralAndError(minX,maxX,minY,maxY,error);
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
  fOutSelBkgTH2DHists.resize(fNTH2D);
}

void ABCDMethod::InitVariables(){
  // 1D histograms of interest
  fTH1DNames.push_back("mgg");
  fTH1DNames.push_back("t1pfmet");

  // 2D histograms of interest
  fTH2DNames.push_back("t1pfmet_mgg");
}

