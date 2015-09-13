#include "ABCDMethod.hh"

ABCDMethod::ABCDMethod( SamplePairVec Samples, const Double_t inLumi, const TString outdir){

  // load RooFit
  gSystem->Load("libRooFit");
  //using namespace RooFit;

  lumi = inLumi;
  fInDir = outdir;
  fOutDir = outdir+"ABCD";

  // values for the different bin edges
  mgg_minAB1 = 100.;
  mgg_minCD  = 110.;
  mgg_maxCD  = 130.;
  mgg_maxAB2 = 180.; 
  met_minB   = 150.;
  met_minD   = 250.;
  met_maxD   = 400.;

  // make output txt file with output table
  fOutTableTxtFile.open(Form("%s/ResultsTableForLatex.txt",fOutDir.Data()));

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

  // delete RooRealVar
  delete fRData[0];
  for (UInt_t mc = 0; mc < fNBkg+1; mc++){
    delete fRBkg[mc];
  }
  for (UInt_t mc = 0; mc < fNSig; mc++){
    delete fRSig[mc];
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
  for (UInt_t mc = 0; mc < fNBkg; mc++){
    fInBkgTH2DHists[0][mc]->Scale(lumi);
    //std::cout << "number entries in bkg in " << fInBkgTH2DHists[0][mc]->GetEntries() << std::endl;
    if (mc == 0) fOutBkgTH2DHists[0] = (TH2D*)fInBkgTH2DHists[0][mc]->Clone();
    else fOutBkgTH2DHists[0]->Add(fInBkgTH2DHists[0][mc]);
  } 
  //std::cout << "number entries in bkg " << fOutBkgTH2DHists[0]->GetEntries() << std::endl;

  // just scale the signal by lumi, don't add together 
  for (UInt_t mc = 0; mc < fNSig; mc++){
    fInSigTH2DHists[0][mc]->Scale(lumi);
    //std::cout << "number entries in sig in " << fInSigTH2DHists[0][mc]->GetEntries() << std::endl;
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

    Data_Int[cat][0] = ABCDMethod::ComputeIntAndErr( fOutDataTH2DHists[0], Data_IntErr[cat][0],  min_x[cat], max_x[cat], min_y[cat], max_y[cat], cat); 
    for (UInt_t mc = 0; mc < fNBkg; mc++){ 
      Bkg_Int[cat][mc] = ABCDMethod::ComputeIntAndErr( fInBkgTH2DHists[0][mc], Bkg_IntErr[cat][mc],  min_x[cat], max_x[cat], min_y[cat], max_y[cat], cat); 
    }
    // after finished with bkg samples separately, look at the combined sample
    Bkg_Int[cat][fNBkg] = ABCDMethod::ComputeIntAndErr( fOutBkgTH2DHists[0], Bkg_IntErr[cat][fNBkg],  min_x[cat], max_x[cat], min_y[cat], max_y[cat], cat);
    for (UInt_t mc = 0; mc < fNSig; mc++){ 
      Sig_Int[cat][mc] = ABCDMethod::ComputeIntAndErr( fInSigTH2DHists[0][mc], Sig_IntErr[cat][mc],  min_x[cat], max_x[cat], min_y[cat], max_y[cat], cat); 
    } 

    //std::cout << "Data " << Data_Int[cat][0] << " " << Data_IntErr[cat][0] << std::endl;
    //for (UInt_t mc = 0; mc < fNBkg+1; mc++){ std::cout << fBkgNames[mc] << " in " << cat << " reg "<< Bkg_Int[cat][mc] << " " << Bkg_IntErr[cat][mc] << std::endl; }
    //for (UInt_t mc = 0; mc < fNSig; mc++){   std::cout << fSigNames[mc] << " in " << cat << " reg "<< Sig_Int[cat][mc] << " " << Sig_IntErr[cat][mc] << std::endl; }

  }// end cat loop over A1,B1,A2,B2,C,D  

  // calculate integral and error for full range
  fFullData_Int.resize(1);
  fFullBkg_Int.resize(fNBkg+1);
  fFullSig_Int.resize(fNSig);
  fFullData_IntErr.resize(1);
  fFullBkg_IntErr.resize(fNBkg+1);
  fFullSig_IntErr.resize(fNSig);
 
  fFullData_Int[0] = ABCDMethod::ComputeIntAndErr( fOutDataTH2DHists[0], fFullData_IntErr[0], mgg_minAB1, mgg_maxAB2, met_minB, met_maxD, fNCat); 
  for (UInt_t mc = 0; mc < fNBkg; mc++){ 
    fFullBkg_Int[mc] = ABCDMethod::ComputeIntAndErr( fInBkgTH2DHists[0][mc], fFullBkg_IntErr[mc], mgg_minAB1, mgg_maxAB2, met_minB, met_maxD, fNCat); 
  }
  fFullBkg_Int[fNBkg] = ABCDMethod::ComputeIntAndErr( fOutBkgTH2DHists[0], fFullBkg_IntErr[fNBkg], mgg_minAB1, mgg_maxAB2, met_minB, met_maxD, fNCat); 
  for (UInt_t mc = 0; mc < fNSig; mc++){ 
    fFullSig_Int[mc] = ABCDMethod::ComputeIntAndErr( fInSigTH2DHists[0][mc], fFullSig_IntErr[mc], mgg_minAB1, mgg_maxAB2, met_minB, met_maxD, fNCat); 
  } 

  //std::cout << "All Data int " << fFullData_Int[0] << " error " << fFullData_IntErr[0] << std::endl;
  //for (UInt_t mc = 0; mc < fNBkg+1; mc++){ std::cout << "All " << fBkgNames[mc] << " int " << fFullBkg_Int[mc] << " " << fFullBkg_IntErr[mc] << std::endl; }
  //for (UInt_t mc = 0; mc < fNSig; mc++){   std::cout << "All " << fSigNames[mc] << " int " << fFullSig_Int[mc] << " " << fFullSig_IntErr[mc] << std::endl; }


  ABCDMethod::GetFinalValuesForABCDReg(); // merge A1&A2->A and B1&B2->B
  //for (UInt_t reg=0; reg <4; reg++){// print out A,B,D,C values
  //  std::cout << "Data in reg " << reg << " has int = " << fData_Int[reg][0] << " and err = " << fData_IntErr[reg][0] << std::endl;
  //  for (UInt_t mc = 0; mc < fNBkg+1; mc++){ std::cout << fBkgNames[mc] << " in reg " << reg << " has int " << fBkg_Int[reg][mc] << " and err " << fBkg_IntErr[reg][mc] << std::endl; }
  //  for (UInt_t mc = 0; mc < fNSig; mc++){   std::cout << fSigNames[mc] << " in reg " << reg << " has int " << fSig_Int[reg][mc] << " and err " << fSig_IntErr[reg][mc] << std::endl; }
  //}


  ABCDMethod::DoABCDCalculations(); // calculate corr & diff values

  fExpData.resize(1);
  fExpBkg.resize(fNBkg+1);
  fExpSig.resize(fNSig);
  fExpErrData.resize(1);
  fExpErrBkg.resize(fNBkg+1);
  fExpErrSig.resize(fNSig);

  fExpData[0]=ABCDMethod::FindExpectedValuesInD(fData_Int[0][0],fData_Int[1][0],fData_Int[3][0],fData_IntErr[0][0],fData_IntErr[1][0],fData_IntErr[3][0],fExpErrData[0]);
  std::cout << "Data: Exp D = " << fExpData[0] << " Exp D err " << fExpErrData[0] << std::endl;
  for (UInt_t mc = 0; mc < fNBkg+1; mc++){
    fExpBkg[mc]=ABCDMethod::FindExpectedValuesInD(fBkg_Int[0][mc],fBkg_Int[1][mc],fBkg_Int[3][mc],fBkg_IntErr[0][mc],fBkg_IntErr[1][mc],fBkg_IntErr[3][mc],fExpErrBkg[mc]);
    std::cout << fBkgNames[mc] << ": Exp D = " << fExpBkg[mc] << " Exp D err " << fExpErrBkg[mc] << std::endl;
  }
  for (UInt_t mc = 0; mc < fNSig; mc++){
    fExpSig[mc]=ABCDMethod::FindExpectedValuesInD(fSig_Int[0][mc],fSig_Int[1][mc],fSig_Int[3][mc],fSig_IntErr[0][mc],fSig_IntErr[1][mc],fSig_IntErr[3][mc],fExpErrSig[mc]);
    std::cout << fSigNames[mc] <<": Exp D = " << fExpSig[mc] << " Exp D err " << fExpErrSig[mc] << std::endl;
  }

  ABCDMethod::FillTable();

  for (UInt_t mc = 0; mc < fNSig; mc++){
    ABCDMethod::WriteDataCard(fSigNames[mc].Data());
  }
}

void ABCDMethod::GetFinalValuesForABCDReg(){
 
  //combine A1&A2 to A region and B1&B2 to B region

  UInt_t fNReg = 4; // for 4 regions A,B,D,C
  fData_Int.resize(fNReg);
  fData_IntErr.resize(fNReg);
  fBkg_Int.resize(fNReg);
  fBkg_IntErr.resize(fNReg);
  fSig_Int.resize(fNReg);
  fSig_IntErr.resize(fNReg);
  for (UInt_t cat = 0; cat < fNReg; cat++){ // loop over each category
    fData_Int[cat].resize(1); 		// only one group for data since it is lumped together
    fData_IntErr[cat].resize(1);		
    fBkg_Int[cat].resize(fNBkg+1);	// do all Bkg separately and then one where all combined
    fBkg_IntErr[cat].resize(fNBkg+1);
    fSig_Int[cat].resize(fNSig);		// do all Sig separately
    fSig_IntErr[cat].resize(fNSig);

    if (cat == 0 || cat == 1){ // A or B region sum the integral and get error from sqrt(errReg1^2+errReg2^2)
      fData_Int[cat][0] = (Data_Int[cat][0]+Data_Int[cat+2][0]); //cat+2 is the corresponding C and D regions 
      fData_IntErr[cat][0] = (std::sqrt(Data_IntErr[cat][0]*Data_IntErr[cat][0]+ Data_IntErr[cat+2][0]*Data_IntErr[cat+2][0])); 
      for (UInt_t mc = 0; mc < fNSig; mc++){ 
        fSig_Int[cat][mc] = (Sig_Int[cat][mc]+Sig_Int[cat+2][mc]); 
        fSig_IntErr[cat][mc] = (std::sqrt(Sig_IntErr[cat][mc]*Sig_IntErr[cat][mc]+ Sig_IntErr[cat+2][mc]*Sig_IntErr[cat+2][mc])); 
      }
      for (UInt_t mc = 0; mc < fNBkg+1; mc++){ 
        fBkg_Int[cat][mc] = (Bkg_Int[cat][mc]+Bkg_Int[cat+2][mc]); 
        fBkg_IntErr[cat][mc] = (std::sqrt(Bkg_IntErr[cat][mc]*Bkg_IntErr[cat][mc]+ Bkg_IntErr[cat+2][mc]*Bkg_IntErr[cat+2][mc])); 
      }
    }
    else{ // D or C region, just take value from calculations above
      fData_Int[cat][0] = (Data_Int[cat+2][0]); //cat+2 is the corresponding C and D regions 
      fData_IntErr[cat][0] = (Data_IntErr[cat+2][0]); 
      for (UInt_t mc = 0; mc < fNSig; mc++){ 
        fSig_Int[cat][mc] = (Sig_Int[cat+2][mc]); 
        fSig_IntErr[cat][mc] = (Sig_IntErr[cat+2][mc]);
      }
      for (UInt_t mc = 0; mc < fNBkg+1; mc++){ 
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
  for (UInt_t mc = 0; mc < fNSig; mc++){ 
    fCorrSig.push_back(fInSigTH2DHists[0][mc]->GetCorrelationFactor(1,2));
  } 

  //calculate Diff = TMath::Abs((NC*NA/NB-ND)/(NC*NA/NB));
  fDiffData.resize(1);
  fDiffBkg.resize(fNBkg+1);
  fDiffSig.resize(fNSig);

  fDiffData[0] = ABCDMethod::FindDiff(fData_Int[0][0],fData_Int[1][0],fData_Int[3][0],fData_Int[2][0]);
  //std::cout << "Data Diff = " << fDiffData[0] << std::endl;
  for (UInt_t mc = 0; mc < fNBkg+1; mc++){
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
    Diff = TMath::Abs((NC*NA/NB-ND)/(NC*NA/NB)); 
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

     fOutTableTxtFile << "\\begin{table}[bthp]" <<std::endl;
     fOutTableTxtFile << "\\begin{tabular}{cc}" <<std::endl;
     fOutTableTxtFile << "\\hline \\hline" <<std::endl;
     fOutTableTxtFile << Form("$\\sqrt{s}$ = 13 TeV; L = %f $fb^{-1}$",lumi) <<" \\\\" <<std::endl;
     fOutTableTxtFile << "$m_{\\gamma \\gamma}$ in $[110-130]$ and MET $>$ 250 GeV"<<" \\\\" <<std::endl;
     fOutTableTxtFile << "\\hline" <<std::endl;

     ABCDMethod::SetRooVariables();
     fOutTableTxtFile << "Data &  " << *(fRData[0]->format(2,"EXPF")) << std::endl;
     std::cout << "Data &  " << *(fRData[0]->format(2,"EXPF")) << std::endl;

     TString name = "";
     for (UInt_t mc = 0; mc < fNBkg; mc++){
       fOutTableTxtFile << fBkgNames[mc].Data() << " &  " << *(fRBkg[mc]->format(2,"EXPP")) <<  " \\\\" << std::endl;
       std::cout << fBkgNames[mc].Data() <<  " &  " << *(fRBkg[mc]->format(2,"EXPF")) << std::endl;
     }
     fOutTableTxtFile << "\\hline" << std::endl;
       fOutTableTxtFile << "Total Bkg &  " << *(fRBkg[fNBkg]->format(2,"EXPF")) <<  " \\\\" << std::endl;

     for (UInt_t mc = 0; mc < fNSig; mc++){
       fOutTableTxtFile << fSigNames[mc].Data() << " &  " << *(fRSig[mc]->format(2,"EXPF")) <<  " \\\\" << std::endl; 
       std::cout << fSigNames[mc] <<  " &  " << *(fRSig[mc]->format(2,"EXPF")) << std::endl;
     }
 
     fOutTableTxtFile << "\\hline \\hline" <<std::endl;
     fOutTableTxtFile << "\\end{tabular}" <<std::endl;
     fOutTableTxtFile << "\\end{table}" <<std::endl;

/*  double ppExpND = ppNC*ppNA/ppNB;
    double ppExpNDerr = sqrt((ppNCerr*ppNCerr*ppNA*ppNA/(ppNB*ppNB))+(ppNAerr*ppNAerr*ppNC*ppNC/(ppNB*ppNB))+(ppNBerr*ppNBerr*ppNC*ppNC*ppNA*ppNA/(pow(ppNB,4))));
    double ppDiff = TMath::Abs((ppNC*ppNA/ppNB-ppND)/(ppNC*ppNA/ppNB));
    RooRealVar ppA("ppA","ppA",ppNA,""); 	//ppNA is integral pp in A region
    ppA.setError(ppNAerr);			//ppNAerr is error pp in A region
    RooRealVar ppB("ppB","ppB",ppNB,"");
    ppB.setError(ppNBerr);
    RooRealVar ppC("ppC","ppC",ppNC,"");
    ppC.setError(ppNCerr);
    RooRealVar ppD("ppD","ppD",ppND,"");
    ppD.setError(ppNDerr);
    RooRealVar ppDexp("ppDexp","ppDexp",ppExpND,"");
    ppDexp.setError(ppExpNDerr);
    std::cout<<"prompt+prompt & "<<ppcorr<<" & "<<*ppA.format(2,"EXPP")<<"&"<<*ppB.format(2,"EXPP")<<"&"<<*ppC.format(2,"EXPP")<<"&"<<*ppD.format(2,"EXPP")<<"&"<<*ppDexp.format(2,"EXPP")<<"&"<<ppDiff<<"\\\\"<<std::endl;
*/

  }
  else std::cout << "Unable to open ResultsTable Output File" <<std::endl;

}

void ABCDMethod::SetRooVariables(){

  UInt_t fNReg = 4;
  fRooData.resize(fNReg);
  fRooBkg.resize(fNReg);
  fRooSig.resize(fNReg);
  TString cat_name = "";
  TString name = "";

  for (UInt_t cat = 0; cat < fNReg; cat++){// loop over A,B,C,D
    fRooData[cat].resize(1);
    fRooBkg[cat].resize(fNBkg+1);
    fRooSig[cat].resize(fNSig);
//
//    if (cat==0) cat_name = "_A";
//    if (cat==1) cat_name = "_B";
//    if (cat==2) cat_name = "_D";
//    if (cat==3) cat_name = "_C";
// 
//    name = Form("Data%s",cat_name.Data());
//
//    RooRealVar Data(name,name,fData_Int[cat][0]);
//    Data.setError(fData_IntErr[cat][0]);
//    fRooData[cat][0]=Data; 
//    Data.Delete();
// 
//    for (UInt_t mc = 0; mc < fNBkg+1; mc++){
//      name = Form("%s%s",fBkgNames[mc].Data(),cat_name.Data());
//      RooRealVar Bkg(name,name,fBkg_Int[cat][mc]);
//      Bkg.setError(fBkg_IntErr[cat][mc]);
//      fRooBkg[cat][mc]=Bkg;
//      Bkg.Delete();
//    }
//   
//    for (UInt_t mc = 0; mc < fNSig; mc++){
//      name = Form("%s%s",fSigNames[mc].Data(),cat_name.Data());
//      RooRealVar Sig(name,name,fSig_Int[cat][mc]);
//      Sig.setError(fSig_IntErr[cat][mc]);
//      fRooSig[cat][mc]=Sig;
//      Sig.Delete();
//    }
//    
  }// end loop over categories

  fRData.resize(1);
  fRBkg.resize(fNBkg+1);
  fRSig.resize(fNSig);

  fRData[0] = new RooRealVar("Data","Data",fFullData_Int[0]);
  fRData[0]->setError(fFullData_IntErr[0]);

  for (UInt_t mc = 0; mc < fNBkg+1; mc++){
    if (mc==fNBkg+1) name = "TotBkg";
    else name = fBkgNames[mc];    
    fRBkg[mc] = new RooRealVar(name,name,fFullBkg_Int[mc]);
    fRBkg[mc]->setError(fFullBkg_IntErr[mc]);
  }

  for (UInt_t mc = 0; mc < fNSig; mc++){
    fRSig[mc] = new RooRealVar(fSigNames[mc],fSigNames[mc],fFullSig_Int[mc]);
    fRSig[mc]->setError(fFullSig_IntErr[mc]);
  }


}

void ABCDMethod::WriteDataCard( const TString fSigName){

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
    fOutTxtFile <<  "observation  0 "  << std::endl;
    fOutTxtFile << "------------------------------" << std::endl;
    fOutTxtFile << "bin     1             1            1           1            1               1 "<< std::endl;
    fOutTxtFile << "process DM            Xgg          pp          pf           hgg             Vh " << std::endl;
    fOutTxtFile << "process 0             1            2           3            4               5 " << std::endl;
    fOutTxtFile << Form("rate %.3f          0.57         1.98        0.13         0.043           1.2 ",0.01)<< std::endl; //FIXME sigrate
    fOutTxtFile << "--------------------------------" << std::endl;
    fOutTxtFile << "#signal related" << std::endl;
    fOutTxtFile << "lumi_8TeV     lnN     1.026000      -          -          -       1.02600   1.02600" << std::endl;
    fOutTxtFile << "eff_trig      lnN     1.010000      -          -          -       1.01000   1.01000" << std::endl;
    fOutTxtFile << "id_eff_eb     lnN     1.02000       -          -          -       1.02000   1.02000   " << std::endl;    
    fOutTxtFile << "vtxEff        lnN   0.996/1.008     -          -          -     0.996/1.008  0.996/1.008" << std::endl; 
    fOutTxtFile << "#background related" << std::endl;
    fOutTxtFile << "abcd_estimate lnN       -         1.26000    1.26000    1.26000      -         -  " << std::endl;
    fOutTxtFile << "Xgg_norm      gmN 19 -         0.03         -          -          -         -  " << std::endl;
    fOutTxtFile << "pp_norm       gmN 36617 -           -        0.000054     -          -         -  " << std::endl;
    fOutTxtFile << "pf_norm       gmN 24408 -           -          -       0.0000053     -         -  " << std::endl;


   

  }
  else std::cout << "Unable to open DataCard Output File" << std::endl;

  fOutTxtFile.close();
}

Double_t ABCDMethod::ComputeIntAndErr(TH2D *& h, Double_t & error, const Double_t minX, const Double_t maxX, const Double_t minY, const Double_t maxY, const UInt_t isReg ){

  Double_t integral = 0.;

  if(h == (TH2D*) NULL) std::cout << "NULL TH2D" << std::endl;

  //std::cout << "minx = " << minX << " maxX = " << maxX << " minY = " << minY << " maxY = " << maxY << std::endl; 

  Int_t binXmin;
  Int_t binXmax;
  Int_t binYmin;
  Int_t binYmax;    
  if (isReg == 4 || isReg == 5 || isReg == 6){ // if signal find the exact bins
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

