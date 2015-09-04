#include "ABCDMethod.hh"

ABCDMethod::ABCDMethod( SamplePairVec Samples, const Double_t inLumi, const TString outdir){

  lumi = inLumi;
  fOutDir = outdir+"ABCD";
  

  MakeOutDirectory(Form("%s",fOutDir.Data()));
  fOutFile = new TFile(Form("%s/analysis.root",fOutDir.Data()),"RECREATE");
  CheckValidFile(fOutFile, Form("%s/analysis.root",fOutDir.Data())); 

}

ABCDMethod::~ABCDMethod(){
  std::cout << "Finished ABCD & Deleting" << std::endl;
 
  delete fOutFile;

}


void ABCDMethod::DoAnalysis(){

}
