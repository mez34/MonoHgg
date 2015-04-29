#include "HighMassGGFitter.cc"
//

void ProduceWorkspaces(int mass, Double_t width, std::string model){
  
  //gSystem->SetIncludePath("-I/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.00-cms5/include");
  // gSystem->Load("libRooFit");
  // gROOT->ProcessLine(".x RooCBCrujffPdf.cxx+");
  //gROOT->ProcessLine(".x RooCBCBPdf.cxx+");
  cout << "mass = " << mass << endl; 
  // runfits(mass, true);
  runfits(mass, true, width, model);
}
