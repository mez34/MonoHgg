
  // Compile:
  // make clean
  // make
  //
  // Call this script:
  // ./main
  //
  // Arguments to Plotter:
  // 1st : location of input data
  // 2nd : output data location
  // 3rd : name of sample 
  // 4th : lumi of data


#include "Plotter.hh"
#include "Combiner.hh"
#include "ReweightPU.hh"
#include "Style.hh"


#include "TROOT.h"
#include "TStyle.h"

#include <iostream>

typedef std::pair<TString,Double_t>  SampleYieldPair;
typedef std::vector<SampleYieldPair> SampleYieldPairVec;

static bool sortByYield(const SampleYieldPair& mcpair1, const SampleYieldPair& mcpair2) {
  return mcpair1.second<=mcpair2.second;
}

int main(){
  setTDRStyle();

  TString inDir = "data/50ns/";
  TString outDir = "diPhoPlots/50ns/";

  bool doFakeData = false;
  bool doTest = false;
  bool doReweightPU = true;
  bool doPlots = true;
  bool doComb = true;
  
  
  //for CMSSW_7_0_pre9: run with root
  //gROOT->LoadMacro("Plotter.cpp++g");
  //Plotter * test1 = new Plotter("./data/ALL_nosel/diPhotons","./diPhoPlots/ALL_nosel/","DMHtoGG",30);

  if (doTest){
    std::cout << "Working on test sample" << std::endl;
    Plotter * test = new Plotter("./data/50ns/","./diPhoPlots/50ns/","GJets",30);
    test->DoPlots();
    delete test;
    std::cout << "Finished test sample" << std::endl;
  }
  if (doFakeData){
    std::cout << "Working on FakeData sample" << std::endl;
    Plotter * FakeData = new Plotter("./data/50ns/","./diPhoPlots/50ns/","FakeData",30);
    FakeData->DoPlots();
    delete FakeData;
    std::cout << "Finished FakeData sample" << std::endl;
  }

 //--------------------------------------------------
 //
 // Make plots for each sample
 //
 //--------------------------------------------------

 if (doPlots){
    std::cout << "Working on GJets sample" << std::endl;
    Plotter * GJets = new Plotter("./data/50ns/","./diPhoPlots/50ns/","GJets",30);
    GJets->DoPlots();
    delete GJets;
    std::cout << "Finished GJets sample" << std::endl;

    std::cout << "Working on QCD sample" << std::endl;
    Plotter * QCD = new Plotter("./data/50ns/","./diPhoPlots/50ns/","QCD",30);
    QCD->DoPlots();
    delete QCD;
    std::cout << "Finished QCD sample" << std::endl;

    std::cout << "Working on GluGluH sample" << std::endl;
    Plotter * GGHGG = new Plotter("./data/50ns/","./diPhoPlots/50ns/","GluGluHToGG",30);
    GGHGG->DoPlots();
    delete GGHGG;
    std::cout << "Finished GluGluH sample" << std::endl;
  
    std::cout << "Working on DMHgg M1000 sample" << std::endl;
    Plotter * DMH_M1000 = new Plotter("./data/50ns/","./diPhoPlots/50ns/","DMHtoGG_M1000",30);
    DMH_M1000->DoPlots();
    delete DMH_M1000;
    std::cout << "Finished DMHgg M1000 sample" << std::endl;
  
    std::cout << "Working on DMHgg M100 sample" << std::endl;
    Plotter * DMH_M100 = new Plotter("./data/50ns/","./diPhoPlots/50ns/","DMHtoGG_M100",30);
    DMH_M100->DoPlots();
    delete DMH_M100;
    std::cout << "Finished DMHgg M100 sample" << std::endl;
  
    std::cout << "Working on DMHgg M10 sample" << std::endl;
    Plotter * DMH_M10 = new Plotter("./data/50ns/","./diPhoPlots/50ns/","DMHtoGG_M10",30);
    DMH_M10->DoPlots();
    delete DMH_M10;
    std::cout << "Finished DMHgg M10 sample" << std::endl;
  
    std::cout << "Working on DMHgg M1 sample" << std::endl;
    Plotter * DMH_M1 = new Plotter("./data/50ns/","./diPhoPlots/50ns/","DMHtoGG_M1",30);
    DMH_M1->DoPlots();
    delete DMH_M1;
    std::cout << "Finished DMHgg M1 sample" << std::endl;
  }// end doPlots

 //--------------------------------------------------
 //
 // Make comb (stack & overlay) plots w/ all samples 
 //
 //--------------------------------------------------

  //DblVec puweights;
  if (doReweightPU){ 
   
    std::cout << "Doing PU Reweighting" << std::endl;
  /*  ReweightPU * reweight = new ReweightPU(PURWSamples, PURWselection, PURWnjetsselection, lumi, nBins_vtx, outdir)
    puweights = reweight->GetPUWeights();
    delete reweight;
   */
  }// end doReweightPU


  if (doComb){
    ColorMap colorMap;
    colorMap["QCD"] 		= kYellow;
    colorMap["GJets"] 		= kGreen;
    colorMap["GluGluHToGG"]	= kCyan;
    colorMap["DMHtoGG_M1000"]	= kMagenta;
    colorMap["DMHtoGG_M100"]	= kMagenta+1;
    colorMap["DMHtoGG_M10"]	= kRed+1;
    colorMap["DMHtoGG_M1"]	= kRed;

    SamplePairVec Samples; // vector to also be used for stack plots
    Samples.push_back(SamplePair("QCD",1)); 
    Samples.push_back(SamplePair("GJets",1)); 
    Samples.push_back(SamplePair("GluGluHToGG",1)); 
    Samples.push_back(SamplePair("DMHtoGG_M1000",0)); 
    Samples.push_back(SamplePair("DMHtoGG_M100",0)); 
    Samples.push_back(SamplePair("DMHtoGG_M10",0)); 
    Samples.push_back(SamplePair("DMHtoGG_M1",0)); 
    if (doFakeData) Samples.push_back(SamplePair("FakeData",5));
  
    UInt_t nbkg = 0;
    UInt_t nsig = 0;
    UInt_t ndata = 0;
  
    for (SamplePairVecIter iter=Samples.begin(); iter != Samples.end(); ++iter){
      std::cout << "Analyzing Sample: "<< (*iter).first.Data() << std::endl;
      if ((*iter).second == 1) {nbkg++;}
      else if ((*iter).second == 0) {nsig++;}
      else {ndata++;} 
    }
    UInt_t nsamples = nbkg + nsig + ndata;
   
    SamplePairVec BkgSamples;
    SamplePairVec SigSamples;
    SamplePairVec DataSamples;
    for (UInt_t isample = 0; isample < nsamples; isample++){
      if (Samples[isample].second == 0) SigSamples.push_back(Samples[isample]);
      else if (Samples[isample].second == 1) BkgSamples.push_back(Samples[isample]);
      else  DataSamples.push_back(Samples[isample]);
    }
  
    // to sort MC by smallest to largest for nice stacked plots
    SampleYieldPairVec tmp_mcyields;
    for (UInt_t mc = 0; mc < nbkg; mc++) {
        // open mc file first
        TString mcfilename = Form("diPhoPlots/50ns/%s/plots_%s.root",BkgSamples[mc].first.Data(),BkgSamples[mc].first.Data());
        TFile * tmp_mcfile = TFile::Open(mcfilename.Data());
        // open nvtx plot
        TH1D * tmpnvtx = (TH1D*)tmp_mcfile->Get("nvtx_n-1");
        // get yield and push back with corresponding sample name
        tmp_mcyields.push_back(SampleYieldPair(BkgSamples[mc].first,tmpnvtx->Integral()));
  
        delete tmpnvtx;
        delete tmp_mcfile;
     }
  
     std::sort(tmp_mcyields.begin(),tmp_mcyields.end(),sortByYield);
     std::cout << "Finished sorting MC, now put samples in right order to be processed" << std::endl;
     BkgSamples.clear();
      for (UInt_t mc = 0; mc < nbkg; mc++) { // init mc double hists
        BkgSamples.push_back(SamplePair(tmp_mcyields[mc].first,1));
      }
  
    Samples.clear();
    for (UInt_t data = 0; data < ndata; data++ ) {
      Samples.push_back(DataSamples[data]);
    }
    for (UInt_t mc = 0; mc < nbkg; mc++ ) {
      Samples.push_back(BkgSamples[mc]);
    }
    for (UInt_t mc = 0; mc < nsig; mc++) {
      Samples.push_back(SigSamples[mc]);
    }
  
    // make overlayed and stack plots
    // Combiner( Samples, lumi, colorMap , outDir, doNmin1plots )
    Combiner *combAll = new Combiner(Samples,300,colorMap,outDir,false);
    combAll->DoComb();
    delete combAll;   
  
    Combiner *combAlln1 = new Combiner(Samples,300,colorMap,outDir,true);
    combAlln1->DoComb();
    delete combAlln1;   
  
  }// end doComb
}// end main
