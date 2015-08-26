
  // Compile:
  // make clean
  // make
  //
  // Call this script:
  // ./main
  //



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

  TString inDir = "./data/50ns/";
  TString outDir = "./diPhoPlots/50ns/";

  bool doFakeData = false;
  bool doTest = false;
  bool doReweightPU = true;
  bool doPlots = false;
  bool doComb = false;

  Double_t lumi = 300.;  
  UInt_t nBins_vtx = 60; 
  
  //for CMSSW_7_0_pre9: run with root
  //gROOT->LoadMacro("Plotter.cpp++g");
  //Plotter * test1 = new Plotter("./data/ALL_nosel/diPhotons","./diPhoPlots/ALL_nosel/","DMHtoGG",30);

  /////////////////////////////////////////////////////
  //
  // Pile up reweighting
  //
  // Inputs to ReweightPU
  // 1st : MC Sample to weight
  // 2nd : Data Sample
  // 3rd : lumi
  // 4th : number of bins for nvtx
  // 5th : input directory of samples
  // 6th : output directory
  //
  /////////////////////////////////////////////////////

  DblVec	puweights_Data;
  DblVec 	puweights_QCD;
  DblVec 	puweights_GJets;
  DblVec	puweights_GGHGG;
  DblVec	puweights_sig1;
  DblVec	puweights_sig10;
  DblVec	puweights_sig100;
  DblVec	puweights_sig1000;	

  // no puweight for data 
  for (UInt_t i=1; i<=nBins_vtx; i++){
    puweights_Data.push_back(1.0);
    puweights_GJets.push_back(1.0);
    puweights_GGHGG.push_back(1.0);
    puweights_sig1000.push_back(1.0);
    puweights_sig100.push_back(1.0);
    puweights_sig10.push_back(1.0);
    puweights_sig1.push_back(1.0);
  }

  if (doReweightPU){ 
    std::cout << "Doing PU Reweighting" << std::endl;
    ReweightPU * reweight = new ReweightPU("QCD","FakeData",lumi, nBins_vtx, inDir, outDir);
    puweights_QCD = reweight->GetPUWeights();
    delete reweight;
  
  }// end doReweightPU
  else{ // if not doReweightPU, set puweights to 1
    std::cout << "No PU Reweighting applied" << std::endl;
    for (UInt_t i=1; i<=nBins_vtx; i++){
      puweights_QCD.push_back(1.0);
      //puweights_GJets.push_back(1.0);
      //puweights_GGHGG.push_back(1.0);
      //puweights_sig1000.push_back(1.0);
      //puweights_sig100.push_back(1.0);
      //puweights_sig10.push_back(1.0);
      //puweights_sig1.push_back(1.0);
    }
  }  

  std::cout << "Finished PU Reweighting" << std::endl;

  /////////////////////////////////////////////////////
  //
  // Make plots for each sample
  //
  // Arguments to Plotter:
  // 1st : location of input data
  // 2nd : output data location
  // 3rd : name of sample 
  // 4th : lumi of data
  //
  /////////////////////////////////////////////////////

  if (doTest){
    std::cout << "Working on test sample" << std::endl;
    Plotter * test = new Plotter(inDir,outDir,"GJets",puweights_GJets,lumi);
    test->DoPlots();
    delete test;
    std::cout << "Finished test sample" << std::endl;
  }
  if (doFakeData){
    std::cout << "Working on FakeData sample" << std::endl;
    Plotter * FakeData = new Plotter(inDir,outDir,"FakeData",puweights_Data,lumi);
    FakeData->DoPlots();
    delete FakeData;
    std::cout << "Finished FakeData sample" << std::endl;
  }
  if (doPlots){
    std::cout << "Working on GJets sample" << std::endl;
    Plotter * GJets = new Plotter(inDir,outDir,"GJets",puweights_GJets,lumi);
    GJets->DoPlots();
    delete GJets;
    std::cout << "Finished GJets sample" << std::endl;

    std::cout << "Working on QCD sample" << std::endl;
    Plotter * QCD = new Plotter(inDir,outDir,"QCD",puweights_QCD,lumi);
    QCD->DoPlots();
    delete QCD;
    std::cout << "Finished QCD sample" << std::endl;

    std::cout << "Working on GluGluH sample" << std::endl;
    Plotter * GGHGG = new Plotter(inDir,outDir,"GluGluHToGG",puweights_GGHGG,lumi);
    GGHGG->DoPlots();
    delete GGHGG;
    std::cout << "Finished GluGluH sample" << std::endl;
  
    std::cout << "Working on DMHgg M1000 sample" << std::endl;
    Plotter * DMH_M1000 = new Plotter(inDir,outDir,"DMHtoGG_M1000",puweights_sig1000,lumi);
    DMH_M1000->DoPlots();
    delete DMH_M1000;
    std::cout << "Finished DMHgg M1000 sample" << std::endl;
  
    std::cout << "Working on DMHgg M100 sample" << std::endl;
    Plotter * DMH_M100 = new Plotter(inDir,outDir,"DMHtoGG_M100",puweights_sig100,lumi);
    DMH_M100->DoPlots();
    delete DMH_M100;
    std::cout << "Finished DMHgg M100 sample" << std::endl;
  
    std::cout << "Working on DMHgg M10 sample" << std::endl;
    Plotter * DMH_M10 = new Plotter(inDir,outDir,"DMHtoGG_M10",puweights_sig10,lumi);
    DMH_M10->DoPlots();
    delete DMH_M10;
    std::cout << "Finished DMHgg M10 sample" << std::endl;
  
    std::cout << "Working on DMHgg M1 sample" << std::endl;
    Plotter * DMH_M1 = new Plotter(inDir,outDir,"DMHtoGG_M1",puweights_sig1,lumi);
    DMH_M1->DoPlots();
    delete DMH_M1;
    std::cout << "Finished DMHgg M1 sample" << std::endl;
  }// end doPlots

  ////////////////////////////////////////////////////
  //
  // Make comb (stack & overlay) plots w/ all samples 
  //
  // Arguments of Combiner
  // 1st : SamplePairVec (Samples) that has Name,VALUE
  // 2rd : lumi
  // 3rd : ColorMap for samples
  // 4th : output directory
  // 5th : bool do N-1 plots 
  //
  ////////////////////////////////////////////////////

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
    Combiner *combAll = new Combiner(Samples,lumi,colorMap,outDir,false);
    combAll->DoComb();
    delete combAll;   
  
    Combiner *combAlln1 = new Combiner(Samples,lumi,colorMap,outDir,true);
    combAlln1->DoComb();
    delete combAlln1;   
  
  }// end doComb
}// end main
