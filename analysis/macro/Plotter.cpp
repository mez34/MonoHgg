#include "Plotter.hh"
#include "Style.hh"
#include "../../../DataFormats/Math/interface/deltaPhi.h"
//#include "mkPlotsLivia/CMS_lumi.C"

Plotter::Plotter( TString inName, TString outName, TString inSpecies, const Double_t lumi){

  // Get ROOT file
  name = inName;
  species = inSpecies;
  inFile = TFile::Open(Form("%s%s.root",name.Data(),species.Data()));
  CheckValidFile(inFile,Form("%s%s.root",name.Data(),species.Data()));  
  // Open Tree from inFile
  tpho = (TTree*)inFile->Get("DiPhotonTree"); 
  CheckValidTree(tpho,"DiPhotonTree",Form("%s%s.root",name.Data(),species.Data()));

  fLumi = lumi;

  // Make output directory
  fName = outName;
  TString FullPath = fName.Data();
  FullPath+=species.Data();
  FullPath+="/";
  MakeOutDirectory(FullPath.Data());

  // Make output ROOT file
  outFile = new TFile(Form("%s/%s/plots_%s.root",fName.Data(),species.Data(),species.Data()),"RECREATE");
  CheckValidFile(outFile,Form("%s/%s/plots_%s.root",fName.Data(),species.Data(),species.Data()));

  // Make TCanvas
  fTH1Canv = new TCanvas();
  fTH2Canv = new TCanvas();

  // set all the branch addresses appropriately
  Plotter::SetBranchAddresses();

  // Initialize all of the variables
  /*Plotter::InitTreeVar();
  NVARIABLES = varname.size();
  Plotter::InitTreeEffVar();
  N2DVARIABLES = effvar.size();
  Plotter::InitPhotonIDSel();
  NSEL = selvar.size();
  */

}// end Plotter::Plotter


Plotter::~Plotter(){
  std::cout << "Finished & Deleting" <<std::endl;
  std::cout << "Deleting Canvases" <<std::endl;
  delete fTH1Canv;
  delete fTH2Canv;
  std::cout << "Deleting inTree" <<std::endl;
  delete tpho;
  std::cout << "Deleting inFile" <<std::endl;
  delete inFile;
  // Write and Close output ROOT file
  //Plotter::DeleteBranches();
  std::cout << "Delete histos" <<std::endl;
  Plotter::DeleteHists();
  std::cout << "Deleting outFile" <<std::endl;
  delete outFile;
  std::cout << "Finished Deleting" <<std::endl;
}// end Plotter::~Plotter


void Plotter::DoPlots(){
  Plotter::SetUpPlots();
 
  nphotons = tpho->GetEntries(); 
  for (UInt_t entry = 0; entry < nphotons; entry++){
    tpho->GetEntry(entry);

    // calculate the weight
    Double_t Weight = weight;

    //Fill histograms
    fTH1DMap["mgg"]->Fill(mgg,Weight);
    fTH1DMap["nvtx"]->Fill(nvtx,Weight);
    fTH1DMap["ptgg"]->Fill(ptgg,Weight);
    fTH1DMap["pt1"]->Fill(pt1,Weight);
    fTH1DMap["pt2"]->Fill(pt2,Weight);
    fTH1DMap["t1pfmet"]->Fill(t1pfmet,Weight);
    fTH1DMap["t1pfmetphi"]->Fill(t1pfmetphi,Weight);
    fTH1DMap["phi1"]->Fill(phi1,Weight);
    fTH1DMap["phi2"]->Fill(phi2,Weight);
    fTH1DMap["eta1"]->Fill(eta1,Weight);
    fTH1DMap["eta2"]->Fill(eta2,Weight);
    fTH1DMap["chiso1"]->Fill(chiso1,Weight);
    fTH1DMap["chiso2"]->Fill(chiso2,Weight);
    fTH1DMap["neuiso1"]->Fill(neuiso1,Weight);
    fTH1DMap["neuiso2"]->Fill(neuiso2,Weight);
    fTH1DMap["phoiso1"]->Fill(phoiso1,Weight);
    fTH1DMap["phoiso2"]->Fill(phoiso2,Weight);
    fTH1DMap["sieie1"]->Fill(sieie1,Weight);
    fTH1DMap["sieie2"]->Fill(sieie2,Weight);
    fTH1DMap["hoe1"]->Fill(hoe1,Weight);
    fTH1DMap["hoe2"]->Fill(hoe2,Weight);
    fTH1DMap["r91"]->Fill(r91,Weight);
    fTH1DMap["r92"]->Fill(r92,Weight);

    Bool_t passCH1 = false;
    Bool_t passCH2 = false;
    Bool_t passNH1 = false;
    Bool_t passNH2 = false;
    Bool_t passPH1 = false;
    Bool_t passPH2 = false;
    Bool_t passS1 = true;
    Bool_t passS2 = true;
    Bool_t passHE1 = false;
    Bool_t passHE2 = false;
    Bool_t passAll1 = false;
    Bool_t passAll2 = false;
    Bool_t passBoth = false;

    if (passCHiso1==0) passCH1 = true; 
    if (passCHiso2==0) passCH2 = true; 
    if (passNHiso1==0) passNH1 = true;
    if (passNHiso2==0) passNH2 = true;
    if (passPHiso1==0) passPH1 = true;
    if (passPHiso2==0) passPH2 = true;
    if (passSieie1==0) passS1 = true;
    if (passSieie2==0) passS2 = true;
    if (passHoe1==0)   passHE1 = true; 
    if (passHoe2==0)   passHE2 = true; 

    if (passCH1 && passNH1 && passPH1 && passS1 && passHE1) passAll1 = true;
    if (passCH2 && passNH2 && passPH2 && passS2 && passHE2) passAll2 = true;
    if (passAll1 && passAll2) passBoth = true;

    //std::cout << passCH1 <<" "<< passNH1 <<" "<< passPH1 <<" "<< passHE1 <<" "<< passS1 << std::endl; 
    //std::cout << passCH2 <<" "<< passNH2 <<" "<< passPH2 <<" "<< passHE2 <<" "<< passS2 << std::endl; 
    //std::cout << passAll1 <<" "<< passAll2 <<" "<< passBoth << std::endl;

    //fill n-1 plots for the photon ID selection variables
    if (passCH1 && passNH1 && passPH1 && passS1)  fTH1DMap["hoe1_n-1"]->Fill(hoe1,Weight); 
    if (passCH1 && passNH1 && passPH1 && passHE1) fTH1DMap["sieie1_n-1"]->Fill(sieie1,Weight);
    if (passCH1 && passNH1 && passHE1 && passS1)  fTH1DMap["phoiso1_n-1"]->Fill(phoiso1,Weight);
    if (passCH1 && passPH1 && passHE1 && passS1)  fTH1DMap["neuiso1_n-1"]->Fill(neuiso1,Weight);
    if (passPH1 && passNH1 && passHE1 && passS1)  fTH1DMap["chiso1_n-1"]->Fill(chiso1,Weight);

    if (passCH2 && passNH2 && passPH2 && passS2)  fTH1DMap["hoe2_n-1"]->Fill(hoe2,Weight); 
    if (passCH2 && passNH2 && passPH2 && passHE2) fTH1DMap["sieie2_n-1"]->Fill(sieie2,Weight);
    if (passCH2 && passNH2 && passHE2 && passS2)  fTH1DMap["phoiso2_n-1"]->Fill(phoiso2,Weight);
    if (passCH2 && passPH2 && passHE2 && passS2)  fTH1DMap["neuiso2_n-1"]->Fill(neuiso2,Weight);
    if (passPH2 && passNH2 && passHE2 && passS2)  fTH1DMap["chiso2_n-1"]->Fill(chiso2,Weight);

    if (passAll1){// fill pho1 plots if these photons pass phoID
      fTH1DMap["pt1_n-1"]->Fill(pt1,Weight);
      fTH1DMap["r91_n-1"]->Fill(r91,Weight);
      fTH1DMap["phi1_n-1"]->Fill(phi1,Weight);
      fTH1DMap["eta1_n-1"]->Fill(eta1,Weight);
    }
    if (passAll2){// fill pho2 plots if these photons pass phoID
      fTH1DMap["pt2_n-1"]->Fill(pt2,Weight);
      fTH1DMap["r92_n-1"]->Fill(r92,Weight);
      fTH1DMap["phi2_n-1"]->Fill(phi2,Weight);
      fTH1DMap["eta2_n-1"]->Fill(eta2,Weight);
    } 
    if (passBoth){
      fTH1DMap["nvtx_n-1"]->Fill(nvtx,Weight);
      fTH1DMap["mgg_n-1"]->Fill(mgg,Weight);  
      fTH1DMap["ptgg_n-1"]->Fill(ptgg,Weight);  
      fTH1DMap["t1pfmet_n-1"]->Fill(t1pfmet,Weight);  
      fTH1DMap["t1pfmetphi_n-1"]->Fill(t1pfmetphi,Weight);  
    }

  }// end loop over entries in tree

  Plotter::SavePlots();


/* //OLD IMPLEMENTATION:
  Plotter::getTree();
  std::cout << "Here1" << std::endl;
  Plotter::make1DHistos();
  std::cout << "Here2" << std::endl;
  Plotter::make2DHistos();
  outFile->Write(); */


}// end Plotter::DoPlots


void Plotter::SetUpPlots(){
  // fill all plots from tree
  fTH1DMap["nvtx"]	= Plotter::MakeTH1DPlot("nvtx","",60,0.,60.,"nvtx","");
  fTH1DMap["mgg"]	= Plotter::MakeTH1DPlot("mgg","",60,50.,300.,"m_{#gamma#gamma} (GeV)","");  
  fTH1DMap["ptgg"]	= Plotter::MakeTH1DPlot("ptgg","",100,0.,1000.,"p_{T,#gamma#gamma} (GeV)","");
  fTH1DMap["t1pfmet"]	= Plotter::MakeTH1DPlot("t1pfmet","",100,0.,1000,"t1PF MET (GeV)","");
  fTH1DMap["t1pfmetphi"]= Plotter::MakeTH1DPlot("t1pfmetphi","",80,-4.,4.,"MET #phi","");
  fTH1DMap["phi1"]	= Plotter::MakeTH1DPlot("phi1","",80,-4.,4.,"#phi(#gamma1)","");
  fTH1DMap["phi2"]	= Plotter::MakeTH1DPlot("phi2","",80,-4.,4.,"#phi(#gamma2)","");
  fTH1DMap["eta1"]	= Plotter::MakeTH1DPlot("eta1","",100,-5.,5.,"#eta(#gamma1)","");
  fTH1DMap["eta2"]	= Plotter::MakeTH1DPlot("eta2","",100,-5.,5.,"#eta(#gamma2)","");
  fTH1DMap["pt1"]	= Plotter::MakeTH1DPlot("pt1","",50,0.,500.,"p_{T,#gamma1} (GeV)","");
  fTH1DMap["pt2"]	= Plotter::MakeTH1DPlot("pt2","",50,0.,500.,"p_{T,#gamma2} (GeV)","");
  fTH1DMap["chiso1"]	= Plotter::MakeTH1DPlot("chiso1","",150,-5.,15.,"CHiso(#gamma1)","");
  fTH1DMap["chiso2"]	= Plotter::MakeTH1DPlot("chiso2","",150,-5.,15.,"CHiso(#gamma2)","");
  fTH1DMap["neuiso1"]	= Plotter::MakeTH1DPlot("neuiso1","",150,-5.,15.,"NHiso(#gamma1)","");
  fTH1DMap["neuiso2"]	= Plotter::MakeTH1DPlot("neuiso2","",150,-5.,15.,"NHiso(#gamma2)","");
  fTH1DMap["phoiso1"]	= Plotter::MakeTH1DPlot("phoiso1","",150,-5.,15.,"PHiso(#gamma1)",""); 
  fTH1DMap["phoiso2"]	= Plotter::MakeTH1DPlot("phoiso2","",150,-5.,15.,"PHiso(#gamma2)",""); 
  fTH1DMap["sieie1"]	= Plotter::MakeTH1DPlot("sieie1","",300,0.,0.03,"#sigma_{i#eta i#eta}(#gamma1)",""); 
  fTH1DMap["sieie2"]	= Plotter::MakeTH1DPlot("sieie2","",300,0.,0.03,"#sigma_{i#eta i#eta}(#gamma2)",""); 
  fTH1DMap["hoe1"]	= Plotter::MakeTH1DPlot("hoe1","",250,0.,0.025,"H/E(#gamma1)","");
  fTH1DMap["hoe2"]	= Plotter::MakeTH1DPlot("hoe2","",250,0.,0.025,"H/E(#gamma2)","");
  fTH1DMap["r91"]	= Plotter::MakeTH1DPlot("r91","",100,0.,1.1,"R9(#gamma1)","");
  fTH1DMap["r92"]	= Plotter::MakeTH1DPlot("r92","",100,0.,1.1,"R9(#gamma2)","");

  //n minus 1 plots
  fTH1DMap["nvtx_n-1"]		= Plotter::MakeTH1DPlot("nvtx_n-1","",60,0.,60.,"nvtx","");
  fTH1DMap["mgg_n-1"]		= Plotter::MakeTH1DPlot("mgg_n-1","",60,50.,300.,"m_{#gamma#gamma} (GeV)","");  
  fTH1DMap["ptgg_n-1"]		= Plotter::MakeTH1DPlot("ptgg_n-1","",100,0.,1000.,"p_{T,#gamma#gamma} (GeV)","");
  fTH1DMap["t1pfmet_n-1"]	= Plotter::MakeTH1DPlot("t1pfmet_n-1","",100,0.,1000,"t1PF MET (GeV)","");
  fTH1DMap["t1pfmetphi_n-1"]	= Plotter::MakeTH1DPlot("t1pfmetphi_n-1","",80,-4.,4.,"MET #phi","");
  fTH1DMap["phi1_n-1"]		= Plotter::MakeTH1DPlot("phi1_n-1","",80,-4.,4.,"#phi(#gamma1)","");
  fTH1DMap["phi2_n-1"]		= Plotter::MakeTH1DPlot("phi2_n-1","",80,-4.,4.,"#phi(#gamma2)","");
  fTH1DMap["eta1_n-1"]		= Plotter::MakeTH1DPlot("eta1_n-1","",100,-5.,5.,"#eta(#gamma1)","");
  fTH1DMap["eta2_n-1"]		= Plotter::MakeTH1DPlot("eta2_n-1","",100,-5.,5.,"#eta(#gamma2)","");
  fTH1DMap["pt1_n-1"]		= Plotter::MakeTH1DPlot("pt1_n-1","",50,0.,500.,"p_{T,#gamma1} (GeV)","");
  fTH1DMap["pt2_n-1"]		= Plotter::MakeTH1DPlot("pt2_n-1","",50,0.,500.,"p_{T,#gamma2} (GeV)","");
  fTH1DMap["chiso1_n-1"]	= Plotter::MakeTH1DPlot("chiso1_n-1","",150,-5.,15.,"CHiso(#gamma1)","");
  fTH1DMap["chiso2_n-1"]	= Plotter::MakeTH1DPlot("chiso2_n-1","",150,-5.,15.,"CHiso(#gamma2)","");
  fTH1DMap["neuiso1_n-1"]	= Plotter::MakeTH1DPlot("neuiso1_n-1","",150,-5.,15.,"NHiso(#gamma1)","");
  fTH1DMap["neuiso2_n-1"]	= Plotter::MakeTH1DPlot("neuiso2_n-1","",150,-5.,15.,"NHiso(#gamma2)","");
  fTH1DMap["phoiso1_n-1"]	= Plotter::MakeTH1DPlot("phoiso1_n-1","",150,-5.,15.,"PHiso(#gamma1)",""); 
  fTH1DMap["phoiso2_n-1"]	= Plotter::MakeTH1DPlot("phoiso2_n-1","",150,-5.,15.,"PHiso(#gamma2)",""); 
  fTH1DMap["sieie1_n-1"]	= Plotter::MakeTH1DPlot("sieie1_n-1","",300,0.,0.03,"#sigma_{i#eta i#eta}(#gamma1)",""); 
  fTH1DMap["sieie2_n-1"]	= Plotter::MakeTH1DPlot("sieie2_n-1","",300,0.,0.03,"#sigma_{i#eta i#eta}(#gamma2)",""); 
  fTH1DMap["hoe1_n-1"]		= Plotter::MakeTH1DPlot("hoe1_n-1","",250,0.,0.025,"H/E(#gamma1)","");
  fTH1DMap["hoe2_n-1"]		= Plotter::MakeTH1DPlot("hoe2_n-1","",250,0.,0.025,"H/E(#gamma2)","");
  fTH1DMap["r91_n-1"]		= Plotter::MakeTH1DPlot("r91_n-1","",100,0.,1.1,"R9(#gamma1)","");
  fTH1DMap["r92_n-1"]		= Plotter::MakeTH1DPlot("r92_n-1","",100,0.,1.1,"R9(#gamma2)","");

}// end Plotter::SetUpPlots

TH1D * Plotter::MakeTH1DPlot(const TString hname, const TString htitle, const Int_t nbins, const Double_t xlow, const Double_t xhigh, const TString xtitle, const TString ytitle){
  TH1D * hist = new TH1D(hname.Data(),htitle.Data(),nbins,xlow,xhigh);
  hist->GetXaxis()->SetTitle(xtitle.Data());
  hist->GetYaxis()->SetTitle(ytitle.Data());
  return hist;
}// end Plotter::MakeTH1DPlot

void Plotter::DoAnalysis(){

}// end Plotter::DoAnalysis

void Plotter::SavePlots(){
  outFile->cd();

  TCanvas * canv = new TCanvas();

  for (TH1DMapIter mapiter = fTH1DMap.begin(); mapiter != fTH1DMap.end(); mapiter++){

    if ((*mapiter).second == (TH1D*) NULL)	{std::cout << "TH1D Null" << std::endl;}
    if (outFile == (TFile*) NULL)		{std::cout << "OutFile Null" << std::endl;}
    if (canv == (TCanvas*) NULL)		{std::cout << "Canvas Null" << std::endl;}

    (*mapiter).second->Write(); // save histos to root file 
    canv->cd();
    (*mapiter).second->Draw("HIST");

    canv->SetLogy(0);
    canv->SaveAs(Form("%s%s/%s.png",fName.Data(),species.Data(),(*mapiter).first.Data()));

    canv->SetLogy(1);
    canv->SaveAs(Form("%s%s/%s_log.png",fName.Data(),species.Data(),(*mapiter).first.Data())); 
  }// end of loop over mapiter
  delete canv;

}// end Plotter::SavePlots

void Plotter::DeleteHists(){
  for (TH1DMapIter mapiter = fTH1DMap.begin(); mapiter != fTH1DMap.end(); mapiter++){
    delete ((*mapiter).second);
  }
  fTH1DMap.clear();
}// end Plotter::DeleteHists

void Plotter::SetBranchAddresses(){
  tpho->SetBranchAddress("weight", &weight,  &b_weight);
  tpho->SetBranchAddress("nvtx",   &nvtx,    &b_nvtx);
  tpho->SetBranchAddress("mgg",    &mgg,     &b_mgg);
  tpho->SetBranchAddress("ptgg",   &ptgg,    &b_ptgg);
  tpho->SetBranchAddress("t1pfmet", &t1pfmet, &b_t1pfmet);   
  tpho->SetBranchAddress("t1pfmetPhi", &t1pfmetphi, &b_t1pfmetPhi);   
  tpho->SetBranchAddress("pt1", &pt1, &b_pt1);   
  tpho->SetBranchAddress("pt2", &pt2, &b_pt2);   
  tpho->SetBranchAddress("chiso1", &chiso1, &b_chiso1);   
  tpho->SetBranchAddress("chiso2", &chiso2, &b_chiso2);   
  tpho->SetBranchAddress("neuiso1", &neuiso1, &b_neuiso1);   
  tpho->SetBranchAddress("neuiso2", &neuiso2, &b_neuiso2);   
  tpho->SetBranchAddress("phoiso1", &phoiso1, &b_phoiso1);   
  tpho->SetBranchAddress("phoiso2", &phoiso2, &b_phoiso2);   
  tpho->SetBranchAddress("sieie1", &sieie1, &b_sieie1);   
  tpho->SetBranchAddress("sieie2", &sieie2, &b_sieie2);   
  tpho->SetBranchAddress("hoe1", &hoe1, &b_hoe1);   
  tpho->SetBranchAddress("hoe2", &hoe2, &b_hoe2);   
  tpho->SetBranchAddress("r91", &r91, &b_r91);   
  tpho->SetBranchAddress("r92", &r92, &b_r92);   
  tpho->SetBranchAddress("phi1", &phi1, &b_phi1);   
  tpho->SetBranchAddress("phi2", &phi2, &b_phi2);   
  tpho->SetBranchAddress("eta1", &eta1, &b_eta1);   
  tpho->SetBranchAddress("eta2", &eta2, &b_eta2);   
  tpho->SetBranchAddress("passCHiso1", &passCHiso1, &b_passCHiso1);   
  tpho->SetBranchAddress("passCHiso2", &passCHiso2, &b_passCHiso2);   
  tpho->SetBranchAddress("passNHiso1", &passNHiso1, &b_passNHiso1);   
  tpho->SetBranchAddress("passNHiso2", &passNHiso2, &b_passNHiso2);   
  tpho->SetBranchAddress("passPHiso1", &passPHiso1, &b_passNHiso1);   
  tpho->SetBranchAddress("passPHiso2", &passPHiso2, &b_passNHiso2);   
  tpho->SetBranchAddress("passSieie1", &passSieie1, &b_passSieie1);
  tpho->SetBranchAddress("passSieie2", &passSieie2, &b_passSieie2);
  tpho->SetBranchAddress("passHoe1", &passHoe1, &b_passHoe1);
  tpho->SetBranchAddress("passHoe2", &passHoe2, &b_passHoe2);

  //tpho->SetBranchAddress("", &, &b_);
  
}// end Plotter::SetBranchAddresses


void Plotter::DeleteBranches(){
  delete b_weight;
  delete b_nvtx;
  delete b_mgg;
  delete b_ptgg;
  delete b_pt1;
  delete b_pt2;
}// end Plotter::DeleteBranches














// DELETE AFTER THIS FIXME

void Plotter::getTree(){

  // Load variables from Tree
  variable[NVARIABLES];    //= {-1000}; // float for most variables 
  intvariable[NVARIABLES]; //= {-1000}; // int for other variables (eleveto, sel, nvtx)

  selvarPair[NSEL]; 
  for(int z=0; z<NVARIABLES; ++z){
    if(z < NSEL){
      tpho->SetBranchAddress(selvar[z].Data(),&selvarPair[z]);
    }

    if(z==8 || z==16 || z>=25 ){ //eleveto, sel, nvtx
       tpho->SetBranchAddress(varname[z].Data(),&intvariable[z]);
     }
     else{
       tpho->SetBranchAddress(varname[z].Data(),&variable[z]);
     }
  }
  nphotons = (int)tpho->GetEntries();
}// end Plotter::getTree


void Plotter::make1DHistos(){
  TH1F *hVar[NVARIABLES];
  TH1F *hVarNmin1[NVARIABLES];
  for (int z=0; z<NVARIABLES; ++z){
    hVar[z] = new TH1F(Form("%s_%s",varname[z].Data(),species.Data()),Form("%s_%s",varname[z].Data(),species.Data()),nbins[z],range[z][0],range[z][1]);
    hVar[z]->GetXaxis()->SetTitle(xaxisLabel[z]);
    hVarNmin1[z] = new TH1F(Form("%s_%s_n-1",varname[z].Data(),species.Data()),Form("%s_%s_n-1",varname[z].Data(),species.Data()),nbins[z],range[z][0],range[z][1]);
    hVarNmin1[z]->GetXaxis()->SetTitle(xaxisLabel[z]);
  }// one histogram for each variable in Tree 

  // additional histograms
  TH1F *hPhi[3]; // phi of the higgs, MET and delta phi H,MET
  hPhi[0] = new TH1F(Form("phi_H_%s",species.Data()),Form("phi_H_%s",species.Data()),80,-4,4);
  hPhi[1] = new TH1F(Form("phi_MET_%s",species.Data()),Form("phi_MET_%s",species.Data()),80,-4,4);
  hPhi[2] = new TH1F(Form("phi_HMET_%s",species.Data()),Form("phi_HMET_%s",species.Data()),80,-4,4);

/*  TH1F *hEff[3]; // eff of pho ID vs PU, pt, eta
  hEff[0] = new TH1F(Form("Eff_PHID_PU_%s",species.Data()),Form("Eff_PHID_PU_%s",species.Data()),60,0,60);
  hEff[1] = new TH1F(Form("Eff_PHID_pt_%s",species.Data()),Form("Eff_PHID_pt_%s",species.Data()),60,0,600);
  hEff[2] = new TH1F(Form("Eff_PHID_eta_%s",species.Data()),Form("Eff_PHID_eta_%s",species.Data()),60,-3,3);
*/
  //plot Mgg and MET with some photonIDsel and cut on other var
  TH1F *hMET = new TH1F(Form("t1pfMet_selMgg_%s",species.Data()),Form("t1pfMet_selMgg_%s",species.Data()),100,0,1000);
  TH1F *hMgg = new TH1F(Form("mgg_selt1pfMet_%s",species.Data()),Form("mgg_selt1pfMet_%s",species.Data()),50,50,300);

  Float_t phiH[nphotons];
  //phiH[nphotons] = {-1000}; 
  Float_t phiHMET[nphotons];
  //phiHMET[nphotons] = {-1000};



  int nphotonsPass=0;
  //int Eff[6][60]={0};

  for (int i=0; i<nphotons; ++i){   
    tpho->GetEntry(i);

    for (UInt_t z=0; z<NVARIABLES; z++){
      if (z==8 || z==16 || z>=25) hVar[z]->Fill(intvariable[z],variable[18]);
      else hVar[z]->Fill(variable[z],variable[18]);
    }

    Bool_t passCHiso1 = false;
    Bool_t passCHiso2 = false;
    Bool_t passNHiso1 = false;
    Bool_t passNHiso2 = false;
    Bool_t passPHiso1 = false;
    Bool_t passPHiso2 = false;
    Bool_t passSieie1 = true;
    Bool_t passSieie2 = true;
    Bool_t passHoe1 = false;
    Bool_t passHoe2 = false;
    Bool_t passAll1 = false;
    Bool_t passAll2 = false;
    Bool_t passBoth = false;

    if (selvarPair[0]==0) passCHiso1 = true; 
    if (selvarPair[1]==0) passCHiso2 = true; 
    if (selvarPair[2]==0) passNHiso1 = true;
    if (selvarPair[3]==0) passNHiso2 = true;
    if (selvarPair[4]==0) passPHiso1 = true;
    if (selvarPair[5]==0) passPHiso2 = true;
    if (selvarPair[6]==0) passSieie1 = true;
    if (selvarPair[7]==0) passSieie2 = true;
    if (selvarPair[8]==0) passHoe1 = true; 
    if (selvarPair[9]==0) passHoe2 = true; 

    if (passCHiso1 && passNHiso1 && passPHiso1 && passSieie1 && passHoe1) passAll1 = true;
    if (passCHiso2 && passNHiso2 && passPHiso2 && passSieie2 && passHoe2) passAll2 = true;
    if (passAll1 && passAll2) passBoth = true;

    if (passCHiso1 && passNHiso1 && passPHiso1 && passHoe1)   hVarNmin1[3]->Fill(variable[3],variable[18]);
    if (passCHiso1 && passNHiso1 && passPHiso1 && passSieie1) hVarNmin1[4]->Fill(variable[4],variable[18]);
    if (passNHiso1 && passPHiso1 && passSieie1 && passHoe1)   hVarNmin1[5]->Fill(variable[5],variable[18]);
    if (passCHiso1 && passNHiso1 && passSieie1 && passHoe1)   hVarNmin1[6]->Fill(variable[6],variable[18]);
    if (passCHiso1 && passPHiso1 && passSieie1 && passHoe1)   hVarNmin1[7]->Fill(variable[7],variable[18]);
    if (passCHiso2 && passNHiso2 && passPHiso2 && passHoe2)   hVarNmin1[11]->Fill(variable[11],variable[18]);
    if (passCHiso2 && passNHiso2 && passPHiso2 && passSieie2) hVarNmin1[12]->Fill(variable[12],variable[18]);
    if (passNHiso2 && passPHiso2 && passSieie2 && passHoe2)   hVarNmin1[13]->Fill(variable[13],variable[18]);
    if (passCHiso2 && passNHiso2 && passSieie2 && passHoe2)   hVarNmin1[14]->Fill(variable[14],variable[18]);
    if (passCHiso2 && passPHiso2 && passSieie2 && passHoe2)   hVarNmin1[15]->Fill(variable[15],variable[18]);

    for (UInt_t z=0; z<NVARIABLES; z++){
      if (z==1 || z==2 || z==8 || z==21 || z==23 || z==25 || z==27){// plot pho1 variables
        if (passAll1) hVarNmin1[z]->Fill(variable[z],variable[18]);
      }
      if (z==9 || z==10 || z==16 || z==22 || z==24 || z==26 || z==28){// plot pho2 variables
        if (passAll2) hVarNmin1[z]->Fill(variable[z],variable[18]);
      }
      if (z==0 || z==17 || z==18 || z==19 || z==20 || z==29){// plot pho1&2 variables
        if (passBoth) hVarNmin1[z]->Fill(variable[z],variable[18]);
      }
    }

    if (variable[0] > 120.0 && variable[0] < 130.0) hMET->Fill(variable[17],variable[18]);
    if (variable[17] > 100.0) hMgg->Fill(variable[0],variable[18]); 
 
   /*
   if (passAny) nphotonsPass++;
    for (int x=0; x<60; x++){
      if (intvariable[29]==x){ // evnts with nvtx = x
        Eff[0][x]++;
        if (passAny) Eff[1][x]++;
      } 
      if ((variable[9]>=10*x && variable[9]<10*(x+1)) || ((variable[1]>=10*x && variable[1]<10*(x+1)))){ //pt bins = 10GeV
        Eff[2][x]++;
        if (passAny) Eff[3][x]++;
      }
      if (((variable[24] >= (-3+(x/10)) && variable[24] < (-3+(x+1)/10)) )||(variable[23] >= (-3+(x/10)) && variable[23]< (-3+(x+1)/10))){ // eta bins = 1/10
        Eff[4][x]++;
        if (passAny) Eff[5][x]++;
      } 
    }// loop over bins for eff plots      
*/   

    phiH[i]=TMath::ATan( (variable[1]*TMath::Sin(variable[21]) - variable[9]*TMath::Sin(variable[22])) / (variable[1]*TMath::Cos(variable[21]) - variable[9]*TMath::Cos(variable[22])) );
    phiHMET[i]=deltaPhi(phiH[i],variable[20]);
    //phiHMET[i]=phiH[i]-variable[20];

    hPhi[0]->Fill(phiH[i],variable[18]);
    hPhi[1]->Fill(variable[20],variable[18]);
    hPhi[2]->Fill(phiHMET[i],variable[18]);


  }// end loop over all photons


  /*Float_t eff[3][60]={0};
  for (int x=0; x<60; x++){
    if(Eff[0][x]!=0) eff[0][x]=(Float_t)Eff[1][x]/(Float_t)Eff[0][x];
    if(Eff[2][x]!=0) eff[1][x]=(Float_t)Eff[3][x]/(Float_t)Eff[2][x];
    if(Eff[4][x]!=0) eff[2][x]=(Float_t)Eff[5][x]/(Float_t)Eff[4][x];
    
    //std::cout << Eff[5][x] << "  " << Eff[4][x] << "  " <<eff[2]<<std::endl;
 
    hEff[0]->Fill(x,eff[0][x]);   
    hEff[1]->Fill(10*x,eff[1][x]);   
    hEff[2]->Fill(-3+(x/10),eff[2][x]);   
  }
  */

  hMgg->GetXaxis()->SetTitle("mgg");
  Plotter::DrawWriteSave1DPlot(hMgg,"mgg_selt1pfMet",true);
  hMET->GetXaxis()->SetTitle("t1pfmet");
  Plotter::DrawWriteSave1DPlot(hMET,"t1pfMet_selMgg",true);

  for (int z=0; z<NVARIABLES; ++z){
    Plotter::DrawWriteSave1DPlot(hVar[z],varname[z],true);
    Plotter::DrawWriteSave1DPlot(hVarNmin1[z],varname[z]+"_n-1",true);
  } 

  hPhi[0]->GetXaxis()->SetTitle("#phi_H");
  hPhi[1]->GetXaxis()->SetTitle("#phi_MET");
  hPhi[2]->GetXaxis()->SetTitle("#Delta#phi(H,MET)");
  Plotter::DrawWriteSave1DPlot(hPhi[0],"phiH",true);
  Plotter::DrawWriteSave1DPlot(hPhi[1],"phiMET",true);
  Plotter::DrawWriteSave1DPlot(hPhi[2],"phiHMET",true);

  /*hEff[0]->GetXaxis()->SetTitle("nvtx");
  hEff[1]->GetXaxis()->SetTitle("p_{T}");
  hEff[2]->GetXaxis()->SetTitle("#eta");
  hEff[0]->GetYaxis()->SetTitle("Efficiency to Pass PhoID");
  hEff[1]->GetYaxis()->SetTitle("Efficiency to Pass PhoID");
  hEff[2]->GetYaxis()->SetTitle("Efficiency to Pass PhoID");
  Plotter::DrawWriteSave1DPlot(hEff[0],"Eff_PHOID_PU",false); 
  Plotter::DrawWriteSave1DPlot(hEff[1],"Eff_PHOID_pt",false); 
  Plotter::DrawWriteSave1DPlot(hEff[2],"Eff_PHOID_eta",false); 
*/
}// end Plotter::make1DHistos

void Plotter::make2DHistos(){

  Int_t range2D[N2DVARIABLES][3]; //nbins,min,max for each 2D variable to plot
  range2D[0][0]= nbins[4]; 
  range2D[0][1]= range[4][0];
  range2D[0][2]= range[4][1];
  range2D[1][0]= nbins[4];   
  range2D[1][1]= range[4][0];
  range2D[1][2]= range[4][1];
  range2D[2][0]= nbins[2];    
  range2D[2][1]= range[2][0];
  range2D[2][2]= range[2][1];
  range2D[3][0]= nbins[2];   
  range2D[3][1]= range[2][0];
  range2D[3][2]= range[2][1];
  range2D[4][0]= nbins[3];   
  range2D[4][1]= range[3][0];
  range2D[4][2]= range[3][1];
  range2D[5][0]= nbins[3];   
  range2D[5][1]= range[3][0];
  range2D[5][2]= range[3][1];
  range2D[6][0]= nbins[5];   
  range2D[6][1]= range[5][0];
  range2D[6][2]= range[5][1];
  range2D[7][0]= nbins[5];   
  range2D[7][1]= range[5][0];
  range2D[7][2]= range[5][1];
  range2D[8][0]= nbins[7];   
  range2D[8][1]= range[7][0];
  range2D[8][2]= range[7][1];
  range2D[9][0]= nbins[7];   
  range2D[9][1]= range[7][0];
  range2D[9][2]= range[7][1];
  range2D[10][0]= nbins[6];   
  range2D[10][1]= range[6][0];
  range2D[10][2]= range[6][1];
  range2D[11][0]= nbins[6];   
  range2D[11][1]= range[6][0];
  range2D[11][2]= range[6][1];


  TH2F *hvPU[N2DVARIABLES];
  TH2F *hvPt[N2DVARIABLES];
  TH2F *hvEta[N2DVARIABLES];    

  for (int z=0; z<N2DVARIABLES; ++z){
    hvPU[z]  = new TH2F(Form("%s_PU_%s",effvar[z].Data(),species.Data()),Form("%s_PU_%s",effvar[z].Data(),species.Data()),60,0,60,range2D[z][0],range2D[z][1],range2D[z][2]);
    hvPt[z]  = new TH2F(Form("%s_pt_%s",effvar[z].Data(),species.Data()),Form("%s_pt_%s",effvar[z].Data(),species.Data()),60,0,600,range2D[z][0],range2D[z][1],range2D[z][2]);
    hvEta[z] = new TH2F(Form("%s_eta_%s",effvar[z].Data(),species.Data()),Form("%s_eta_%s",effvar[z].Data(),species.Data()),60,-3,3,range2D[z][0],range2D[z][1],range2D[z][2]);
  }

  TH2F *hMggvMet = new TH2F(Form("mgg_t1pfMet_%s",species.Data()),Form("mgg_t1pfMet_%s",species.Data()),100,0,1000,50,50,300);

  Float_t var2D[N2DVARIABLES];

  for (int i=0; i<nphotons; ++i){
    tpho->GetEntry(i);
    
    var2D[0]= variable[4];
    var2D[1]= variable[12];
    var2D[2]= variable[2];
    var2D[3]= variable[10];
    var2D[4]= variable[3];
    var2D[5]= variable[11];
    var2D[6]= variable[5];
    var2D[7]= variable[13];
    var2D[8]= variable[7];
    var2D[9]= variable[15];
    var2D[10]= variable[6];
    var2D[11]= variable[14];

    hMggvMet->Fill(variable[17],variable[0],variable[18]);

    for (int z=0; z<N2DVARIABLES; ++z){
      hvPU[z]->Fill(intvariable[29],var2D[z],variable[18]);
      if (z==0 || z== 2 || z==4 || z==6 || z==8 || z==10){// first photon
        hvPt[z]->Fill(variable[1],var2D[z],variable[18]);
        hvEta[z]->Fill(variable[23],var2D[z],variable[18]);
      }
      else{ // second photon
        hvPt[z]->Fill(variable[9],var2D[z],variable[18]);
        hvEta[z]->Fill(variable[24],var2D[z],variable[18]);
      }
    }  
  }


  for (int z=0; z<N2DVARIABLES; ++z){
    hvPU[z]->GetYaxis()->SetTitle(effvar[z].Data());
    hvPt[z]->GetYaxis()->SetTitle(effvar[z].Data());
    hvEta[z]->GetYaxis()->SetTitle(effvar[z].Data());
    hvPU[z]->GetXaxis()->SetTitle("nvtx");
    hvPt[z]->GetXaxis()->SetTitle("p_{T}");
    hvEta[z]->GetXaxis()->SetTitle("#eta");
    Plotter::DrawWriteSave2DPlot(hvPU[z],"PU",effvar[z].Data()); 
    Plotter::DrawWriteSave2DPlot(hvPt[z],"Pt",effvar[z].Data()); 
    Plotter::DrawWriteSave2DPlot(hvEta[z],"Eta",effvar[z].Data()); 
  }

  hMggvMet->GetXaxis()->SetTitle("t1pfMet [GeV]");
  hMggvMet->GetYaxis()->SetTitle("M(#gamma#gamma) [GeV]");
  Plotter::DrawWriteSave2DPlot(hMggvMet,"t1pfMet","mgg");
 
}// end Plotter::make2DHistos


void Plotter::DrawWriteSave1DPlot(TH1F *& h, const TString plotName, const Bool_t DrawNorm){

  fTH1Canv->cd();

  gStyle->SetOptStat(1110);
  gStyle->SetStatBorderSize(0);
  gStyle->SetStatFont(42);
  gStyle->SetStatX(0.84);
  gStyle->SetStatY(0.93);

  h->SetTitle(0);

  fTH1Canv->SetLogy(0);

  if(DrawNorm){
    //if( h->Integral()!=0){ 
    //  h->DrawNormalized();
    //}
    h->Draw();
    Plotter::FindMinAndMax(h,0);
  }
  else{
    h->Draw();
    h->SetMaximum(1.10);
  }

  CMSLumi(fTH1Canv,0,fLumi);
  //Style * cmsLumi = new Style(0.3);
  //cmsLumi->CMSLumi(fTH1Canv,1);
 

//  h->GetXaxis()->SetTitleFont(42);
//  h->GetXaxis()->SetLabelSize(0.03);
//  h->GetXaxis()->SetTitleSize(0.75);

  //std::cout << "Hereiam" << std::endl;
  if (h == (TH1F*) NULL) {std::cout << "WOW" << std::endl;}
  if (outFile == (TFile*) NULL) {std::cout << "WOW" << std::endl;}
  //for (int i = 1; i <= h->GetNbinsX(); i++){
  //  std::cout << h->GetBinContent(i) << std::endl;
  //}
  h->Write();
  //std::cout << "Here" << std::endl;
  if (fTH1Canv == (TCanvas*) NULL) {std::cout << "WOW" << std::endl;}
  //std::cout << Form("%s%s/%s_%s.png",fName.Data(),species.Data(),plotName.Data(),species.Data()) << std::endl;
  fTH1Canv->cd();
  fTH1Canv->SaveAs(Form("%s%s/%s_%s.png",fName.Data(),species.Data(),plotName.Data(),species.Data()));
//  delete cmsLumi;

  //std::cout << "Here" << std::endl;

  fTH1Canv->SetLogy(1);
  if (DrawNorm){
    if (h->Integral() != 0){ h->DrawNormalized(); }
    else { h->Draw(); }
    Plotter::FindMinAndMax(h,1);
  }
  else{
    h->Draw();
    h->SetMaximum(1.1*h->GetMaximum());
  }
  h->Write();
  fTH1Canv->SaveAs(Form("%s%s/%s_%s_log.png",fName.Data(),species.Data(),plotName.Data(),species.Data()));
}// end Plotter::DrawWriteSave1DPlot


void Plotter::DrawWriteSave2DPlot(TH2F *& h, const TString varX, const TString varY){
  gStyle->SetOptStat(0);
  fTH2Canv->cd();
  fTH2Canv->SetLogy(0);
  h->Draw("colz");
  h->SetTitle(0);
  h->Write();
  //Plotter::CMS_Lumi(0.3)
  fTH2Canv->SaveAs(Form("%s%s/%s_%s_%s.png",fName.Data(),species.Data(),varY.Data(),varX.Data(),species.Data()));

  // FIXME problem with log plots (min always is zero)
/*  fTH2Canv->SetLogy(1);
  fTH2Canv->SetLogx(1);  
  h->Draw("colz");
  h->Write();
  fTH2Canv->SaveAs(Form("%s%s/%s_%s_%s_log.png",fName.Data(),species.Data(),varY.Data(),varX.Data(),species.Data()));*/
}// end Plotter::DrawWriteSave2DPlot

void Plotter::FindMinAndMax(TH1F *& h, int plotLog){
  Float_t max = h->GetMaximum();
  if (plotLog==1) h->SetMaximum(10*max);
  if (plotLog==0) h->SetMaximum(2*max);

  Float_t min = 1000;
  Bool_t newmin = false;

  for (Int_t bin=1; bin <= h->GetNbinsX(); bin++){
    Float_t tmpmin = h->GetBinContent(bin);
    if ((tmpmin < min) && (tmpmin > 0)){
      min = tmpmin;
      newmin = true;
    }
  }

  if (newmin){
    h->SetMinimum(0.90*min);
  }
}// end Plotter::FindMinAndMax






void Plotter::InitTreeVar(){
  varname.push_back("mgg");
  varname.push_back("pt1");
  varname.push_back("r91");
  varname.push_back("sieie1");
  varname.push_back("hoe1");
  varname.push_back("chiso1");
  varname.push_back("phoiso1");
  varname.push_back("neuiso1");
  varname.push_back("eleveto1");
  varname.push_back("pt2");
  varname.push_back("r92");
  varname.push_back("sieie2");
  varname.push_back("hoe2");
  varname.push_back("chiso2");
  varname.push_back("phoiso2");
  varname.push_back("neuiso2");
  varname.push_back("eleveto2");
  varname.push_back("t1pfmet");
  varname.push_back("weight");
  varname.push_back("ptgg");
  varname.push_back("t1pfmetPhi");
  varname.push_back("phi1");
  varname.push_back("phi2");
  varname.push_back("eta1");
  varname.push_back("eta2");
  varname.push_back("presel1");
  varname.push_back("presel2");
  varname.push_back("sel1");
  varname.push_back("sel2");
  varname.push_back("nvtx");

  nbins[0]=50; 		// mgg
  nbins[1]=50; 		// pt1
  nbins[2]=100; 	// r91
  nbins[3]=300; 	// sieie1
  nbins[4]=250; 	// hoe1
  nbins[5]=150; 	// chiso1
  nbins[6]=150; 	// phoiso1
  nbins[7]=150; 	// neuiso1
  nbins[8]=100; 	// eleveto1
  nbins[9]=nbins[1];  	// pt2
  nbins[10]=nbins[2]; 	// r92
  nbins[11]=nbins[3]; 	// sieie2
  nbins[12]=nbins[4]; 	// hoe2
  nbins[13]=nbins[5]; 	// chiso2
  nbins[14]=nbins[6]; 	// phoiso2
  nbins[15]=nbins[7]; 	// neuiso2
  nbins[16]=nbins[8]; 	// eleveto2
  nbins[17]=100;     	// t1pfmet
  nbins[18]=100; 	// weight
  nbins[19]=100;	// ptgg
  nbins[20]=80;		// t1pfmetphi
  nbins[21]=nbins[20];	// phi1
  nbins[22]=nbins[20];  // phi2
  nbins[23]=100;	// eta1
  nbins[24]=nbins[23];	// eta2
  nbins[25]=2;		// presel1
  nbins[26]=nbins[25];	// presel2
  nbins[27]=nbins[25];	// sel1
  nbins[28]=nbins[25];	// sel2
  nbins[29]=60;		// nvtx

  range[0][0]=50.; 		// mgg
  range[1][0]=0.; 		// pt1
  range[2][0]=0.; 		// r91
  range[3][0]=0.; 		// sieie1
  range[4][0]=0.; 		// hoe1
  range[5][0]=-5.; 		// chiso1
  range[6][0]=-5.; 		// phoiso1
  range[7][0]=-5.; 		// neuiso1
  range[8][0]=-1.; 		// eleveto1
  range[9][0]=range[1][0]; 	// pt2
  range[10][0]=range[2][0]; 	// r92
  range[11][0]=range[3][0]; 	// sieie2
  range[12][0]=range[4][0]; 	// hoe2
  range[13][0]=range[5][0]; 	// chiso2
  range[14][0]=range[6][0]; 	// phoiso2
  range[15][0]=range[7][0]; 	// neuiso2
  range[16][0]=range[8][0]; 	// eleveto2
  range[17][0]=0.; 		// t1pfmet
  range[18][0]=0.; 		// weight
  range[19][0]=0.;		// ptgg
  range[20][0]=-4.;		// t1pfmetphi
  range[21][0]=range[20][0];	// phi1
  range[22][0]=range[20][0];	// phi2
  range[23][0]=-5.;		// eta1
  range[24][0]=range[23][0];	// eta2
  range[25][0]=0.;		// presel1
  range[26][0]=range[25][0];	// presel2
  range[27][0]=range[25][0];	// sel1
  range[28][0]=range[25][0];	// sel2
  range[29][0]=0.;		// nvtx

  range[0][1]=300.; 		// mgg
  range[1][1]=500.; 		// pt1
  range[2][1]=1.1; 		// r91
  range[3][1]=0.03; 		// sieie1
  range[4][1]=0.025; 		// hoe1
  range[5][1]=15.; 		// chiso1
  range[6][1]=15.; 		// phoiso1
  range[7][1]=15.; 		// neuiso1
  range[8][1]=1.; 		// eleveto1
  range[9][1]=range[1][1]; 	// pt2
  range[10][1]=range[2][1]; 	// r92
  range[11][1]=range[3][1]; 	// sieie2
  range[12][1]=range[4][1]; 	// hoe2
  range[13][1]=range[5][1]; 	// chiso2
  range[14][1]=range[6][1]; 	// phoiso2
  range[15][1]=range[7][1]; 	// neuiso2
  range[16][1]=range[8][1]; 	// eleveto2
  range[17][1]=1000.; 		// t1pfmet
  range[18][1]=100.;		// weight
  range[19][1]=1000.;		// ptgg
  range[20][1]=4.;		// t1pfmetphi
  range[21][1]=range[20][1];	// phi1
  range[22][1]=range[20][1];	// phi2
  range[23][1]=-5.;		// eta1
  range[24][1]=range[23][1];	// eta2
  range[25][1]=2.;		// presel1
  range[26][1]=range[25][1];	// presel2
  range[27][1]=range[25][1];	// sel1
  range[28][1]=range[25][1];	// sel2
  range[29][1]=60.;		// nvtx

  xaxisLabel.push_back("m(#gamma#gamma) [GeV]");
  xaxisLabel.push_back("p_{T}(#gamma1) [GeV]");
  xaxisLabel.push_back("R9(#gamma1)");
  xaxisLabel.push_back("#sigma_{i#eta i#eta}(#gamma1)");
  xaxisLabel.push_back("H/E(#gamma1)");
  xaxisLabel.push_back("CHiso(#gamma1)");
  xaxisLabel.push_back("PHOiso(#gamma1)");
  xaxisLabel.push_back("NEUiso(#gamma1)");
  xaxisLabel.push_back("ele veto(#gamma1)");
  xaxisLabel.push_back("p_{T}(#gamma2) [GeV]");
  xaxisLabel.push_back("R9(#gamma2)");
  xaxisLabel.push_back("#sigma_{i#eta i#eta}(#gamma2)");
  xaxisLabel.push_back("H/E(#gamma2)");
  xaxisLabel.push_back("CHiso(#gamma2)");
  xaxisLabel.push_back("PHOiso(#gamma2)");
  xaxisLabel.push_back("NEUiso(#gamma2)");
  xaxisLabel.push_back("ele veto(#gamma2)");
  xaxisLabel.push_back("type 1 PF MET(#gamma #gamma)");
  xaxisLabel.push_back("weight");  
  xaxisLabel.push_back("p_{T}(#gamma#gamma)");
  xaxisLabel.push_back("type 1 PF MET #phi");
  xaxisLabel.push_back("#phi(#gamma1)");
  xaxisLabel.push_back("#phi(#gamma2)");
  xaxisLabel.push_back("#eta(#gamma1)");
  xaxisLabel.push_back("#eta(#gamma2)");
  xaxisLabel.push_back("presel(#gamma1)");
  xaxisLabel.push_back("presel(#gamma2)");
  xaxisLabel.push_back("sel(#gamma1)");
  xaxisLabel.push_back("sel(#gamma2)");
  xaxisLabel.push_back("nvtx");

}// end Plotter::InitTreeVar

void Plotter::InitTreeEffVar(){
  effvar.push_back("hoe1");
  effvar.push_back("hoe2");
  effvar.push_back("r91");
  effvar.push_back("r92");
  effvar.push_back("sieie1");
  effvar.push_back("sieie2");
  effvar.push_back("chiso1");
  effvar.push_back("chiso2");
  effvar.push_back("neuiso1");
  effvar.push_back("neuiso2");
  effvar.push_back("phoiso1");
  effvar.push_back("phoiso2");
}// end Plotter::InitTreeEffVar

void Plotter::InitPhotonIDSel(){
  selvar.push_back("passCHiso1"); 
  selvar.push_back("passCHiso2");
  selvar.push_back("passNHiso1"); 
  selvar.push_back("passNHiso2");
  selvar.push_back("passPHiso1"); 
  selvar.push_back("passPHiso2");
  selvar.push_back("passSieie1");
  selvar.push_back("passSieie2");
  selvar.push_back("passHoe1");
  selvar.push_back("passHoe2");
}// end Plotter::InitPhotonIDSel




