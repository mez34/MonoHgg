#include "Plotter.hh"
#include "mkPlotsLivia/CMS_lumi.C"

Plotter::Plotter( TString inName, TString outName, TString inSpecies){


  // Get ROOT file
  name = inName;
  species = inSpecies;
  inFile = TFile::Open(Form("%s_%s_nosel.root",name.Data(),species.Data()));
  
  fName = outName;
  // Make output directory
  FileStat_t dummyFileState;
  TString FullPath = fName.Data();
  FullPath+=species.Data();
  FullPath+="/";
  if (gSystem->GetPathInfo(FullPath.Data(), dummyFileState) == 1){
    TString mkDir = "mkdir -p ";
    mkDir += FullPath.Data();
    gSystem->Exec(mkDir.Data());
  }

  // Make output ROOT file
  outFile = new TFile(Form("%s/%s/plots_%s.root",fName.Data(),species.Data(),species.Data()),"RECREATE");

  // Make TCanvas
  fTH1Canv = new TCanvas();
  fTH2Canv = new TCanvas();

  NVARIABLES = 30;
  N2DVARIABLES = 12;

}// end Plotter::Plotter

Plotter::~Plotter(){
  // Write and Close output ROOT file
  outFile->Write();
  outFile->Close(); 

  std::cout<<"Finished & Deleting"<<std::endl;
  delete inFile;
  delete outFile;
  delete fTH1Canv;
  delete fTH2Canv;
}// end Plotter::~Plotter

void Plotter::getTree(){
  // Open Tree from inFile
  tpho = (TTree*)inFile->Get("DiPhotonTree"); 

  // Load variables from Tree
  variable[NVARIABLES]    = {-1000}; // float for most variables 
  intvariable[NVARIABLES] = {-1000}; // int for other variables (eleveto, sel, nvtx)

  varname[0]="mgg";
  varname[1]="pt1";
  varname[2]="r91";
  varname[3]="sieie1";
  varname[4]="hoe1";
  varname[5]="chiso1";
  varname[6]="phoiso1";
  varname[7]="neuiso1";
  varname[8]="eleveto1";
  varname[9]="pt2";
  varname[10]="r92";
  varname[11]="sieie2";
  varname[12]="hoe2";
  varname[13]="chiso2";
  varname[14]="phoiso2";
  varname[15]="neuiso2";
  varname[16]="eleveto2";
  varname[17]="t1pfmet";
  varname[18]="weight";
  varname[19]="ptgg";
  varname[20]="t1pfmetPhi";
  varname[21]="phi1";
  varname[22]="phi2";
  varname[23]="eta1";
  varname[24]="eta2";
  varname[25]="presel1";
  varname[26]="presel2";
  varname[27]="sel1";
  varname[28]="sel2";
  varname[29]="nvtx";

  for(int z=0; z<NVARIABLES; ++z){
    if(z==8 || z==16 || z==25 || z==26 || z==27 || z==28 || z==29){ //eleveto, sel, nvtx
       tpho->SetBranchAddress(varname[z],&intvariable[z]);
     }
     else{
       tpho->SetBranchAddress(varname[z],&variable[z]);
     }
  }
 
  nphotons = (int)tpho->GetEntries();

}// end Plotter::getTree


void Plotter::make1DHistos(){

  nbins[0]=60; 		// mgg
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

  range[0][0]=0.; 		// mgg
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

  TString xaxisLabel[NVARIABLES];
  xaxisLabel[0]="m(#gamma#gamma) [GeV]";
  xaxisLabel[1]="p_{T}(#gamma1) [GeV]";
  xaxisLabel[2]="R9(#gamma1)";
  xaxisLabel[3]="#sigma_{i#eta i#eta}(#gamma1)";
  xaxisLabel[4]="H/E(#gamma1)";
  xaxisLabel[5]="CHiso(#gamma1)";
  xaxisLabel[6]="PHOiso(#gamma1)";
  xaxisLabel[7]="NEUiso(#gamma1)";
  xaxisLabel[8]="ele veto(#gamma1)";
  xaxisLabel[9]="p_{T}(#gamma2) [GeV]";
  xaxisLabel[10]="R9(#gamma2)";
  xaxisLabel[11]="#sigma_{i#eta i#eta}(#gamma2)";
  xaxisLabel[12]="H/E(#gamma2)";
  xaxisLabel[13]="CHiso(#gamma2)";
  xaxisLabel[14]="PHOiso(#gamma2)";
  xaxisLabel[15]="NEUiso(#gamma2)";
  xaxisLabel[16]="ele veto(#gamma2)";
  xaxisLabel[17]="type 1 PF MET(#gamma #gamma)";
  xaxisLabel[18]="weight";  
  xaxisLabel[19]="p_{T}(#gamma#gamma)";
  xaxisLabel[20]="type 1 PF MET #phi";
  xaxisLabel[21]="#phi(#gamma1)";
  xaxisLabel[22]="#phi(#gamma2)";
  xaxisLabel[23]="#eta(#gamma1)";
  xaxisLabel[24]="#eta(#gamma2)";
  xaxisLabel[25]="presel(#gamma1)";
  xaxisLabel[26]="presel(#gamma2)";
  xaxisLabel[27]="sel(#gamma1)";
  xaxisLabel[28]="sel(#gamma2)";
  xaxisLabel[29]="nvtx";


  TH1F *hVar[NVARIABLES];
  for (int z=0; z<NVARIABLES; ++z){
    hVar[z] = new TH1F(Form("%s_%s",varname[z].Data(),species.Data()),Form("%s_%s",varname[z].Data(),species.Data()),nbins[z],range[z][0],range[z][1]);
    hVar[z]->GetXaxis()->SetTitle(xaxisLabel[z]);
  }// one histogram for each variable in Tree 

  // additional histograms
  TH1F *hPhi[3]; // phi of the higgs, MET and delta phi H,MET
  hPhi[0] = new TH1F(Form("phi_H_%s",species.Data()),Form("phi_H_%s",species.Data()),80,-4,4);
  hPhi[1] = new TH1F(Form("phi_MET_%s",species.Data()),Form("phi_MET_%s",species.Data()),80,-4,4);
  hPhi[2] = new TH1F(Form("phi_HMET_%s",species.Data()),Form("phi_HMET_%s",species.Data()),80,-4,4);

  TH1F *hEff[3]; // eff of pho ID vs PU, pt, eta
  hEff[0] = new TH1F(Form("Eff_PHID_PU_%s",species.Data()),Form("Eff_PHID_PU_%s",species.Data()),60,0,60);
  hEff[1] = new TH1F(Form("Eff_PHID_pt_%s",species.Data()),Form("Eff_PHID_pt_%s",species.Data()),60,0,600);
  hEff[2] = new TH1F(Form("Eff_PHID_eta_%s",species.Data()),Form("Eff_PHID_eta_%s",species.Data()),60,-3,3);


  Float_t phiH[nphotons]    = {-1000}; 
  Float_t phiHMET[nphotons] = {-1000};

  for (int i=0; i<nphotons; ++i){
    tpho->GetEntry(i);
    for (int z=0; z<NVARIABLES; ++z){
      hVar[z]->Fill(variable[z],variable[18]);
    }// end loop over the variables in the Tree
    phiH[i]=TMath::ATan( (variable[1]*TMath::Sin(variable[21]) - variable[9]*TMath::Sin(variable[22])) / (variable[1]*TMath::Cos(variable[21]) - variable[9]*TMath::Cos(variable[22])) );
    phiHMET[i]=phiH[i]-variable[20];

    hPhi[0]->Fill(phiH[i],variable[18]);
    hPhi[0]->GetXaxis()->SetTitle("#phi_H");
    hPhi[1]->Fill(variable[20],variable[18]);
    hPhi[1]->GetXaxis()->SetTitle("#phi_MET");
    hPhi[2]->Fill(phiHMET[i],variable[18]);
    hPhi[2]->GetXaxis()->SetTitle("#Delta#phi(H,MET)");

  }// end loop over all photons
 
  for (int z=0; z<NVARIABLES; ++z){
    Plotter::DrawWriteSave1DPlot(hVar[z],varname[z]);
  } 
  Plotter::DrawWriteSave1DPlot(hPhi[0],"phiH");
  Plotter::DrawWriteSave1DPlot(hPhi[1],"phiMET");
  Plotter::DrawWriteSave1DPlot(hPhi[2],"phiHMET");

}// end Plotter::make1DHistos

void Plotter::make2DHistos(){

  Float_t evar[N2DVARIABLES];
  
  TString effvar[N2DVARIABLES];
  effvar[0]="hoe1";
  effvar[1]="hoe2";
  effvar[2]="r91";
  effvar[3]="r92";
  effvar[4]="sieie1";
  effvar[5]="sieie2";
  effvar[6]="chiso1";
  effvar[7]="chiso2";
  effvar[8]="neuiso1";
  effvar[9]="neuiso2";
  effvar[10]="phoiso1";
  effvar[11]="phoiso2";

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
    hvPt[z]  = new TH2F(Form("%s_pt_%s",effvar[z].Data(),species.Data()),Form("%s_pt_%s",effvar[z].Data(),species.Data()),60,0,60,range2D[z][0],range2D[z][1],range2D[z][2]);
    hvEta[z] = new TH2F(Form("%s_eta_%s",effvar[z].Data(),species.Data()),Form("%s_eta_%s",effvar[z].Data(),species.Data()),60,0,60,range2D[z][0],range2D[z][1],range2D[z][2]);
  }

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

    for (int z=0; z<N2DVARIABLES; ++i){
      hvPU[z]->Fill(variable[29],var2D[z],variable[18]);
      if (z==0 || z== 2 || z==4 || z==6 || z==8 || z==10){// first photon
        hvPt[z]->Fill(variable[1],var2D[z],variable[18]);
        hvEta[z]->Fill(variable[23],var2D[z],variable[18]);
      }
      else{// second photon
        hvPt[z]->Fill(variable[9],var2D[z],variable[18]);
        hvEta[z]->Fill(variable[24],var2D[z],variable[18]);
      }
    }  
  }


  for (int z=0; z<N2DVARIABLES; ++z){
    Plotter::DrawWriteSave2DPlot(hvPU[z],"PU",effvar[z]); 
    Plotter::DrawWriteSave2DPlot(hvPt[z],"Pt",effvar[z]); 
    Plotter::DrawWriteSave2DPlot(hvEta[z],"Eta",effvar[z]); 
  } 
}// end Plotter::make2DHistos


void Plotter::Fill1DHistos(){
  for (Int_t i=0; i<nphotons; i++){
    for (Int_t z=0; z<NVARIABLES; z++){
      //hVar[z]->Fill(variable[z],variable[18]); //Fill by weight
    } 
  }
}// end Plotter::Fill1DHistos

void Plotter::DrawWriteSave1DPlot(TH1F *& h, const TString plotName){
  fTH1Canv->cd();
  fTH1Canv->SetLogy(0);
  h->DrawNormalized();
  Plotter::FindMinAndMax(h,0);
  CMS_lumi( (TPad*)fTH1Canv->cd(),true,0);
  h->Write();
  fTH1Canv->SaveAs(Form("%s%s/%s_%s.png",fName.Data(),species.Data(),plotName.Data(),species.Data()));

  fTH1Canv->SetLogy(1);
  Plotter::FindMinAndMax(h,1);
  h->DrawNormalized();
  CMS_lumi( (TPad*)fTH1Canv->cd(),true,0);
  h->Write();
  fTH1Canv->SaveAs(Form("%s%s/%s_%s_log.png",fName.Data(),species.Data(),plotName.Data(),species.Data()));
}// end Plotter::DrawWriteSave1DPlot


void Plotter::DrawWriteSave2DPlot(TH2F *& h, const TString varX, const TString varY){
  fTH2Canv->cd();
  fTH2Canv->SetLogy(0);
  h->Draw("colz");
  CMS_lumi( (TPad*)fTH2Canv->cd(),true,0);
  h->Write();
  fTH2Canv->SaveAs(Form("%s%s/%s_%s_%s.png",fName.Data(),species.Data(),varY.Data(),varX.Data(),species.Data()));

  fTH2Canv->SetLogy(1);
  h->Draw("colz");
  CMS_lumi( (TPad*)fTH2Canv->cd(),true,0);
  h->Write();
  fTH2Canv->SaveAs(Form("%s%s/%s_%s_%s_log.png",fName.Data(),species.Data(),varY.Data(),varX.Data(),species.Data()));
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


