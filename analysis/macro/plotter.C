#include "TMath.h"
#include "TTree.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TFile.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "mkPlotsLivia/CMS_lumi.h"
#include "mkPlotsLivia/CMS_lumi.C"

#include <iostream>

#define NSPECIES 1
#define NVARIABLES 18

void plotter( char * name ){
//  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1110);
  gStyle->SetStatBorderSize(0);
  gStyle->SetStatFont(42);
  gStyle->SetStatX(0.84);  
  gStyle->SetStatY(0.93);
//  gROOT->ProcessLine(".x ./DiphotonStyle.C");

  TString suffix[NSPECIES];
  suffix[0]=name;

  TString infile = "50kSamples";

  // open file and tree
  TFile *f = new TFile("data/"+infile+"/diPhotons_"+suffix[0]+".root");
  TTree *tev = (TTree*)f->Get("diPhoAna");
  TTree *tpho = (TTree*)f->Get("diPhoAna/DiPhotonTree");
  // create output file
  TFile *fOut = new TFile("diPhoPlots/"+infile+"/"+suffix[0]+"/diPhotHistos_"+suffix[0]+".root","RECREATE");
  
  // variables
  float variable[NVARIABLES];
  int intvariable[NVARIABLES];
  
  // initialize the variables
  for(int i=0; i<NVARIABLES; i++){
    variable[i] = -1000.;
    intvariable[i] = -1000;
  }
  
  TString varname[NVARIABLES];
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

  int nbins[NVARIABLES];
  nbins[0]=60; 		// mgg
  nbins[1]=50; 		// pt1
  nbins[2]=100; 	// r91
  nbins[3]=300; 	// sieie1
  nbins[4]=250; 	// hoe1
  nbins[5]=200; 	// chiso1
  nbins[6]=200; 	// phoiso1
  nbins[7]=200; 	// neuiso1
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

  float min[NVARIABLES];
  min[0]=0.; 		// mgg
  min[1]=0.; 		// pt1
  min[2]=0.; 		// r91
  min[3]=0.; 		// sieie1
  min[4]=0.; 		// hoe1
  min[5]=-10.; 		// chiso1
  min[6]=-10.; 		// phoiso1
  min[7]=-10.; 		// neuiso1
  min[8]=-1.; 		// eleveto1
  min[9]=min[1]; 	// pt2
  min[10]=min[2]; 	// r92
  min[11]=min[3]; 	// sieie2
  min[12]=min[4]; 	// hoe2
  min[13]=min[5]; 	// chiso2
  min[14]=min[6]; 	// phoiso2
  min[15]=min[7]; 	// neuiso2
  min[16]=min[8]; 	// eleveto2
  min[17]=0.; 		// t1pfmet

  float max[NVARIABLES];
  max[0]=300.; 		// mgg
  max[1]=500.; 		// pt1
  max[2]=1.1; 		// r91
  max[3]=0.03;	 	// sieie1
  max[4]=0.025;	 	// hoe1
  max[5]=10.; 		// chiso1
  max[6]=10.; 		// phoiso1
  max[7]=10.; 		// neuiso1
  max[8]=1.; 		// eleveto1
  max[9]=max[1]; 	// pt2
  max[10]=max[2]; 	// r92
  max[11]=max[3]; 	// sieie2
  max[12]=max[4]; 	// hoe2
  max[13]=max[5]; 	// chiso2
  max[14]=max[6]; 	// phoiso2
  max[15]=max[7]; 	// neuiso2
  max[16]=max[8]; 	// eleveto2
  max[17]=1000.; 	// t1pfmet

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
  


  // make histograms
  TH1F *h[NVARIABLES];
  for (int z=0; z<NVARIABLES; ++z){
     h[z] = new TH1F(varname[z]+"_"+suffix[0],varname[z]+"_"+suffix[0],nbins[z],min[z],max[z]);
     if(z==8 || z==16 ){
       tpho->SetBranchAddress(varname[z],&intvariable[z]);
     }
     else{
       tpho->SetBranchAddress(varname[z],&variable[z]);
     }
  }// one histogram for each variable loop

/*  TH1F *isoEB[6];
  isoEB[0]=new TH1F(variable[5]+"_ptadjustEB_"+suffix[0],varname[5]+"_"+suffix[0],nbins[5],min[5],max[5]);
  isoEB[1]=new TH1F(variable[6]+"_ptadjustEB_"+suffix[0],varname[6]+"_"+suffix[0],nbins[6],min[6],max[6]);
  isoEB[2]=new TH1F(variable[7]+"_ptadjustEB_"+suffix[0],varname[7]+"_"+suffix[0],nbins[7],min[7],max[7]);
  isoEB[3]=new TH1F(variable[13]+"_ptadjustEB_"+suffix[0],varname[13]+"_"+suffix[0],nbins[13],min[13],max[13]);
  isoEB[4]=new TH1F(variable[14]+"_ptadjustEB_"+suffix[0],varname[14]+"_"+suffix[0],nbins[14],min[14],max[14]);
  isoEB[5]=new TH1F(variable[15]+"_ptadjustEB_"+suffix[0],varname[15]+"_"+suffix[0],nbins[15],min[15],max[15]);
  TH1F *isoEE[6];
  isoEE[0]=new TH1F(variable[5]+"_ptadjustEE_"+suffix[0],varname[5]+"_"+suffix[0],nbins[5],min[5],max[5]);
  isoEE[1]=new TH1F(variable[6]+"_ptadjustEE_"+suffix[0],varname[6]+"_"+suffix[0],nbins[6],min[6],max[6]);
  isoEE[2]=new TH1F(variable[7]+"_ptadjustEE_"+suffix[0],varname[7]+"_"+suffix[0],nbins[7],min[7],max[7]);
  isoEE[3]=new TH1F(variable[13]+"_ptadjustEE_"+suffix[0],varname[13]+"_"+suffix[0],nbins[13],min[13],max[13]);
  isoEE[4]=new TH1F(variable[14]+"_ptadjustEE_"+suffix[0],varname[14]+"_"+suffix[0],nbins[14],min[14],max[14]);
  isoEE[5]=new TH1F(variable[15]+"_ptadjustEE_"+suffix[0],varname[15]+"_"+suffix[0],nbins[15],min[15],max[15]);
*/
 
  int nphotons = (int)tpho->GetEntries();
  for (int i=0; i<nphotons; i++){
    tpho->GetEntry(i);
    if (variable[1] > variable[0]/3 && variable[9] > variable[0]/4)
    { // pt1 > mgg/3 & pt2 > mgg/4 selection
      for (int z=0; z<NVARIABLES; z++){
        if (z==8 || z==16){ // eleveto1 & eleveto2
          h[z]->Fill(intvariable[z]);
        }
        else{
          h[z]->Fill(variable[z]);
        }

	// do pt adjustment from the pho id cut
        /*if (z==5){ // chiso1
	}
	if (z==6){ // neuiso1
        }
	if (z==7){ // phoiso1
        }
	if (z==13){ // chiso2
        }
	if (z==14){ // neuiso2
        }
	if (z==15){ // phoiso2
        }*/

      }// loop over variables
    }// loose pt selection
  }// loop over photons

  double max1;

  for (int z=0; z<NVARIABLES; z++){
    fOut->cd();
    //make two TCanvas one regular one logy 
    TCanvas* c1 = new TCanvas("c1","",1200,800);
    c1->SetLogy(0);
    h[z]->DrawNormalized();
    max1 = h[z]->GetMaximum();
    h[z]->SetMaximum(10*max1);
    h[z]->GetXaxis()->SetTitle(xaxisLabel[z]);
    //h[z]->GetXaxis()->SetTitleSize(0.75);
    //h[z]->GetXaxis()->SetTitleFont(42);
    CMS_lumi( (TPad*)c1->cd(),true,0);   
 
    TCanvas* c2 = new TCanvas("c2","",1200,800);
    c2->SetLogy(1);
    h[z]->DrawNormalized();
    h[z]->SetMaximum(10*max1);
    h[z]->GetXaxis()->SetTitle(xaxisLabel[z]);
    //h[z]->GetXaxis()->SetTitleSize(30);
    //h[z]->GetXaxis()->SetTitleFont(42);
    CMS_lumi( (TPad*)c2->cd(),true,0);
   

    c1->SaveAs("diPhoPlots/"+infile+"/"+suffix[0]+"/"+varname[z]+"_"+suffix[0]+".png");
    c2->SaveAs("diPhoPlots/"+infile+"/"+suffix[0]+"/"+varname[z]+"_"+suffix[0]+"_log.png");
  }// loop over variables


  fOut->Write();
  fOut->Close();
}




