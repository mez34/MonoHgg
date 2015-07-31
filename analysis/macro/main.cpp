#include "Plotter.hh"
#include "TROOT.h"

#include <iostream>


int main(){
  //gROOT->LoadMacro("Plotter.cpp++g");
  
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
 
//  Plotter test1("./data/50ns/","./diPhoPlots/50ns/","GJets");

  Plotter * test1 = new Plotter("./data/ALL_nosel/diPhotons","./diPhoPlots/ALL_nosel/","DMHtoGG",30);
  test1->DoPlots();
  delete test1;

/*  Plotter test2("./data/50ns/","./diPhoPlots/50ns/","QCD");
  test2.getTree();
  test2.make1DHistos();
  test2.make2DHistos();
  test2.Fill1DHistos();
*/

}
