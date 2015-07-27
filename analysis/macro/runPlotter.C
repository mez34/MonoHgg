void runPlotter(){
  gROOT->LoadMacro("Plotter.cpp++g");
  
  // Call this script:
  // root -l -b -q 'runPlotter()'
  //
  // Arguments to Plotter:
  // 1st : location of input data
  // 2nd : output data location
  // 3rd : name of sample 
 
  Plotter test("data/ALL_nosel/diPhotons","diPhoPlots/ALL_nosel/","DMHtoGG");
  test.getTree();
  test.make1DHistos();
  test.make2DHistos();
  test.Fill1DHistos();
}
