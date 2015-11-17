#! /bin/sh                                                                                                    

lumi=$1      # in pb
echo "Adding weights for " $lumi " pb-1..."

root -l -b <<EOF
.L addWeightsToTree.cc++  
 
addWeights("data/25ns_v7_noEV_3fb/GJet_Pt-20to40.root", $lumi);
addWeights("data/25ns_v7_noEV_3fb/GJet_Pt-40toInf.root", $lumi);
addWeights("data/25ns_v7_noEV_3fb/QCD_Pt-30to40.root", $lumi);
addWeights("data/25ns_v7_noEV_3fb/QCD_Pt-30toInf.root", $lumi);
addWeights("data/25ns_v7_noEV_3fb/QCD_Pt-40toInf.root", $lumi);
addWeights("data/25ns_v7_noEV_3fb/GluGluHToGG.root", $lumi);
addWeights("data/25ns_v7_noEV_3fb/DiPhoton.root", $lumi);
addWeights("data/25ns_v7_noEV_3fb/VH.root", $lumi);
addWeights("data/25ns_v7_noEV_3fb/DYJetsToLL.root", $lumi);

addWeights("data/25ns_v7_noEV_3fb/2HDM_mZP600.root", $lumi);
addWeights("data/25ns_v7_noEV_3fb/2HDM_mZP800.root", $lumi);
addWeights("data/25ns_v7_noEV_3fb/2HDM_mZP1000.root", $lumi);
addWeights("data/25ns_v7_noEV_3fb/2HDM_mZP1200.root", $lumi);
addWeights("data/25ns_v7_noEV_3fb/2HDM_mZP1400.root", $lumi);
addWeights("data/25ns_v7_noEV_3fb/2HDM_mZP1700.root", $lumi);
addWeights("data/25ns_v7_noEV_3fb/2HDM_mZP2000.root", $lumi);
addWeights("data/25ns_v7_noEV_3fb/2HDM_mZP2500.root", $lumi);

addWeights("data/25ns_v7_noEV_3fb/DoubleEG_p.root", $lumi);
addWeights("data/25ns_v7_noEV_3fb/DoubleEG_0.root", $lumi);
addWeights("data/25ns_v7_noEV_3fb/DoubleEG_1.root", $lumi);
addWeights("data/25ns_v7_noEV_3fb/DoubleEG_2.root", $lumi);

.q

EOF

echo "done weighting."

#addWeights("data/50ns_betaV4/DoubleEG.root", $lumi);

#previous betaV2 version Z,W+,W- are separate:
#addWeights("data/50ns/ZH.root", $lumi);
#addWeights("data/50ns/WplusH.root", $lumi);
#addWeights("data/50ns/WminusH.root", $lumi);

#signal samples weren't remade 
#addWeights("data/50ns/DMHtoGG_M1000.root", $lumi);
#addWeights("data/50ns/DMHtoGG_M100.root", $lumi);
#addWeights("data/50ns/DMHtoGG_M10.root", $lumi);
#addWeights("data/50ns/DMHtoGG_M1.root", $lumi);



#addWeights("data/merged/GGJets_M-200To500.root", $lumi);
#addWeights("data/merged/GGJets_M-500To1000.root", $lumi);
#addWeights("data/merged/GGJets_M-1000To2000.root", $lumi);
#addWeights("data/merged/GGJets_M-2000To4000.root", $lumi);
#addWeights("data/merged/GGJets_M-4000To8000.root", $lumi);
#addWeights("data/merged/GGJets_M-8000To13000.root", $lumi);
#addWeights("data/merged/GJets_HT-100to200.root", $lumi);
#addWeights("data/merged/GJets_HT-200to400.root", $lumi);
#addWeights("data/merged/GJets_HT-400to600.root", $lumi);
#addWeights("data/merged/GJets_HT-600toInf.root", $lumi);
#addWeights("data/merged/QCD_HT-100To250.root", $lumi);    
#addWeights("data/merged/QCD_HT-250To500.root", $lumi);   
#addWeights("data/merged/QCD_HT-500To1000.root", $lumi);   
#addWeights("data/merged/QCD_HT-1000ToInf.root", $lumi);   
#addWeights("data/merged/RSGravToGG_kMpl-001_M-750.root",  $lumi, 750);
#addWeights("data/merged/RSGravToGG_kMpl-001_M-1500.root", $lumi, 1500);
#addWeights("data/merged/RSGravToGG_kMpl-001_M-5000.root", $lumi, 5000);
#addWeights("data/merged/RSGravToGG_kMpl-01_M-1500.root", $lumi, 1500);
#addWeights("data/merged/RSGravToGG_kMpl-01_M-3000.root", $lumi, 3000);
#addWeights("data/merged/RSGravToGG_kMpl-02_M-1500.root", $lumi, 1500);
#addWeights("data/merged/RSGravToGG_kMpl-02_M-3000.root", $lumi, 3000);
#addWeights("data/merged/RSGravToGG_kMpl-02_M-5000.root", $lumi, 5000);
