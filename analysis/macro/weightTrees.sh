#! /bin/sh                                                                                                    

lumi=$1      # in pb
echo "Adding weights for " $lumi " pb-1..."

root -l -b <<EOF
.L addWeightsToTree.cc+   
addWeights("data/merged/GGJets_M-200To500.root", $lumi);
addWeights("data/merged/GGJets_M-500To1000.root", $lumi);
addWeights("data/merged/GGJets_M-1000To2000.root", $lumi);
addWeights("data/merged/GGJets_M-2000To4000.root", $lumi);
addWeights("data/merged/GGJets_M-4000To8000.root", $lumi);
addWeights("data/merged/GGJets_M-8000To13000.root", $lumi);
addWeights("data/merged/GJets_HT-100to200.root", $lumi);
addWeights("data/merged/GJets_HT-200to400.root", $lumi);
addWeights("data/merged/GJets_HT-400to600.root", $lumi);
addWeights("data/merged/GJets_HT-600toInf.root", $lumi);
addWeights("data/merged/RSGravToGG_kMpl-01_M-1500.root", $lumi);
addWeights("data/merged/RSGravToGG_kMpl-01_M-3000.root", $lumi);
.q

EOF

echo "done weighting."
