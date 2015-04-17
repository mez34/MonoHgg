#! /bin/sh 
# this scripts creates a merged root file in the self-created merged

mkdir -p data/merged

hadd data/merged/GGJets_M-200To500.root    data/GGJets_M-200To500/GGJets_M-200To500_*root
hadd data/merged/GGJets_M-500To1000.root   data/GGJets_M-500To1000/GGJets_M-500To1000_*root
hadd data/merged/GGJets_M-1000To2000.root  data/GGJets_M-1000To2000/GGJets_M-1000To2000_*root
hadd data/merged/GGJets_M-2000To4000.root  data/GGJets_M-2000To4000/GGJets_M-2000To4000_*root
hadd data/merged/GGJets_M-4000To8000.root  data/GGJets_M-4000To8000/GGJets_M-4000To8000_*root
hadd data/merged/GGJets_M-8000To13000.root data/GGJets_M-8000To13000/GGJets_M-8000To13000_*root
#
hadd data/merged/GJets_HT-100to200.root data/GJets_HT-100to200/GJets_HT-100to200*root
hadd data/merged/GJets_HT-200to400.root data/GJets_HT-200to400/GJets_HT-200to400*root
hadd data/merged/GJets_HT-400to600.root data/GJets_HT-400to600/GJets_HT-400to600*root
hadd data/merged/GJets_HT-600toInf.root data/GJets_HT-600toInf/GJets_HT-600toInf*root
#
hadd data/merged/RSGravToGG_kMpl-01_M-1500.root data/RSGravToGG_kMpl-01_M-1500/RSGravToGG_kMpl-01_M-1500*root
hadd data/merged/RSGravToGG_kMpl-01_M-3000.root data/RSGravToGG_kMpl-01_M-3000/RSGravToGG_kMpl-01_M-3000*root
