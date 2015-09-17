#! /bin/sh 
# this scripts creates a merged root file in the self-created merged

mkdir -p data/50ns/

hadd data/50ns/DMHtoGG_M1000.root   ../../output/job_201594_17733/Higgs_scalar_nohdecay_gg_1000GeV/Higgs_scalar_nohdecay_gg_1000GeV_*.root  
hadd data/50ns/DMHtoGG_M100.root    ../../output/job_201595_125340/Higgs_scalar_nohdecay_gg_100GeV/Higgs_scalar_nohdecay_gg_100GeV_*.root
hadd data/50ns/DMHtoGG_M10.root     ../../output/job_201595_125345/Higgs_scalar_nohdecay_gg_10GeV/Higgs_scalar_nohdecay_gg_10GeV_*.root
hadd data/50ns/DMHtoGG_M1.root      ../../output/job_201595_125349/Higgs_scalar_nohdecay_gg_1GeV/Higgs_scalar_nohdecay_gg_1GeV_*.root

hadd data/50ns/GJet_Pt-20to40.root  ../../output/job_201595_125355/GJet_Pt-20to40/GJet_Pt-20to40*.root
hadd data/50ns/GJet_Pt-40toInf.root ../../output/job_201595_125454/GJet_Pt-40toInf/GJet_Pt-40toInf*.root

hadd data/50ns/QCD_Pt-30to40.root   ../../output/job_201594_171656/QCD_Pt-30to40/QCD_Pt-30to40*.root 
hadd data/50ns/QCD_Pt-30toInf.root  ../../output/job_201594_17176/QCD_Pt-30toInf/QCD_Pt-30toInf*.root
hadd data/50ns/QCD_Pt-40toInf.root  ../../output/job_201594_171744/QCD_Pt-40toInf/QCD_Pt-40toInf*.root

hadd data/50ns/GluGluHToGG.root     ../../output/job_201595_125557/GluGluHToGG/GluGluHToGG*.root

hadd data/50ns/DiPhoton.root	    ../../output/job_2015917_12130/DiPhoton/DiPhoton*.root

hadd data/50ns/ZH.root		    ../../output/job_2015917_112043/ZH/ZH*.root
#hadd data/50ns/ZH.root		    ../../output/job_201595_12522/ZH/ZH*.root
hadd data/50ns/WplusH.root	    ../../output/job_201595_125252/WplusH/WplusH*.root
hadd data/50ns/WminusH.root	    ../../output/job_201595_12528/WminusH/WminusH*.root

hadd data/50ns/DoubleEG.root	    ../../output/job_2015910_12489/DoubleEG/DoubleEG*.root
#hadd data/50ns/DoubleEG.root	    ../../output/job_201594_171326/DoubleEG/DoubleEG*.root


# 50ns sample without triggers
#hadd data/50ns/DMHtoGG_M1000.root   ../../output/job_2015819_123249/Higgs_scalar_nohdecay_gg_1000GeV/Higgs_scalar_nohdecay_gg_1000GeV_*.root
#hadd data/50ns/DMHtoGG_M100.root    ../../output/job_2015821_111846/Higgs_scalar_nohdecay_gg_100GeV/Higgs_scalar_nohdecay_gg_100GeV_*.root
#hadd data/50ns/DMHtoGG_M10.root     ../../output/job_2015821_11192/Higgs_scalar_nohdecay_gg_10GeV/Higgs_scalar_nohdecay_gg_10GeV_*.root
#hadd data/50ns/DMHtoGG_M1.root      ../../output/job_2015823_172439/Higgs_scalar_nohdecay_gg_1GeV/Higgs_scalar_nohdecay_gg_1GeV_*.root
#
#hadd data/50ns/GJet_Pt-20to40.root  ../../output/job_2015821_101936/GJet_Pt-20to40/GJet_Pt-20to40*.root
#hadd data/50ns/GJet_Pt-40toInf.root ../../output/job_2015821_101950/GJet_Pt-40toInf/GJet_Pt-40toInf*.root
#
#hadd data/50ns/QCD_Pt-30to40.root   ../../output/job_2015821_102018/QCD_Pt-30to40/QCD_Pt-30to40*.root
#hadd data/50ns/QCD_Pt-30toInf.root  ../../output/job_2015821_102027/QCD_Pt-30toInf/QCD_Pt-30toInf*.root
#hadd data/50ns/QCD_Pt-40toInf.root  ../../output/job_2015821_102044/QCD_Pt-40toInf/QCD_Pt-40toInf*.root
#
#hadd data/50ns/GluGluHToGG.root     ../../output/job_2015821_102239/GluGluHToGG/GluGluHToGG*.root
#
#hadd data/50ns/DoubleEG.root	    ../../output/job_2015825_95650/DoubleEG/DoubleEG*.root


#hadd data/50ns/GJet_Pt-20to40.root  ../../output/job_2015813_12337/GJet_Pt-20to40/GJet_Pt-20to40*.root
#hadd data/50ns/GJet_Pt-40toInf.root ../../output/job_2015813_123336/GJet_Pt-40toInf/GJet_Pt-40toInf*.root
##
#hadd data/50ns/QCD_Pt-30to40.root   ../../output/job_2015813_123240/QCD_Pt-30to40/QCD_Pt-30to40*.root
#hadd data/50ns/QCD_Pt-30toInf.root  ../../output/job_2015813_12314/QCD_Pt-30toInf/QCD_Pt-30toInf*.root
#hadd data/50ns/QCD_Pt-40toInf.root  ../../output/job_2015813_123037/QCD_Pt-40toInf/QCD_Pt-40toInf*.root
#
#hadd data/50ns/GluGluHToGG.root     ../../output/job_2015813_12309/GluGluHToGG/GluGluHToGG*.root
#
#hadd data/50ns/DMHtoGG_M1000.root   ../../output/job_2015813_122639/Higgs_scalar_nohdecay_gg_1000GeV/Higgs_scalar_nohdecay_gg_1000GeV_*.root
#hadd data/50ns/DMHtoGG_M100.root    ../../output/job_2015813_12271/Higgs_scalar_nohdecay_gg_100GeV/Higgs_scalar_nohdecay_gg_100GeV_*.root
#hadd data/50ns/DMHtoGG_M10.root     ../../output/job_2015813_122720/Higgs_scalar_nohdecay_gg_10GeV/Higgs_scalar_nohdecay_gg_10GeV_*.root
#hadd data/50ns/DMHtoGG_M1.root      ../../output/job_2015813_122740/Higgs_scalar_nohdecay_gg_1GeV/Higgs_scalar_nohdecay_gg_1GeV_*.root




#../../output/job_201582_161620/Higgs_scalar_nohdecay_gg_1000GeV/Higgs_scalar_nohdecay_gg_1000GeV_*.root  
#
#hadd data/merged/GGJets_M-200To500.root    data/GGJets_M-200To500/GGJets_M-200To500_*root
#hadd data/merged/GGJets_M-500To1000.root   data/GGJets_M-500To1000/GGJets_M-500To1000_*root
#hadd data/merged/GGJets_M-1000To2000.root  data/GGJets_M-1000To2000/GGJets_M-1000To2000_*root
#hadd data/merged/GGJets_M-2000To4000.root  data/GGJets_M-2000To4000/GGJets_M-2000To4000_*root
#hadd data/merged/GGJets_M-4000To8000.root  data/GGJets_M-4000To8000/GGJets_M-4000To8000_*root
#hadd data/merged/GGJets_M-8000To13000.root data/GGJets_M-8000To13000/GGJets_M-8000To13000_*root
#
#hadd data/merged/GJets_HT-100to200.root data/GJets_HT-100to200/GJets_HT-100to200*root
#hadd data/merged/GJets_HT-200to400.root data/GJets_HT-200to400/GJets_HT-200to400*root
#hadd data/merged/GJets_HT-400to600.root data/GJets_HT-400to600/GJets_HT-400to600*root
#hadd data/merged/GJets_HT-600toInf.root data/GJets_HT-600toInf/GJets_HT-600toInf*root
##
#hadd data/merged/QCD_HT-100To250.root  data/QCD_HT-100To250/QCD_HT-100To250*root
#hadd data/merged/QCD_HT-250To500.root  data/QCD_HT-250To500/QCD_HT-250To500*root
#hadd data/merged/QCD_HT-500To1000.root data/QCD_HT-500To1000/QCD_HT-500To1000*root
#hadd data/merged/QCD_HT-1000ToInf.root data/QCD_HT-1000ToInf/QCD_HT-1000ToInf*root
#
#hadd data/merged/RSGravToGG_kMpl-001_M-750.root  data/RSGravToGG_kMpl-001_M-750/RSGravToGG_kMpl-001_M-750*root
#hadd data/merged/RSGravToGG_kMpl-001_M-1500.root data/RSGravToGG_kMpl-001_M-1500/RSGravToGG_kMpl-001_M-1500*root
#hadd data/merged/RSGravToGG_kMpl-001_M-5000.root data/RSGravToGG_kMpl-001_M-5000/RSGravToGG_kMpl-001_M-5000*root
#hadd data/merged/RSGravToGG_kMpl-01_M-1500.root data/RSGravToGG_kMpl-01_M-1500/RSGravToGG_kMpl-01_M-1500*root
#hadd data/merged/RSGravToGG_kMpl-01_M-3000.root data/RSGravToGG_kMpl-01_M-3000/RSGravToGG_kMpl-01_M-3000*root
#hadd data/merged/RSGravToGG_kMpl-02_M-1500.root data/RSGravToGG_kMpl-02_M-1500/RSGravToGG_kMpl-02_M-1500*root
#hadd data/merged/RSGravToGG_kMpl-02_M-3000.root data/RSGravToGG_kMpl-02_M-3000/RSGravToGG_kMpl-02_M-3000*root
#hadd data/merged/RSGravToGG_kMpl-02_M-5000.root data/RSGravToGG_kMpl-02_M-5000/RSGravToGG_kMpl-02_M-5000*root
