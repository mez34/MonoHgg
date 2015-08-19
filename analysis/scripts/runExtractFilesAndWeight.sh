#! /bin/sh
#run extractJSONS.py for all samples in file

python extractFilesAndWeight.py -i lists_sig/Higgs_scalar_nohdecay_gg_1000GeV.json	-o Higgs_scalar_nohdecay_gg_1000GeV	-d lists_sig

##signal MC
#python extractFilesAndWeight.py -i lists_sig1/Higgs_scalar_nohdecay_gg_1000GeV.json	-o Higgs_scalar_nohdecay_gg_1000GeV	-d lists_sig1
#python extractFilesAndWeight.py -i lists_sig1/Higgs_scalar_nohdecay_gg_100GeV.json	-o Higgs_scalar_nohdecay_gg_100GeV	-d lists_sig1
#python extractFilesAndWeight.py -i lists_sig1/Higgs_scalar_nohdecay_gg_10GeV.json	-o Higgs_scalar_nohdecay_gg_10GeV	-d lists_sig1
#python extractFilesAndWeight.py -i lists_sig1/Higgs_scalar_nohdecay_gg_1GeV.json	-o Higgs_scalar_nohdecay_gg_1GeV	-d lists_sig1
#
##MC
#python extractFilesAndWeight.py -i lists_50ns/MC/GJet_Pt-20to40.json	-o GJet_Pt-20to40	-d lists_50ns
#python extractFilesAndWeight.py -i lists_50ns/MC/GJet_Pt-40toInf.json	-o GJet_Pt-40toInf 	-d lists_50ns
#python extractFilesAndWeight.py -i lists_50ns/MC/QCD_Pt-30to40.json	-o QCD_Pt-30to40  	-d lists_50ns
#python extractFilesAndWeight.py -i lists_50ns/MC/QCD_Pt-30toInf.json	-o QCD_Pt-30toInf 	-d lists_50ns
#python extractFilesAndWeight.py -i lists_50ns/MC/QCD_Pt-40toInf.json	-o QCD_Pt-40toInf 	-d lists_50ns
#python extractFilesAndWeight.py -i lists_50ns/MC/GluGluHToGG.json	-o GluGluHToGG  	-d lists_50ns
#python extractFilesAndWeight.py -i lists_50ns/MC/ttHJetToGG_M120.json	-o ttHJetToGG_M120 	-d lists_50ns
#python extractFilesAndWeight.py -i lists_50ns/MC/ttHJetToGG_M125.json	-o ttHJetToGG_M125 	-d lists_50ns
#python extractFilesAndWeight.py -i lists_50ns/MC/ttHJetToGG_M130.json	-o ttHJetToGG_M130 	-d lists_50ns
##Data
#python extractFilesAndWeight.py -i lists_50ns/Data/SinglePhoton.json	-o SinglePhoton 	-d lists_50ns
#python extractFilesAndWeight.py -i lists_50ns/Data/DoubleEG.json	-o DoubleEG	 	-d lists_50ns
#python extractFilesAndWeight.py -i lists_50ns/Data/EGamma.json 	-o EGamma	 	-d lists_50ns
