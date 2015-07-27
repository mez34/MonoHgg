#! /bin/sh
#run extractJSONS.py for all samples in file


#signal MC
python extractFilesAndWeight.py -i list_sig1/Higgs_scalar_nohdecay_gg_1000GeV.json	-o Higgs_scalar_nohdecay_gg_1000GeV	-d lists_sig1
python extractFilesAndWeight.py -i list_sig1/Higgs_scalar_nohdecay_gg_100GeV.json	-o Higgs_scalar_nohdecay_gg_100GeV	-d lists_sig1
python extractFilesAndWeight.py -i list_sig1/Higgs_scalar_nohdecay_gg_10GeV.json	-o Higgs_scalar_nohdecay_gg_10GeV	-d lists_sig1
python extractFilesAndWeight.py -i list_sig1/Higgs_scalar_nohdecay_gg_1GeV.json		-o Higgs_scalar_nohdecay_gg_1GeV	-d lists_sig1

#MC
python extractFilesAndWeight.py -i lists_50ns/GJet_Pt-20to40.json	-o GJet_Pt-20to40	-d lists_50ns
python extractFilesAndWeight.py -i lists_50ns/GJet_Pt-40toInf.json	-o GJet_Pt-40toInf 	-d lists_50ns
python extractFilesAndWeight.py -i lists_50ns/QCD_Pt-30to40.json	-o QCD_Pt-30to40  	-d lists_50ns
python extractFilesAndWeight.py -i lists_50ns/QCD_Pt-30toInf.json	-o QCD_Pt-30toInf 	-d lists_50ns
python extractFilesAndWeight.py -i lists_50ns/QCD_Pt-40toInf.json	-o QCD_Pt-40toInf 	-d lists_50ns
python extractFilesAndWeight.py -i lists_50ns/ttHJetToGG_M120.json	-o ttHJetToGG_M120 	-d lists_50ns
python extractFilesAndWeight.py -i lists_50ns/ttHJetToGG_M125.json	-o ttHJetToGG_M125 	-d lists_50ns
python extractFilesAndWeight.py -i lists_50ns/ttHJetToGG_M130.json	-o ttHJetToGG_M130 	-d lists_50ns
#Data?
python extractFilesAndWeight.py -i lists_50ns/SinglePhoton.json		-o SinglePhoton 	-d lists_50ns
python extractFilesAndWeight.py -i lists_50ns/DoubleEG.json		-o DoubleEG	 	-d lists_50ns
python extractFilesAndWeight.py -i lists_50ns/EGamma.json 		-o EGamma	 	-d lists_50ns
