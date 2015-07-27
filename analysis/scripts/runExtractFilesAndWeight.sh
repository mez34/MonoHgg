#! /bin/sh
#run extractJSONS.py for all samples in file

#MC
python extractFilesAndWeight.py -i lists_50ns/GJet_Pt-20to40.json	-o GJet_Pt-20to40	-d lists_50ns
python extractFilesAndWeight.py -i lists_50ns/GJet_Pt-40toInf.json	-o GJet_Pt-40toInf 	-d lists_50ns
python extractFilesAndWeight.py -i lists_50ns/QCD_Pt-30to40.json	-o QCD_Pt-30to40  	-d lists_50ns
python extractFilesAndWeight.py -i lists_50ns/QCD_Pt-30toInf.json	-o QCD_Pt-30toInf 	-d lists_50ns
python extractFilesAndWeight.py -i lists_50ns/QCD_Pt-40toInf.json	-o QCD_Pt-40toInf 	-d lists_50ns
python extractFilesAndWeight.py -i lists_50ns/ttHJetToGG_M120.json	-o ttHJetToGG_M120 	-d lists_50ns
python extractFilesAndWeight.py -i lists_50ns/ttHJetToGG_M125.json	-o ttHJetToGG_M125 	-d lists_50ns
python extractFilesAndWeight.py -i lists_50ns/ttHJetToGG_M120.json	-o ttHJetToGG_M120 	-d lists_50ns
#Data?
python extractFilesAndWeight.py -i lists_50ns/SinglePhoton.json		-o SinglePhoton 	-d lists_50ns
python extractFilesAndWeight.py -i lists_50ns/DoubleEG.json		-o DoubleEG	 	-d lists_50ns
python extractFilesAndWeight.py -i lists_50ns/EGamma.json 		-o EGamma	 	-d lists_50ns
