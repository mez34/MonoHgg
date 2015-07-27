#! /bin/sh
#run extractJSONS.py for all samples in file

#signal MC
python extractJSONS.py -i DMHggDataset.json	-o Higgs_scalar_nohdecay_gg_1000GeV 	-d lists_sig1
python extractJSONS.py -i DMHggDataset.json	-o Higgs_scalar_nohdecay_gg_100GeV 	-d lists_sig1
python extractJSONS.py -i DMHggDataset.json	-o Higgs_scalar_nohdecay_gg_10GeV 	-d lists_sig1
python extractJSONS.py -i DMHggDataset.json	-o Higgs_scalar_nohdecay_gg_1GeV 	-d lists_sig1



#MC
python extractJSONS.py -i datasets.json -o GJet_Pt-20to40	-d lists_50ns
python extractJSONS.py -i datasets.json -o GJet_Pt-40toInf 	-d lists_50ns
python extractJSONS.py -i datasets.json -o QCD_Pt-30to40  	-d lists_50ns
python extractJSONS.py -i datasets.json -o QCD_Pt-30toInf 	-d lists_50ns
python extractJSONS.py -i datasets.json -o QCD_Pt-40toInf 	-d lists_50ns
python extractJSONS.py -i datasets.json -o ttHJetToGG_M120 	-d lists_50ns
python extractJSONS.py -i datasets.json -o ttHJetToGG_M125 	-d lists_50ns
python extractJSONS.py -i datasets.json -o ttHJetToGG_M130 	-d lists_50ns
#Data
python extractJSONS.py -i datasets.json -o SinglePhoton 	-d lists_50ns
python extractJSONS.py -i datasets.json -o DoubleEG	 	-d lists_50ns
python extractJSONS.py -i datasets.json -o EGamma	 	-d lists_50ns
