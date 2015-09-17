#! /bin/sh
#run extractJSONS.py for all samples in file

python extractJSONS.py -i datssets_746_dipho.json -o DiPhoton				-d lists_50ns/MC

python extractJSONS.py -i datasets_746_WZH.json -o WplusH				-d lists_50ns/MC
python extractJSONS.py -i datasets_746_WZH.json -o WminusH				-d lists_50ns/MC
python extractJSONS.py -i datasets_746_WZH.json -o ZH 					-d lists_50ns/MC


python extractJSONS.py -i datasets_746_sig.json -o Higgs_scalar_nohdecay_gg_1000GeV	-d lists_50ns/MC
python extractJSONS.py -i datasets_746_sig.json -o Higgs_scalar_nohdecay_gg_100GeV	-d lists_50ns/MC
python extractJSONS.py -i datasets_746_sig.json -o Higgs_scalar_nohdecay_gg_10GeV	-d lists_50ns/MC
python extractJSONS.py -i datasets_746_sig.json -o Higgs_scalar_nohdecay_gg_1GeV	-d lists_50ns/MC
##signal MC
#python extractJSONS.py -i datasets_746sig.json	-o Higgs_scalar_nohdecay_gg_1000GeV 	-d lists_sig1
#python extractJSONS.py -i datasets_746sig.json	-o Higgs_scalar_nohdecay_gg_100GeV 	-d lists_sig1
#python extractJSONS.py -i datasets_746sig.json	-o Higgs_scalar_nohdecay_gg_10GeV 	-d lists_sig1
#python extractJSONS.py -i datasets_746sig.json	-o Higgs_scalar_nohdecay_gg_1GeV 	-d lists_sig1
#
##MC
python extractJSONS.py -i datasets.json 	-o GJet_Pt-20to40			-d lists_50ns/MC
python extractJSONS.py -i datasets.json 	-o GJet_Pt-40toInf 			-d lists_50ns/MC
python extractJSONS.py -i datasets.json 	-o QCD_Pt-30to40  			-d lists_50ns/MC
python extractJSONS.py -i datasets.json 	-o QCD_Pt-30toInf 			-d lists_50ns/MC
python extractJSONS.py -i datasets.json 	-o QCD_Pt-40toInf 			-d lists_50ns/MC
#
python extractJSONS.py -i datasets_746.json	-o GluGluHToGG				-d lists_50ns/MC
#
#python extractJSONS.py -i datasets.json 	-o ttHJetToGG_M120 			-d lists_50ns/MC
#python extractJSONS.py -i datasets.json 	-o ttHJetToGG_M125 			-d lists_50ns/MC
#python extractJSONS.py -i datasets.json 	-o ttHJetToGG_M130 			-d lists_50ns/MC
##Data
#python extractJSONS.py -i datasets.json 	-o SinglePhoton 			-d lists_50ns/Data
#python extractJSONS.py -i datasets_746dEG.json	-o DoubleEG	 			-d lists_50ns/Data
#python extractJSONS.py -i datasets.json 	-o EGamma	 			-d lists_50ns/Data
