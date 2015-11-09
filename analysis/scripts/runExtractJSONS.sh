#! /bin/sh
#run extractJSONS.py for all samples in file

# All FLASHgg version Spring15BetaV7:
python extractJSONS.py -i datasets/datasets_7415_betaV7_sig.json	-o 2HDM_MZP600		-d lists_25ns_v7/MC
python extractJSONS.py -i datasets/datasets_7415_betaV7_sig.json	-o 2HDM_MZP800		-d lists_25ns_v7/MC
python extractJSONS.py -i datasets/datasets_7415_betaV7_sig.json	-o 2HDM_MZP1000		-d lists_25ns_v7/MC
python extractJSONS.py -i datasets/datasets_7415_betaV7_sig.json	-o 2HDM_MZP1200		-d lists_25ns_v7/MC
python extractJSONS.py -i datasets/datasets_7415_betaV7_sig.json	-o 2HDM_MZP1400		-d lists_25ns_v7/MC
python extractJSONS.py -i datasets/datasets_7415_betaV7_sig.json	-o 2HDM_MZP1700		-d lists_25ns_v7/MC
python extractJSONS.py -i datasets/datasets_7415_betaV7_sig.json	-o 2HDM_MZP2000		-d lists_25ns_v7/MC
python extractJSONS.py -i datasets/datasets_7415_betaV7_sig.json	-o 2HDM_MZP2500		-d lists_25ns_v7/MC
#python extractJSONS.py -i datasets/datasets_7415_betaV7.json	-o DiPhoton		-d lists_25ns_v7/MC
#python extractJSONS.py -i datasets/datasets_7415_betaV7.json	-o DYJetsToLL		-d lists_25ns_v7/MC
#python extractJSONS.py -i datasets/datasets_7415_betaV7.json	-o QCD_Pt-30to40	-d lists_25ns_v7/MC
#python extractJSONS.py -i datasets/datasets_7415_betaV7.json	-o QCD_Pt-40toInf	-d lists_25ns_v7/MC
#python extractJSONS.py -i datasets/datasets_7415_betaV7.json	-o QCD_Pt-30toInf	-d lists_25ns_v7/MC
#python extractJSONS.py -i datasets/datasets_7415_betaV7.json	-o GJet_Pt-20to40	-d lists_25ns_v7/MC
#python extractJSONS.py -i datasets/datasets_7415_betaV7.json	-o GJet_Pt-40toInf	-d lists_25ns_v7/MC
#python extractJSONS.py -i datasets/datasets_7415_betaV7.json	-o GluGluHToGG_M-125	-d lists_25ns_v7/MC
#python extractJSONS.py -i datasets/datasets_7415_betaV7.json	-o VHToGG_M125		-d lists_25ns_v7/MC
#
#python extractJSONS.py -i datasets/datasets_7415_betaV7metaV1.json	-o DoubleEG 	-d lists_25ns_v7/Data
# FLASHgg v7 data separate by runD & runC:
#python extractJSONS.py -i datasets/datasets_7415_betaV7metaV1.json	-o DoubleEG/RunIISpring15-ReMiniAOD-BetaV7-25ns-Spring15BetaV7-v0-Run2015D 	-d lists_25ns_v7/Data
#python extractJSONS.py -i datasets/datasets_7415_betaV7.json	-o DoubleEG/RunIISpring15-ReMiniAOD-BetaV7-25ns-Spring15BetaV7-v0-Run2015C	-d lists_25ns_v7/Data




## with mix of Flashgg versions (spring15 betaV2 and betaV4):

#python extractJSONS.py -i datasets_746_dipho.json -o DiPhoton				-d lists_50ns/MC
#python extractJSONS.py -i datasets_746_WZH.json -o WplusH				-d lists_50ns/MC
#python extractJSONS.py -i datasets_746_WZH.json -o WminusH				-d lists_50ns/MC
#python extractJSONS.py -i datasets_746_WZH.json -o ZH 					-d lists_50ns/MC
#python extractJSONS.py -i datasets_746_sig.json -o Higgs_scalar_nohdecay_gg_1000GeV	-d lists_50ns/MC
#python extractJSONS.py -i datasets_746_sig.json -o Higgs_scalar_nohdecay_gg_100GeV	-d lists_50ns/MC
#python extractJSONS.py -i datasets_746_sig.json -o Higgs_scalar_nohdecay_gg_10GeV	-d lists_50ns/MC
#python extractJSONS.py -i datasets_746_sig.json -o Higgs_scalar_nohdecay_gg_1GeV	-d lists_50ns/MC
#python extractJSONS.py -i datasets.json 	-o GJet_Pt-20to40			-d lists_50ns/MC
#python extractJSONS.py -i datasets.json 	-o GJet_Pt-40toInf 			-d lists_50ns/MC
#python extractJSONS.py -i datasets.json 	-o QCD_Pt-30to40  			-d lists_50ns/MC
#python extractJSONS.py -i datasets.json 	-o QCD_Pt-30toInf 			-d lists_50ns/MC
#python extractJSONS.py -i datasets.json 	-o QCD_Pt-40toInf 			-d lists_50ns/MC
#python extractJSONS.py -i datasets_746.json	-o GluGluHToGG				-d lists_50ns/MC
#python extractJSONS.py -i datasets.json 	-o ttHJetToGG_M125 			-d lists_50ns/MC
##Data
#python extractJSONS.py -i datasets.json 	-o SinglePhoton 			-d lists_50ns/Data
#python extractJSONS.py -i datasets_746dEG.json	-o DoubleEG	 			-d lists_50ns/Data
#python extractJSONS.py -i datasets.json 	-o EGamma	 			-d lists_50ns/Data
