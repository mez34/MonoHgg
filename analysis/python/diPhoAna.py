import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.PythonUtilities.LumiList as LumiList  
import FWCore.ParameterSet.Types as CfgTypes  

isMC = True;
#shouldn't need to change the bools below:
is25ns = True;
is2015D =  True;
is2015DFromChiara = True;

process = cms.Process("diPhoAna")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag.globaltag = 'POSTLS170_V5::All'     # Phys14 samples
#process.GlobalTag.globaltag = 'MCRUN2_74_V9A'         # 50ns

if ((isMC==False and is2015D)):
    process.GlobalTag = GlobalTag(process.GlobalTag, '74X_dataRun2_Prompt_v2', '')
    print "74X_dataRun2_Prompt_v2"
elif (isMC):
    process.GlobalTag = GlobalTag(process.GlobalTag, 'MCRUN2_74_V9', '')
    print "MCRUN2_74_V9"

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( 5000 ) )

process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring(

	# Spring15 DATA
#"/store/group/phys_higgs/cmshgg/mdonega/flashgg/RunIISpring15-50ns/Spring15BetaV2/DoubleEG/RunIISpring15-50ns-Spring15BetaV2-v0-Run2015B-PromptReco-v1/150716_154839/0000/myMicroAODOutputFile_1.root ",
#"/store/group/phys_higgs/cmshgg/mdonega/flashgg/RunIISpring15-50ns/Spring15BetaV2/DoubleEG/RunIISpring15-50ns-Spring15BetaV2-v0-Run2015B-PromptReco-v1/150716_154839/0000/myMicroAODOutputFile_10.root ",
#"/store/group/phys_higgs/cmshgg/mdonega/flashgg/RunIISpring15-50ns/Spring15BetaV2/DoubleEG/RunIISpring15-50ns-Spring15BetaV2-v0-Run2015B-PromptReco-v1/150716_154839/0000/myMicroAODOutputFile_100.root "


#Spring15 MC
#"/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISpring15-50ns/Spring15BetaV2/GluGluHToGG_M-125_13TeV_powheg_pythia8/RunIISpring15-50ns-Spring15BetaV2-v0-RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/150722_193245/0000/myMicroAODOutputFile_1.root", 
#"/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISpring15-50ns/Spring15BetaV2/GluGluHToGG_M-125_13TeV_powheg_pythia8/RunIISpring15-50ns-Spring15BetaV2-v0-RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/150722_193245/0000/myMicroAODOutputFile_10.root", 
#"/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISpring15-50ns/Spring15BetaV2/GluGluHToGG_M-125_13TeV_powheg_pythia8/RunIISpring15-50ns-Spring15BetaV2-v0-RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/150722_193245/0000/myMicroAODOutputFile_2.root"

	#"/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISpring15-50ns/Spring15BetaV4/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15-50ns-Spring15BetaV4-v0-RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/150813_120650/0000/myMicroAODOutputFile_106.root"
	#"/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISpring15-50ns/Spring15BetaV4/GluGluHToGG_M-125_13TeV_powheg_pythia8/RunIISpring15-50ns-Spring15BetaV4-v0-RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/150812_181349/0000/myMicroAODOutputFile_1.root"
	#"/store/user/mzientek/RunIISpring15-50ns/Higgs_scalar/Higgs_scalar_nohdecay_gg_1000GeV_13TeV_RunIISpring15-50ns-Spring15BetaV1_MetaV3-v0-RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/150818_162723/0000/myMicroAODOutputFile_1.root"
	#"file:myMicroAODOutputFile.root"
	#"/store/group/phys_higgs/cmshgg/musella/flashgg/EXOSpring15_7412_v2/diphotons_7412_v1/SingleElectron/EXOSpring15_7412_v2-diphotons_7412_v1-v2-Run2015D-PromptReco-v3/151006_034715/0000/diphotonsMicroAOD_1.root"
	#"/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISpring15-ReMiniAOD-BetaV7-25ns/Spring15BetaV7/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/RunIISpring15-ReMiniAOD-BetaV7-25ns-Spring15BetaV7-v0-RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151021_152953/0000/myMicroAODOutputFile_1.root",
	#"/store/group/phys_higgs/cmshgg/musella/flashgg/EXOSpring15_7415_v2/diphotons_7415_v2/DoubleEG/EXOSpring15_7415_v2-diphotons_7415_v2-v0-Run2015D-05Oct2015-v1/151019_005512/0000/diphotonsMicroAOD_1.root"
	#"/store/group/phys_higgs/cmshgg/mdonega/flashgg/RunIISpring15-50ns/Spring15BetaV2/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/RunIISpring15-50ns-Spring15BetaV2-v0-RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/150716_155016/0000/myMicroAODOutputFile_1.root" 
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia2/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia2-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_1000GeV_13TeV_MINIAODSIM_v11-7d492cb64f2cdaff326f939f96e45c96/150724_112944/0000/myMicroAODOutputFile_1.root"
	#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_1.root",
	#"/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISpring15-ReMiniAOD-BetaV7-25ns/Spring15BetaV7/DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa/RunIISpring15-ReMiniAOD-BetaV7-25ns-Spring15BetaV7-v0-RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151021_151551/0000/myMicroAODOutputFile_10.root" 
	#"/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISpring15-ReMiniAOD-BetaV7-25ns/Spring15BetaV7/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/RunIISpring15-ReMiniAOD-BetaV7-25ns-Spring15BetaV7-v0-RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151021_151907/0000/myMicroAODOutputFile_1.root" 
	#"/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISpring15-ReMiniAOD-BetaV7-25ns/Spring15BetaV7/QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/RunIISpring15-ReMiniAOD-BetaV7-25ns-Spring15BetaV7-v0-RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151021_152445/0000/myMicroAODOutputFile_1.root" 
	"/store/user/mzientek/MonoHgg_2HDM_MZP600_A0300_13TeV/mzientek_test_MonoH_2HDM_MZP600_A0300_RunIISpring15-25ns-Spring15BetaV7/151104_142438/0000/myMicroAODOutputFile_1.root"
        )
                            )
if (isMC==False and is2015DFromChiara):
    print "applying 2015D json"                                
    process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())  
    JSONfile = '/afs/cern.ch/user/c/crovelli/public/json2015/doubleEG/processedAndGolden_2015D_oct25.json'
    myLumis = LumiList.LumiList(filename = JSONfile).getCMSSWString().split(',')  
    process.source.lumisToProcess.extend(myLumis)                              
    print myLumis 

process.load("flashgg/MicroAOD/flashggPhotons_cfi")
process.load("flashgg/MicroAOD/flashggDiPhotons_cfi")

process.TFileService = cms.Service("TFileService",fileName = cms.string("diPhotons.root"))

process.diPhoAna = cms.EDAnalyzer('NewDiPhoAnalyzer',
                                  VertexTag = cms.untracked.InputTag('offlineSlimmedPrimaryVertices'),
				  METTag=cms.untracked.InputTag('slimmedMETs'),
                                  genPhotonExtraTag = cms.InputTag("flashggGenPhotonsExtra"),    
                                  DiPhotonTag = cms.untracked.InputTag('flashggDiPhotons'),
                                  #reducedBarrelRecHitCollection = cms.InputTag('reducedEgamma','reducedEBRecHits'),
                                  #reducedEndcapRecHitCollection = cms.InputTag('reducedEgamma','reducedEERecHits'),
                                  PileUpTag = cms.untracked.InputTag('slimmedAddPileupInfo'),
                                  generatorInfo = cms.InputTag("generator"),
                                  dopureweight = cms.untracked.int32(0),
                                  bits         = cms.InputTag('TriggerResults::HLT'),
                                  #sampleIndex  = cms.untracked.int32(101),   
				  sampleIndex  = cms.untracked.int32(13),
                                  puWFileName  = cms.string('/afs/cern.ch/user/m/mzientek/public/pileupWeights___processedAndGolden_2015D_oct25.root'),   # chiara  
                                  xsec         = cms.untracked.double(1), #pb
                                  kfac         = cms.untracked.double(1.),
                                  sumDataset   = cms.untracked.double(100.0)   # chiara
                                  )

process.p = cms.Path(process.diPhoAna)
