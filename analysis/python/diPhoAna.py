import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("diPhoAna")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'POSTLS170_V5::All'
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
#process.GlobalTag.globaltag = 'POSTLS170_V5'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( 50000 ) )

process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring(

        # Phys14, GG+jets, AnV1 
        #"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GGJets_M-500To1000_Pt-50_13TeV-sherpa/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141301/0000/diphotonsMicroAOD_1.root"

	# DM Hgg, Livia
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_1.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_2.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_3.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_4.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_5.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_6.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_7.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_8.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_9.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_10.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_11.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_12.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_13.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_14.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_15.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_16.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_17.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_18.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_19.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_20.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_21.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_22.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_23.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_24.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_25.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_26.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_27.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_28.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_29.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_30.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_31.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_32.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_33.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_34.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_35.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_36.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_37.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_38.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_39.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_40.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_41.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_42.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_43.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_44.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_45.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_46.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_47.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_48.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_49.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_50.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_51.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_52.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_53.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_54.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_55.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_56.root",
	#"/store/group/phys_higgs/soffi/flashgg/testMonoHLivia/Phys14MicroAODV3-55-gc1f8d91/Higgs_scalar/testMonoHLivia-Phys14MicroAODV3-55-gc1f8d91-v0-soffi-Higgs_scalar_nohdecay_gg_100GeV_13TeV_MINIAODSIM_v6-7d492cb64f2cdaff326f939f96e45c96/150626_174708/0000/myMicroAODOutputFile_57.root",

	# WZHtoGG, Livia
	#"/store/group/phys_higgs/soffi/MonoX/MonoH/MicroAOD/test/MicroAOD_WZHToGG_M-125_13TeV.root"


	#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_21.root"
      
        #Phys14, GGH from Seth
        #"/store/group/phys_higgs/soffi/MonoX/MonoH/MicroAOD/test/MicroAOD_GluGluToHToGG_M-125_13TeV.root"

        # Phys14, RS 3TeV, AnV1 
        #"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/RSGravToGG_kMpl-01_M-3000_Tune4C_13TeV-pythia8/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU30bx50_PHYS14_25_V1-v1/150330_141554/0000/diphotonsMicroAOD_2.root"

	# Background samples:
	#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_21.root" 


	#QCD samples
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_1.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_2.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_3.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_4.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_5.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_6.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_7.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_8.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_9.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_10.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_11.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_12.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_13.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_14.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_15.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_16.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_17.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_18.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_19.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_20.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_21.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_22.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_23.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_24.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_25.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_26.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_27.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_28.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_29.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_30.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_31.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_32.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_33.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_34.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_35.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_36.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_37.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_38.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_39.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_40.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_41.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_42.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_43.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_44.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_45.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_46.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_47.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_48.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_49.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_50.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_51.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_52.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_53.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_54.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_55.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_56.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_57.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_58.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_59.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_60.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_61.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_62.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_63.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_64.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_65.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_66.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_67.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_68.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_69.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_70.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_71.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_72.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_73.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_74.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_75.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_76.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_77.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_78.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_79.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_80.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_81.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_82.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_83.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_84.root",
"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/QCD_HT-100To250_13TeV-madgraph/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141436/0000/diphotonsMicroAOD_85.root",

   #GJets Sample
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_1.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_2.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_3.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_4.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_5.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_6.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_7.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_8.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_9.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_10.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_11.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_12.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_13.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_14.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_15.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_16.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_17.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_18.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_19.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_20.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_21.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_22.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_23.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_24.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_25.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_26.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_27.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_28.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_29.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_3.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_30.root"
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_31.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_32.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_33.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_34.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_35.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_36.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_37.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_38.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_39.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_40.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_41.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_42.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_43.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_44.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_45.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_46.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_47.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_48.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_49.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_50.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_51.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_52.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_53.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_54.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_55.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_56.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_57.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_58.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_59.root",
#"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14ANv1/diphotonsPhys14AnV1/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/ExoPhys14ANv1-diphotonsPhys14AnV1-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150330_141300/0000/diphotonsMicroAOD_60.root",


        )
                            )

process.load("flashgg/MicroAOD/flashggPhotons_cfi")
process.load("flashgg/MicroAOD/flashggDiPhotons_cfi")

process.TFileService = cms.Service("TFileService",fileName = cms.string("diPhotons_QCD_nosel.root"))

process.diPhoAna = cms.EDAnalyzer('DiPhoAnalyzer',
                                  VertexTag = cms.untracked.InputTag('offlineSlimmedPrimaryVertices'),
                                  METTag=cms.untracked.InputTag('slimmedMETs'),
                                  genPhotonExtraTag = cms.InputTag("flashggGenPhotonsExtra"),    
                                  DiPhotonTag = cms.untracked.InputTag('flashggDiPhotons'),
                                  reducedBarrelRecHitCollection = cms.InputTag('reducedEgamma','reducedEBRecHits'),
                                  reducedEndcapRecHitCollection = cms.InputTag('reducedEgamma','reducedEERecHits'),
                                  PileupTag = cms.untracked.InputTag('addPileupInfo'),
                                  generatorInfo = cms.InputTag("generator"),
                                  dopureweight = cms.untracked.int32(0),
                                  sampleIndex  = cms.untracked.int32(5),
                                  puWFileName  = cms.string('xxx'),   # chiara  
                                  xsec         = cms.untracked.double(28730000.),
                                  kfac         = cms.untracked.double(1.),
                                  sumDataset   = cms.untracked.double(49253.)
                                  #sumDataset   = cms.untracked.double(49972.0)
                                  )

process.p = cms.Path(process.diPhoAna)

