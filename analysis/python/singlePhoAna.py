import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("singlePhoAna")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'POSTLS170_V5::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( -1 ) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 100 )
process.MessageLogger.cerr.threshold = 'ERROR' # can't get suppressWarning to work: disable all warnings for now

process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring(
        # phys14, kMpl-01, M-3000
        "/store/group/phys_higgs/cmshgg/flashgg/ExoPhys14/diphotons-phys14-v1/RSGravToGG_kMpl-01_M-3000_Tune4C_13TeV-pythia8/ExoPhys14-diphotons-phys14-v1-v0-Phys14DR-PU30bx50_PHYS14_25_V1-v1/141202_121354/0000/myOutputFile_1.root",
        "/store/group/phys_higgs/cmshgg/flashgg/ExoPhys14/diphotons-phys14-v1/RSGravToGG_kMpl-01_M-3000_Tune4C_13TeV-pythia8/ExoPhys14-diphotons-phys14-v1-v0-Phys14DR-PU30bx50_PHYS14_25_V1-v1/141202_121354/0000/myOutputFile_2.root",
        "/store/group/phys_higgs/cmshgg/flashgg/ExoPhys14/diphotons-phys14-v1/RSGravToGG_kMpl-01_M-3000_Tune4C_13TeV-pythia8/ExoPhys14-diphotons-phys14-v1-v0-Phys14DR-PU30bx50_PHYS14_25_V1-v1/141202_121354/0000/myOutputFile_3.root"
        ## 5TeV                                         
        #"root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/flashgg/ExoPhys14/alphaV1-96-g5e4dc54/RSGravToGG_kMpl-02_M-5000_Tune4C_13TeV-pythia8/ExoPhys14-alphaV1-96-g5e4dc54-v1-Phys14DR-PU20bx25_PHYS14_25_V1-v1/141112_020912/0000/myOutputFile_1.root",
        #"root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/flashgg/ExoPhys14/alphaV1-96-g5e4dc54/RSGravToGG_kMpl-02_M-5000_Tune4C_13TeV-pythia8/ExoPhys14-alphaV1-96-g5e4dc54-v1-Phys14DR-PU20bx25_PHYS14_25_V1-v1/141112_020912/0000/myOutputFile_2.root",
        #"root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/flashgg/ExoPhys14/alphaV1-96-g5e4dc54/RSGravToGG_kMpl-02_M-5000_Tune4C_13TeV-pythia8/ExoPhys14-alphaV1-96-g5e4dc54-v1-Phys14DR-PU20bx25_PHYS14_25_V1-v1/141112_020912/0000/myOutputFile_3.root"
        )
                            )

process.load("flashgg/MicroAODProducers/flashggPhotons_cfi")

process.eventCount = cms.EDProducer("EventCountProducer")

from flashgg.MicroAODProducers.flashggMicroAODOutputCommands_cff import microAODDefaultOutputCommand

process.TFileService = cms.Service("TFileService",fileName = cms.string("singlePhotonTree.root"))
process.singlePhoAna = cms.EDAnalyzer('SinglePhoAnalyzer',
                                      reducedBarrelRecHitCollection = cms.InputTag('reducedEgamma','reducedEBRecHits'),
                                      reducedEndcapRecHitCollection = cms.InputTag('reducedEgamma','reducedEERecHits')
                                      )

process.p = cms.Path(process.singlePhoAna)

