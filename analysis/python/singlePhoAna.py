import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("singlePhoAna")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'POSTLS170_V5::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( 100 ) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 100 )
process.MessageLogger.cerr.threshold = 'ERROR' # can't get suppressWarning to work: disable all warnings for now

process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring(
        # Phys14, RSGrav, kMpl001, mG=5k    -- new production
        #"/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14_v2/diphotonsPhys14V2/RSGravToGG_kMpl001_M_5000_Tune4C_13TeV_pythia8/ExoPhys14_v2-diphotonsPhys14V2-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150128_133931/0000/myOutputFile_1.root"

        # Phys14, RSGrav, kMpl01, mG=5k   -- old production
        "/store/group/phys_higgs/cmshgg/flashgg/ExoPhys14/diphotonsPhys14V1/RSGravToGG_kMpl01_M-5000_Tune4C_13TeV-pythia8/ExoPhys14-diphotonsPhys14V1-v0-Phys14DR-PU40bx25_PHYS14_25_V1-v1/141205_002243/0000/myOutputFile_1.root"

        # Phys14, RSGrav, kMpl01, mG=3k   -- old production
        #"/store/group/phys_higgs/cmshgg/flashgg/ExoPhys14/diphotons-phys14-v1/RSGravToGG_kMpl-01_M-3000_Tune4C_13TeV-pythia8/ExoPhys14-diphotons-phys14-v1-v0-Phys14DR-PU30bx50_PHYS14_25_V1-v1/141202_121354/0000/myOutputFile_1.root",
        #"/store/group/phys_higgs/cmshgg/flashgg/ExoPhys14/diphotons-phys14-v1/RSGravToGG_kMpl-01_M-3000_Tune4C_13TeV-pythia8/ExoPhys14-diphotons-phys14-v1-v0-Phys14DR-PU30bx50_PHYS14_25_V1-v1/141202_121354/0000/myOutputFile_2.root",
        #"/store/group/phys_higgs/cmshgg/flashgg/ExoPhys14/diphotons-phys14-v1/RSGravToGG_kMpl-01_M-3000_Tune4C_13TeV-pythia8/ExoPhys14-diphotons-phys14-v1-v0-Phys14DR-PU30bx50_PHYS14_25_V1-v1/141202_121354/0000/myOutputFile_3.root" 

        # Phys14, RSGrav, kMpl01, mG=1.5k  -- old production      
        #"/store/group/phys_higgs/cmshgg/flashgg/ExoPhys14/diphotons-phys14-v1/RSGravToGG_kMpl-01_M-1500_Tune4C_13TeV-pythia8/ExoPhys14-diphotons-phys14-v1-v0-P#hys14DR-PU30bx50_PHYS14_25_V1-v1/141202_121328/0000/myOutputFile_1.root", 
        #"/store/group/phys_higgs/cmshgg/flashgg/ExoPhys14/diphotons-phys14-v1/RSGravToGG_kMpl-01_M-1500_Tune4C_13TeV-pythia8/ExoPhys14-diphotons-phys14-v1-v0-Phys14DR-PU30bx50_PHYS14_25_V1-v1/141202_121328/0000/myOutputFile_2.root", 
        #"/store/group/phys_higgs/cmshgg/flashgg/ExoPhys14/diphotons-phys14-v1/RSGravToGG_kMpl-01_M-1500_Tune4C_13TeV-pythia8/ExoPhys14-diphotons-phys14-v1-v0-Phys14DR-PU30bx50_PHYS14_25_V1-v1/141202_121328/0000/myOutputFile_3.root"
        )
                            )

process.load("flashgg/MicroAODProducers/flashggPhotons_cfi")

process.eventCount = cms.EDProducer("EventCountProducer")

from flashgg.MicroAODProducers.flashggMicroAODOutputCommands_cff import microAODDefaultOutputCommand

#process.out = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string('myOutputFile.root'),
#                               outputCommands = microAODDefaultOutputCommand
#                               )
#process.out.outputCommands.append("keep *_*_*_*")


process.TFileService = cms.Service("TFileService",fileName = cms.string("singlePhotonTree.root"))
process.singlePhoAna = cms.EDAnalyzer('SinglePhoAnalyzer',
                                      packedGenParticles = cms.InputTag('flashggGenPhotons'),
                                      reducedBarrelRecHitCollection = cms.InputTag('reducedEgamma','reducedEBRecHits'),
                                      reducedEndcapRecHitCollection = cms.InputTag('reducedEgamma','reducedEERecHits')
                                      )

process.p = cms.Path(process.singlePhoAna)

#process.e = cms.EndPath(process.out)
