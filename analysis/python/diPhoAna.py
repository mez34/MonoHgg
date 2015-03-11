import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("diPhoAna")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'POSTLS170_V5::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( 200 ) )

process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring(
        # Phys14, RSGrav, kMpl001, mG=5k    -- new production
        "/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14_v2/diphotonsPhys14V2/RSGravToGG_kMpl001_M_5000_Tune4C_13TeV_pythia8/ExoPhys14_v2-diphotonsPhys14V2-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150128_133931/0000/myOutputFile_1.root"
        )
                            )

# chiara
#inputlist = FileUtils.loadListFromFile('fileList_ggToH_125_13TeV.txt')
#readFiles = cms.untracked.vstring( *inputlist)
#process.source = cms.Source("PoolSource", fileNames = readFiles)

process.load("flashgg/MicroAODProducers/flashggPhotons_cfi")
process.load("flashgg/MicroAODProducers/flashggDiPhotons_cfi")

process.TFileService = cms.Service("TFileService",fileName = cms.string("diPhotons.root"))

process.diPhoAna = cms.EDAnalyzer('DiPhoAnalyzer',
                                  VertexTag = cms.untracked.InputTag('offlineSlimmedPrimaryVertices'),
                                  packedGenParticles = cms.untracked.InputTag('flashggGenPhotons'),
                                  DiPhotonTag = cms.untracked.InputTag('flashggDiPhotons'),
                                  PileupTag = cms.untracked.InputTag('addPileupInfo'),
                                  dopureweight = cms.untracked.int32(0),
                                  sampleIndex  = cms.untracked.int32(0),
                                  puWFileName  = cms.string('xxx'),   # chiara                                                            
                                  xsec         = cms.untracked.double(1.),
                                  kfac         = cms.untracked.double(1.)
                                  )

process.p = cms.Path(process.diPhoAna)

