import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.PythonUtilities.LumiList as LumiList  
import FWCore.ParameterSet.Types as CfgTypes  

isMC = True;
#should actually not need to change the bools below
is25ns = True;
is2015D = True;
is2015DFromChiara = True;

process = cms.Process("diPhoAna")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag import GlobalTag

if ((isMC==False and is2015D)):
    process.GlobalTag = GlobalTag(process.GlobalTag, '74X_dataRun2_Prompt_v2', '')
    print "74X_dataRun2_Prompt_v2"
elif (isMC):
    process.GlobalTag = GlobalTag(process.GlobalTag, 'MCRUN2_74_V9', '')
    print "MCRUN2_74_V9"


#process.GlobalTag.globaltag = 'MCRUN2_74_V9A' 		#50ns
#process.GlobalTag.globaltag = 'POSTLS170_V5::All' 	#Phys14

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( -1 ) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(

#'/store/group/phys_higgs/cmshgg/musella/flashgg/ExoPhys14_v2/diphotonsPhys14V2/RSGravToGG_kMpl001_M_5000_Tune4C_13TeV_pythia8/ExoPhys14_v2-diphotonsPhys14V2-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150128_133931/0000/myOutputFile_1.root'
                                                              )
if (isMC==False and is2015DFromChiara):
    print "applying 2015D json"                                
    process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())  
    JSONfile = '/afs/cern.ch/user/c/crovelli/public/json2015/doubleEG/processedAndGolden_2015D_oct25.json'
    myLumis = LumiList.LumiList(filename = JSONfile).getCMSSWString().split(',')  
    process.source.lumisToProcess.extend(myLumis)  

process.load("flashgg/MicroAOD/flashggPhotons_cfi")
process.load("flashgg/MicroAOD/flashggDiPhotons_cfi")

process.TFileService = cms.Service("TFileService",fileName = cms.string("OUTPUT"))

process.diPhoAna = cms.EDAnalyzer('NewDiPhoAnalyzer',
                                  VertexTag = cms.untracked.InputTag('offlineSlimmedPrimaryVertices'),
				  METTag=cms.untracked.InputTag('slimmedMETs'),
                                  genPhotonExtraTag = cms.InputTag("flashggGenPhotonsExtra"),   
                                  reducedBarrelRecHitCollection = cms.InputTag('reducedEgamma','reducedEBRecHits'),
                                  reducedEndcapRecHitCollection = cms.InputTag('reducedEgamma','reducedEERecHits'),
                                  DiPhotonTag = cms.untracked.InputTag('flashggDiPhotons'),
                                  PileUpTag = cms.untracked.InputTag('slimmedAddPileupInfo'),
                                  generatorInfo = cms.InputTag("generator"),
				  bits	        = cms.InputTag('TriggerResults::HLT'),
                                  dopureweight = PU, 
                                  sampleIndex  = SI,
                                  puWFileName  = weights,
                                  xsec         = XS,
                                  kfac         = KF,
                                  sumDataset   = SDS
                                  )

process.p = cms.Path(process.diPhoAna)

