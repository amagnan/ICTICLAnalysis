import FWCore.ParameterSet.Config as cms
from Configuration.ProcessModifiers.convertHGCalDigisSim_cff import convertHGCalDigisSim

from Configuration.Eras.Era_Phase2_cff import Phase2

# options to customise the production of the ntuples
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing('analysis')

options.register('fillTriplets',
                 0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Enable fill triplets")

options.register('debug',
                 0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Enable debug printout")

#options.register('inputTracksters',
#                 ('ticlTrackstersEM1','ticlTrackstersEM3'),
#                 VarParsing.VarParsing.multiplicity.list,
#                 VarParsing.VarParsing.varType.string,
#                 "List of Input tracksters collections")


options.parseArguments()

process = cms.Process("ticlTree")

from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
process = cms.Process('PROD',Phase2C9)
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
#process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.MessageLogger.cerr.FwkReport.reportEvery = 50

process.TFileService = cms.Service("TFileService",
        fileName = cms.string(options.outputFile),
        closeFileFast = cms.untracked.bool(True)
        )

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles)
)

process.load('SimCalorimetry.HGCalSimProducers.hgcHitAssociation_cfi')
process.load('SimCalorimetry.HGCalAssociatorProducers.LCToSCAssociation_cfi')
process.load('RecoLocalCalo.HGCalRecProducers.hgcalRecHitMapProducer_cfi')

process.sim_task = cms.Task(
    process.hgcalRecHitMapProducer,
    process.scAssocByEnergyScoreProducer,
    process.layerClusterSimClusterAssociation)

process.load("ICTICLAnalysis.TiCLTreeProducer.TiCLTreeProducer_cfi")
process.ticlTree.FillTripletsInfo = cms.int32(options.fillTriplets)
process.ticlTree.Debug = cms.int32(options.debug)
process.ticlTree.trksterVec          = cms.VInputTag(
    #cms.InputTag("ticlSimTracksters"  ,  "", "RECO"),
    #cms.InputTag("ticlTrackstersDummy1"    , ""                 , "TICL" ),
    #cms.InputTag("ticlTrackstersDummy2"    , ""                 , "TICL" ),
    cms.InputTag("ticlTrackstersDummy3"    , ""                 , "TICL" ),
    #cms.InputTag("ticlTrackstersEM1"       , ""                 , "TICL" ),
    #cms.InputTag("ticlTrackstersEM2"       , ""                 , "TICL" ),
    cms.InputTag("ticlTrackstersEM3"       , ""                 , "TICL" ),
#    cms.InputTag("ticlTrackstersEM3a"      , ""                 , "TICL" ),
#    cms.InputTag("ticlTrackstersEM3b"      , ""                 , "TICL" ),
#    cms.InputTag("ticlTrackstersEM3c"      , ""                 , "TICL" ),
    #cms.InputTag("ticlTrackstersHAD1"      , ""                 , "TICL" ),
    #cms.InputTag("ticlTrackstersHAD2"      , ""                 , "TICL" ),
    cms.InputTag("ticlTrackstersHAD3"      , ""                 , "TICL" ),
    #cms.InputTag("ticlTrackstersTRK1"      , ""                 , "TICL" ),
    #cms.InputTag("ticlTrackstersTRK2"      , ""                 , "TICL" ),
    cms.InputTag("ticlTrackstersTRK3"      , ""                 , "TICL" ),
    cms.InputTag("ticlTrackstersTrkEM"     , ""                 , "TICL" ),
    cms.InputTag("ticlTrackstersEM"        , ""                 , "TICL" ),
    cms.InputTag("ticlTrackstersTrk"       , ""                 , "TICL" ),
    cms.InputTag("ticlTrackstersHAD"       , ""                 , "TICL" ),
    cms.InputTag("ticlSimTracksters"       , ""                 , "TICL" )
)
process.ticlTree.iterTypeVec = cms.vstring(
    #"Dummy1","Dummy2",
    "Dummy3",
    #"EM1","EM2",
    "EM3",
 #   "EM3a","EM3b","EM3c",
    #"HAD1","HAD2",
    "HAD3",
    #"TRK1","TRK2",
    "TRK3",
    "TrkEM","EM","Trk","HAD",
    "Sim"
)

#process.pid.trksterVec = cms.VInputTag(options.inputTracksters)
process.ticl_seq = cms.Sequence(
    process.sim_task
)

process.p = cms.Path(process.ticl_seq*process.ticlTree)
