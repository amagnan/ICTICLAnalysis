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

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9

process = cms.Process("ticlTree",Phase2C17I13M9)

process.load('Configuration.Geometry.GeometryExtended2026D98Reco_cff')

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
#process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T25', '')

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

process.load("ICTICLAnalysis.TiCLTreeProducer.CellAreaChecker_cfi")
process.areaCheck.Debug = cms.int32(options.debug)

#process.ticlTree.trksterVec          = cms.VInputTag(
    #cms.InputTag("ticlTrackstersTrk"       , ""                 , "TICL" ),
    #cms.InputTag("ticlTrackstersHAD"       , ""                 , "TICL" ),
    #cms.InputTag("ticlSimTracksters"       , ""                 , "TICL" )
#    cms.InputTag("ticlSimTracksters", "fromCPs", "TICL"),
    # ORIG
#    cms.InputTag("ticlTrackstersTrkEM"     , ""                 , "RECO" ),
#    cms.InputTag("ticlTrackstersEM"        , ""                 , "RECO" ),
    #RERECO
#    cms.InputTag("ticlTrackstersEM3"       , ""                 , "TICL" ),
#    cms.InputTag("ticlTrackstersCLUE3D3"   , ""                 , "TICL" ),
    
#
#)
#process.ticlTree.iterTypeVec = cms.vstring(
 #    "SimFromCP",
 #   "TrkEM","EM",
 #   "EM3","CLUE3D3" #@@ TICL !!
#)

#process.pid.trksterVec = cms.VInputTag(options.inputTracksters)
process.ticl_seq = cms.Sequence(
    process.sim_task
)

#process.p = cms.Path(process.ticl_seq*process.ticlTree*process.areaCheck)
process.p = cms.Path(process.ticl_seq*process.ticlTree)
#process.p = cms.Path(process.ticl_seq*process.areaCheck)
