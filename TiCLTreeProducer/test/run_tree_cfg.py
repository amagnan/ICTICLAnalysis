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

process.load("ICTICLAnalysis.TiCLTreeProducer.TiCLTreeProducer_cfi")
process.ticlTree.FillTripletsInfo = cms.int32(options.fillTriplets)

#process.pid.trksterVec = cms.VInputTag(options.inputTracksters)

process.p = cms.Path(process.ticlTree)
