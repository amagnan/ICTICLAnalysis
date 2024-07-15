import FWCore.ParameterSet.Config as cms

areaCheck = cms.EDAnalyzer('CellAreaChecker',
                     hgcalRecHitsEE     = cms.InputTag("HGCalRecHit"        , "HGCEERecHits"),
                     hgcalRecHitsFH     = cms.InputTag("HGCalRecHit"        , "HGCHEFRecHits"),
                     hgcalRecHitsBH     = cms.InputTag("HGCalRecHit"        , "HGCHEBRecHits"),
                     Debug = cms.int32(0)
                           )
