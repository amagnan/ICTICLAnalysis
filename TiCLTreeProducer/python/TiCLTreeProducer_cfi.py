import FWCore.ParameterSet.Config as cms

ticlTree = cms.EDAnalyzer('TiCLTreeProducer',
                     caloParticles      = cms.InputTag("mix"                , "MergedCaloTruth"),
                     hgcalRecHitsEE     = cms.InputTag("HGCalRecHit"        , "HGCEERecHits"),
                     hgcalRecHitsFH     = cms.InputTag("HGCalRecHit"        , "HGCHEFRecHits"),
                     hgcalRecHitsBH     = cms.InputTag("HGCalRecHit"        , "HGCHEBRecHits"),
                     hgcalLayerClusters = cms.InputTag("hgcalLayerClusters" , ""                 , "RECO"),
                     layerClusterTime   = cms.InputTag("hgcalLayerClusters" , "timeLayerCluster" , "RECO"),
                     trksterVec          = cms.VInputTag( 
                         cms.InputTag("ticlTrackstersDummy1"       , ""                 , "TICL" ),
                         cms.InputTag("ticlTrackstersDummy2"       , ""                 , "TICL" ),
                         cms.InputTag("ticlTrackstersDummy3"       , ""                 , "TICL" ),
                         cms.InputTag("ticlTrackstersEM1"       , ""                 , "TICL" ),
                         cms.InputTag("ticlTrackstersEM2"       , ""                 , "TICL" ),
                         cms.InputTag("ticlTrackstersEM3"       , ""                 , "TICL" ),
                         cms.InputTag("ticlTrackstersEM3relax"       , ""                 , "TICL" ),
                         cms.InputTag("ticlTrackstersEMDef"       , ""                 , "TICL" )),
                         iterTypeVec = cms.vstring(
                             "Dummy1","Dummy2","Dummy3",
                             "EM1","EM2","EM3",
                             "EM3relax","EMDef"),
                         FillTripletsInfo = cms.bool(False),

                 )
