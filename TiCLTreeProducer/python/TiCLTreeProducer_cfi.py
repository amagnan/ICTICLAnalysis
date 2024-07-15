import FWCore.ParameterSet.Config as cms

ticlTree = cms.EDAnalyzer('TiCLTreeProducer',
                     caloParticles      = cms.InputTag("mix"                , "MergedCaloTruth"),
                     genParticles       = cms.InputTag("genParticles"                , "", "HLT"),
                     simTracks          = cms.InputTag("g4SimHits"                , "", "SIM"),
                     simVertices        = cms.InputTag("g4SimHits"                , "", "SIM"),
                     simClusters        = cms.InputTag("mix"                , "MergedCaloTruth"),
                     hgcalRecHitsEE     = cms.InputTag("HGCalRecHit"        , "HGCEERecHits"),
                     hgcalRecHitsFH     = cms.InputTag("HGCalRecHit"        , "HGCHEFRecHits"),
                     hgcalRecHitsBH     = cms.InputTag("HGCalRecHit"        , "HGCHEBRecHits"),
                     layerClusterSimClusterAssociator = cms.untracked.InputTag("layerClusterSimClusterAssociation"),
                     hgcalLayerClusters = cms.InputTag("hgcalMergeLayerClusters" , ""                 , "RECO"),
                     layerClusterTime   = cms.InputTag("hgcalMergeLayerClusters" , "timeLayerCluster" , "RECO"),
                     hgcalLayerClustersNEW = cms.InputTag("hgcalMergeLayerClusters" , ""                 , "TICL"),
                     layerClusterTimeNEW   = cms.InputTag("hgcalMergeLayerClusters" , "timeLayerCluster" , "TICL"),
                     trksterVec          = cms.VInputTag( 
                         cms.InputTag("ticlSimTracksters"  ,  "fromCPs", "RECO"),
                         cms.InputTag("ticlTrackstersMerge",  "", "RECO"),
                         cms.InputTag("ticlTrackstersCLUE3DHigh"  ,  "", "RECO"),
                         cms.InputTag("ticlTrackstersCLUE3DHigh"  ,  "", "TICL"),
                     ),
                         #cms.InputTag("ticlTrackstersDummy1"       , ""                 , "TICL" ),
                         #cms.InputTag("ticlTrackstersDummy2"       , ""                 , "TICL" ),
                         #cms.InputTag("ticlTrackstersDummy3"       , ""                 , "TICL" ),
                         #cms.InputTag("ticlTrackstersEM1"       , ""                 , "TICL" ),
                         #cms.InputTag("ticlTrackstersEM2"       , ""                 , "TICL" ),
                         #cms.InputTag("ticlTrackstersEM3"       , ""                 , "TICL" ),
                         #cms.InputTag("ticlTrackstersEM3relax"       , ""                 , "TICL" ),
                         #cms.InputTag("ticlTrackstersEMDef"       , ""                 , "TICL" )),
                     iterTypeVec = cms.vstring(
                         "SimFromCP", "Merge","CLUE3DHigh","CLUE3DNEWLC"
                     ),
                         #"Dummy1","Dummy2","Dummy3",
                         #"EM1","EM2","EM3",
                         #"EM3relax","EMDef",
                     FillTripletsInfo = cms.int32(0),
                     Debug = cms.int32(0)
                 )
