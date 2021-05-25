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
                     hgcalLayerClusters = cms.InputTag("hgcalLayerClusters" , ""                 , "RECO"),
                     layerClusterTime   = cms.InputTag("hgcalLayerClusters" , "timeLayerCluster" , "RECO"),
                     trksterVec          = cms.VInputTag( 
                         cms.InputTag("ticlSimTracksters"  ,  "", "RECO"),
                         cms.InputTag("ticlTrackstersEM"   ,  "", "RECO"),
                         cms.InputTag("ticlTrackstersHAD"  ,  "", "RECO"),
                         cms.InputTag("ticlTrackstersMerge",  "", "RECO"),
                         cms.InputTag("ticlTrackstersTrk"  ,  "", "RECO"),
                         cms.InputTag("ticlTrackstersTrkEM",  "", "RECO")
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
                         "Sim","EM","HAD", "Merge","Trk","TrkEM"
                     ),
                         #"Dummy1","Dummy2","Dummy3",
                         #"EM1","EM2","EM3",
                         #"EM3relax","EMDef",
                     FillTripletsInfo = cms.int32(0),
                     Debug = cms.int32(0)
                 )
