// system include files
#include <memory>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "Geometry/CaloTopology/interface/HGCalTopology.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticleFwd.h"

#include "SimDataFormats/CaloAnalysis/interface/SimClusterFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"

#include "SimDataFormats/Associations/interface/LayerClusterToSimClusterAssociator.h"
#include "SimDataFormats/Associations/interface/LayerClusterToSimClusterAssociatorBaseImpl.h"

#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "RecoLocalCalo/HGCalRecProducers/interface/HGCalClusteringAlgoBase.h"

#include "DataFormats/ParticleFlowReco/interface/HGCalMultiCluster.h"
#include "DataFormats/Common/interface/ValueMap.h"

// from HGC Validator code
#include "Validation/HGCalValidation/interface/HGVHistoProducerAlgo.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "DataFormats/HGCalReco/interface/Trackster.h"

//ROOT includes
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <string>
#include <vector>


#include "ICTICLAnalysis/TiCLTreeProducer/interface/CommonDataFormats.h"


#include "ICTICLAnalysis/TiCLTreeProducer/interface/LightTree.h"

#include "TSystem.h"
#include "TVector3.h"
#include "TLorentzVector.h"


bool sortLCsByEnergyAndLayer(const layercluster& a, const layercluster& b) {
  return (a.energy_ > b.energy_) || ((a.energy_ == b.energy_) && (a.layer_ < b.layer_));
}
bool sortLCsByLayerAndEnergy(const layercluster& a, const layercluster& b) {
  return (a.layer_ < b.layer_) || ((a.layer_ == b.layer_) && (a.energy_ > b.energy_));
}

double cosTheta(const ROOT::Math::XYZVector & AC,
		const ROOT::Math::XYZVector & AB){
  
  return sqrt( AC.Dot(AB)*AC.Dot(AB) / (AC.Mag2()*AB.Mag2()) );
  
};

class TiCLTreeProducer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit TiCLTreeProducer(const edm::ParameterSet&);
  ~TiCLTreeProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  virtual void fillHitMap(std::map<DetId, const HGCRecHit*>& hitMap,
                          const HGCRecHitCollection& rechitsEE,
                          const HGCRecHitCollection& rechitsFH,
                          const HGCRecHitCollection& rechitsBH) const;
  std::vector<unsigned> getClosestTrackstersToCPByIndex(const caloparticle& cp,
							const reco::CaloClusterCollection& lcs,
							const std::vector<ticl::Trackster>& trkster,
							float maxDrTrksterCP,
							const std::map<DetId, const HGCRecHit*>& hitMap);

  double getTrksterEnFromCP(const ticl::Trackster& trkster,
                           const reco::CaloClusterCollection& lcs,
                           const caloparticle& cp);
  void fillLCvector(std::vector<layercluster> & aLCvec,
		    const reco::BasicCluster& aLC,
		    const unsigned & aIdx,
		    const int & tsMult,
		    const std::map<DetId, const HGCRecHit*>& hitMap);
  void getLCs(const reco::CaloClusterCollection& lcs, 
	      std::vector<layercluster> & out,
	      const std::map<DetId, const HGCRecHit*>& hitMap);
  void getLCsFromTrkster(const ticl::Trackster& trkster, 
			 const reco::CaloClusterCollection& lcs, 
			 std::vector<layercluster> & out,
			 const std::map<DetId, const HGCRecHit*>& hitMap);

  void initialiseLCTreeVariables();

  void fillDoubletsInfo(const ticl::Trackster & thisTrackster,
			const reco::CaloClusterCollection &  layerClusters,
			LightTree & myTree
			);

  std::shared_ptr<hgcal::RecHitTools> recHitTools;
  
  // ----------member data ---------------------------
  edm::EDGetTokenT<std::vector<CaloParticle>> caloParticlesToken_;
  edm::EDGetTokenT<std::vector<SimCluster>> simClustersToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsEEToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsFHToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsBHToken_;

  std::vector<edm::EDGetTokenT<std::vector<ticl::Trackster> > > trksterTokenVec_;

  edm::InputTag associatorLayerClusterSimCluster_;
  edm::EDGetTokenT<hgcal::SimToRecoCollectionWithSimClusters> associatorMapSimToRecoToken_;
  edm::EDGetTokenT<hgcal::RecoToSimCollectionWithSimClusters> associatorMapRecoToSimToken_;  

  edm::EDGetTokenT<reco::CaloClusterCollection> hgcalLayerClustersToken_;
  edm::EDGetTokenT<edm::ValueMap<std::pair<float, float> > > hgcalLayerClusterTimeToken_;

  unsigned debug;

  const unsigned nL = 28;

  unsigned nIters;

  edm::RunNumber_t irun;
  edm::EventNumber_t ievent;
  edm::LuminosityBlockNumber_t ilumiblock;

  std::vector<LightTree> lighttreeVec;
  std::vector<std::string> iterTypeVec;

  int fillTripletsInfo;

  TTree *treeAllLC;

  int nAllLC;
  std::vector<double> all_lc_energy;
  std::vector<double> all_lc_eta;
  std::vector<double> all_lc_phi;
  std::vector<double> all_lc_seedEnergy;
  std::vector<double> all_lc_seedEta;
  std::vector<double> all_lc_seedPhi;
  std::vector<double> all_lc_x;
  std::vector<double> all_lc_y;
  std::vector<double> all_lc_z;
  std::vector<int> all_lc_layer;
  std::vector<int> all_lc_nrechits;
  std::vector<int> all_lc_mult;

};


TiCLTreeProducer::TiCLTreeProducer(const edm::ParameterSet& iConfig):
  caloParticlesToken_(consumes<std::vector<CaloParticle>>(iConfig.getParameter<edm::InputTag>("caloParticles"))),
  simClustersToken_(consumes<std::vector<SimCluster>>(iConfig.getParameter<edm::InputTag>("simClusters"))),
  hgcalRecHitsEEToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsEE"))),
  hgcalRecHitsFHToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsFH"))),
  hgcalRecHitsBHToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsBH"))),
  //dummyTrksterToken_(consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("dummyTrkster"))),
//trkTrksterToken_(consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("trkTrkster"))),
//hadTrksterToken_(consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("hadTrkster"))),
//mipTrksterToken_(consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("mipTrkster"))),
//mergeTrksterToken_(consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("mergeTrkster"))),
  associatorLayerClusterSimCluster_(iConfig.getUntrackedParameter<edm::InputTag>("layerClusterSimClusterAssociator")),
  associatorMapSimToRecoToken_(consumes<hgcal::SimToRecoCollectionWithSimClusters>(associatorLayerClusterSimCluster_)),
  associatorMapRecoToSimToken_(consumes<hgcal::RecoToSimCollectionWithSimClusters>(associatorLayerClusterSimCluster_)),
  hgcalLayerClustersToken_(consumes<reco::CaloClusterCollection>(iConfig.getParameter<edm::InputTag>("hgcalLayerClusters"))),
  hgcalLayerClusterTimeToken_(consumes<edm::ValueMap<std::pair<float, float>>>(iConfig.getParameter<edm::InputTag>("layerClusterTime"))) {

  recHitTools.reset(new hgcal::RecHitTools());
  //now do what ever initialization is needed
  usesResource("TFileService");
  edm::Service<TFileService> file;

  fillTripletsInfo = iConfig.getParameter<int>("FillTripletsInfo");
  debug = iConfig.getParameter<int>("Debug");

  std::cout << " -- Running on Tracksters: " << std::endl;
  std::vector<edm::InputTag> inputVec = iConfig.getParameter<std::vector<edm::InputTag> >("trksterVec");
  nIters = inputVec.size();
  trksterTokenVec_.reserve(nIters);
  lighttreeVec.resize(nIters);
  iterTypeVec = iConfig.getParameter<std::vector<std::string> >("iterTypeVec");
  for (unsigned iT(0); iT<nIters; ++iT){
    trksterTokenVec_.push_back(consumes<std::vector<ticl::Trackster>>(inputVec[iT]));
    std::cout << inputVec[iT] << " " << iterTypeVec[iT] << std::endl;
    lighttreeVec[iT].makeTree(file,iterTypeVec[iT],fillTripletsInfo);
  }
  treeAllLC = file->make<TTree>("treeLC", "treeAllLC");

  treeAllLC->Branch("nAllLC", &nAllLC, "nAllLC/I");
  treeAllLC->Branch("all_lc_energy", &all_lc_energy);
  treeAllLC->Branch("all_lc_eta", &all_lc_eta);
  treeAllLC->Branch("all_lc_phi", &all_lc_phi);
  treeAllLC->Branch("all_lc_seedEnergy", &all_lc_seedEnergy);
  treeAllLC->Branch("all_lc_seedEta", &all_lc_seedEta);
  treeAllLC->Branch("all_lc_seedPhi", &all_lc_seedPhi);
  treeAllLC->Branch("all_lc_layer", &all_lc_layer);
  treeAllLC->Branch("all_lc_nrechits", &all_lc_nrechits);
  treeAllLC->Branch("all_lc_mult", &all_lc_mult);
  treeAllLC->Branch("all_lc_x", &all_lc_x);
  treeAllLC->Branch("all_lc_y", &all_lc_y);
  treeAllLC->Branch("all_lc_z", &all_lc_z);
 
}

TiCLTreeProducer::~TiCLTreeProducer() {}

void TiCLTreeProducer::fillHitMap(std::map<DetId, const HGCRecHit*>& hitMap,
                     const HGCRecHitCollection& rechitsEE,
                     const HGCRecHitCollection& rechitsFH,
                     const HGCRecHitCollection& rechitsBH) const {
  hitMap.clear();
  for (const auto& hit : rechitsEE) {
    hitMap.emplace(hit.detid(), &hit);
  }

  for (const auto& hit : rechitsFH) {
    hitMap.emplace(hit.detid(), &hit);
  }

  for (const auto& hit : rechitsBH) {
    hitMap.emplace(hit.detid(), &hit);
  }
}  // end of TiCLTreeProducer::fillHitMap


std::vector<unsigned> TiCLTreeProducer::getClosestTrackstersToCPByIndex(const caloparticle& cp,
								   const reco::CaloClusterCollection& lcs,
								   const std::vector<ticl::Trackster>& trksters,
								   float maxDrTrksterCP,
								   const std::map<DetId, const HGCRecHit*>& hitMap) {
  std::vector<unsigned> closestTrksters_;
  unsigned idx = 0;
  for (auto const& t : trksters) {
    std::vector<layercluster> lcsFromTrkster;
    getLCsFromTrkster(t, lcs, lcsFromTrkster,hitMap);
    if (lcsFromTrkster.size()==0) continue;
    std::sort(lcsFromTrkster.begin(), lcsFromTrkster.end(), sortLCsByEnergyAndLayer);
      
    float dr = reco::deltaR(cp.eta_, cp.phi_, lcsFromTrkster[0].eta_, lcsFromTrkster[0].phi_);
    if (dr < maxDrTrksterCP) {
      closestTrksters_.push_back(idx);
    }
    idx++;
  }
  return closestTrksters_;
}

double TiCLTreeProducer::getTrksterEnFromCP(const ticl::Trackster& trkster,
					    const reco::CaloClusterCollection& lcs,
					    const caloparticle& cp) {
  //std::cout << " -- IN getTrksterEnFromCP: trackster E = " << trkster.raw_energy() << std::endl;
  double enFromRecHits_ = 0.;

  // get the indices associated to the LCs of this trackster
  const std::vector<unsigned int>& lcIdxs_ = trkster.vertices();
  //if (debug) std::cout << " --- number of LCs " << lcIdxs_.size() << std::endl;
  // loop over these idxs
  for (unsigned int ilc = 0; ilc < lcIdxs_.size(); ++ilc) {
    // get the lc and the corresponding rechits to this lc
    const reco::BasicCluster& lc = lcs[lcIdxs_[ilc]];
    auto const& hf = lc.hitsAndFractions();
    // loop over the rechits of this specific layer cluster
    for (unsigned int j = 0; j < hf.size(); j++) {
      const DetId detid_ = hf[j].first;
      std::map<DetId,double>::const_iterator lFound = cp.rechitsmap_.find(detid_);
      if (lFound != cp.rechitsmap_.end()) {
        enFromRecHits_ += lFound->second;//already multiplied by fraction
      }
    }  // end of looping over the rechits
  }    // end of looping over the idxs

  //std::cout << " OUT E = " << enFromRecHits_ << std::endl;
  return enFromRecHits_;
}  // end of getTrksterEnFromCP

void TiCLTreeProducer::getLCs(const reco::CaloClusterCollection& lcs,
			      std::vector<layercluster> & layerclusters,
			      const std::map<DetId, const HGCRecHit*>& hitMap) {
  
  for (unsigned int ilc = 0; ilc < lcs.size(); ++ilc) {
    const reco::BasicCluster& lc = lcs.at(ilc);
    fillLCvector(layerclusters,lc,ilc,-2,hitMap);
  }

}  // end of getLCs


void TiCLTreeProducer::getLCsFromTrkster(const ticl::Trackster& trkster,
					 const reco::CaloClusterCollection& lcs,
					 std::vector<layercluster> & layerclusters,
					 const std::map<DetId, const HGCRecHit*>& hitMap) {
  const std::vector<unsigned int>& lcIdxs = trkster.vertices();

  for (unsigned int ilc = 0; ilc < lcIdxs.size(); ++ilc) {
    const reco::BasicCluster& lc = lcs.at(lcIdxs[ilc]);
    fillLCvector(layerclusters,lc,lcIdxs[ilc],trkster.vertex_multiplicity(ilc),hitMap);
  }

}  // end of getLCsFromTrkster


void TiCLTreeProducer::fillLCvector(std::vector<layercluster> & aLCvec,
				    const reco::BasicCluster& aLC,
				    const unsigned & aIdx,
				    const int & tsMult,
				    const std::map<DetId, const HGCRecHit*>& hitMap){
  auto const& hf = aLC.hitsAndFractions();
  
  int layer_ = 0;
  for (unsigned int j = 0; j < hf.size(); j++) {
    const DetId detid_ = hf[j].first;
    layer_ = recHitTools->getLayerWithOffset(detid_);
    break;
  }
  //for (unsigned int j = 0; j < hf.size(); j++) {
  //std::cout << aLC.printHitAndFraction(j) << std::endl;
  //}
  
  aLCvec.push_back(layercluster());
  auto& layercluster_ = aLCvec.back();
  layercluster_.energy_ = aLC.energy();///(double)trkster.vertex_multiplicity(ilc);
  layercluster_.eta_ = aLC.eta();
  layercluster_.phi_ = aLC.phi();

  //get seed information
  auto const itcheck = hitMap.find(aLC.seed());
  if (itcheck != hitMap.end()) {
    const HGCRecHit* hit = itcheck->second;
    layercluster_.seedEnergy_ = hit->energy();///(double)trkster.vertex_multiplicity(ilc);
  }
  //get position wrt (0,0,0)
  layercluster_.seedEta_ = recHitTools->getEta(aLC.seed());
  layercluster_.seedPhi_ = recHitTools->getPhi(aLC.seed());


  if (tsMult>-2){
    layercluster_.algo_ = aLC.algo();
    layercluster_.tsMult_ = tsMult;
  }
  
  layercluster_.x_ = aLC.position().x();
  layercluster_.y_ = aLC.position().y();
  layercluster_.z_ = aLC.position().z();
  layercluster_.nrechits_ = aLC.hitsAndFractions().size();
  layercluster_.layer_ = abs(layer_);
  layercluster_.idxTracksterLC_ = aIdx;
}

void TiCLTreeProducer::initialiseLCTreeVariables(){

  nAllLC = 0;
  all_lc_energy.clear();
  all_lc_eta.clear();
  all_lc_phi.clear();
  all_lc_seedEnergy.clear();
  all_lc_seedEta.clear();
  all_lc_seedPhi.clear();
  all_lc_layer.clear();
  all_lc_nrechits.clear();
  all_lc_mult.clear();
  all_lc_x.clear();
  all_lc_y.clear();
  all_lc_z.clear();

};

void TiCLTreeProducer::fillDoubletsInfo(const ticl::Trackster & thisTrackster,
			   const reco::CaloClusterCollection &  layerClusters,
			   LightTree & myTree
			   ){

  std::vector< std::vector<Triplet> > tripletVec;
  std::vector<Triplet> dummy;
  dummy.reserve(10);
  tripletVec.resize(nL,dummy);
  //loop on doublets: consider it is the outerDoublet.
  for (const auto &edge : thisTrackster.edges()) {//loop on edges
    auto & ic = layerClusters[edge[0]];//B
    auto & oc = layerClusters[edge[1]];//C
    auto const & cl_in = ic.hitsAndFractions()[0].first;
    auto const & cl_out = oc.hitsAndFractions()[0].first;
    auto const layer_in = recHitTools->getLayerWithOffset(cl_in);
    auto const layer_out = recHitTools->getLayerWithOffset(cl_out);
    
    Triplet lTriplet;
    lTriplet.initialise();
    lTriplet.layerB_ = layer_in;
    lTriplet.layerC_ = layer_out;
    lTriplet.eB_ = ic.energy();
    lTriplet.eC_ = oc.energy();
    lTriplet.etaB_ = ic.eta();

    // Alpha angles
    const auto & outer_outer_pos = oc.position();
    const auto & outer_inner_pos = ic.position();
    //const auto & seed = thisTrackster.seedIndex();
    auto seedGlobalPos = math::XYZPoint(0,0,0);
    
    lTriplet.cosAlphaOuter_ = (outer_inner_pos-seedGlobalPos).Dot(outer_outer_pos - outer_inner_pos)/
      sqrt((outer_inner_pos-seedGlobalPos).Mag2()*(outer_outer_pos-outer_inner_pos).Mag2());
    
    // To complete the triplet, another inner loop
    // is therefore needed.
    std::vector<std::array<unsigned int, 2>> innerDoublets;
    std::vector<std::array<unsigned int, 2>> outerDoublets;
    for ( const auto & otherEdge : thisTrackster.edges()) {
      if (otherEdge[1] == edge[0]) {
	innerDoublets.push_back(otherEdge);
      }
      if (edge[1] == otherEdge[0]) {
	outerDoublets.push_back(otherEdge);
      }
    }
    
    lTriplet.outer_in_links_ = innerDoublets.size();
    lTriplet.outer_out_links_ = outerDoublets.size();

    if (innerDoublets.size()==0) {
      tripletVec[layer_in-1].push_back(lTriplet);
    }

    for (const auto & inner : innerDoublets) {
      const auto & inner_ic = layerClusters[inner[0]];
      const auto & inner_inner_pos = inner_ic.position();

      lTriplet.cosBeta_ = (outer_inner_pos-inner_inner_pos).Dot(outer_outer_pos-inner_inner_pos)/
	sqrt((outer_inner_pos-inner_inner_pos).Mag2()*(outer_outer_pos-inner_inner_pos).Mag2());
      
      lTriplet.eA_ = inner_ic.energy();
      lTriplet.layerA_ = recHitTools->getLayerWithOffset(inner_ic.hitsAndFractions()[0].first);

      lTriplet.cosAlphaInner_ = (inner_inner_pos-seedGlobalPos).Dot(outer_inner_pos - inner_inner_pos)/
	sqrt((inner_inner_pos-seedGlobalPos).Mag2()*(outer_inner_pos-inner_inner_pos).Mag2());

      //get the number of links also for this inner doublet
      std::vector<std::array<unsigned int, 2>> IinnerDoublets;
      std::vector<std::array<unsigned int, 2>> IouterDoublets;
      for ( const auto & otherEdge : thisTrackster.edges()) {
	if (otherEdge[1] == edge[0]) {
	  IinnerDoublets.push_back(otherEdge);
	}
	if (edge[1] == otherEdge[0]) {
	  IouterDoublets.push_back(otherEdge);
	}
      }
    
      lTriplet.inner_in_links_ = IinnerDoublets.size();
      lTriplet.inner_out_links_ = IouterDoublets.size();

      tripletVec[layer_in-1].push_back(lTriplet);

   }//loop on inner doublets
  }//loop on edges

  myTree.fillTriplets(tripletVec);
  
}//fillDoubletsInfo



void TiCLTreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  edm::ESHandle<CaloGeometry> geom;
  iSetup.get<CaloGeometryRecord>().get(geom);
  recHitTools->setGeometry(*geom);

  unsigned nLayers = recHitTools->lastLayerBH();

  edm::Handle<HGCRecHitCollection> recHitHandleEE;
  iEvent.getByToken(hgcalRecHitsEEToken_, recHitHandleEE);

  edm::Handle<HGCRecHitCollection> recHitHandleFH;
  iEvent.getByToken(hgcalRecHitsFHToken_, recHitHandleFH);

  edm::Handle<HGCRecHitCollection> recHitHandleBH;
  iEvent.getByToken(hgcalRecHitsBHToken_, recHitHandleBH);


  //edm::Handle<std::vector<ticl::Trackster>> dummyTrksterHandle;
  //iEvent.getByToken(dummyTrksterToken_, dummyTrksterHandle);
  //const std::vector<ticl::Trackster>& dummyTrksters = *dummyTrksterHandle;

  std::vector < const std::vector<ticl::Trackster>* > trkstersVec;
  for (unsigned iT(0); iT<nIters; ++iT){
    edm::Handle<std::vector<ticl::Trackster> > trksterHandle;
    iEvent.getByToken(trksterTokenVec_[iT], trksterHandle);
    const std::vector<ticl::Trackster>& trkstersEle = *trksterHandle;
    trkstersVec.push_back(&trkstersEle);
  }

  edm::Handle<std::vector<CaloParticle>> CaloParticles;
  iEvent.getByToken(caloParticlesToken_, CaloParticles);
  const CaloParticleCollection& cps = *CaloParticles;

  edm::Handle<hgcal::SimToRecoCollectionWithSimClusters> simToRecoCollHandle;
  iEvent.getByToken(associatorMapSimToRecoToken_, simToRecoCollHandle);
  const hgcal::SimToRecoCollectionWithSimClusters& simToRecoColl = *simToRecoCollHandle;

  edm::Handle<hgcal::RecoToSimCollectionWithSimClusters> recoToSimCollHandle;
  iEvent.getByToken(associatorMapRecoToSimToken_, recoToSimCollHandle);
  const hgcal::RecoToSimCollectionWithSimClusters& recoToSimColl = *recoToSimCollHandle;
  auto recSimColl = *recoToSimCollHandle; 

  edm::Handle<reco::CaloClusterCollection> layerClusterHandle;
  iEvent.getByToken(hgcalLayerClustersToken_, layerClusterHandle);
  const reco::CaloClusterCollection& lcs = *layerClusterHandle;

  std::map<DetId, const HGCRecHit*> hitMap;
  fillHitMap(hitMap, *recHitHandleEE, *recHitHandleFH, *recHitHandleBH);

  //edm::Handle<edm::ValueMap<std::pair<float, float>>> lcTimeHandle;
  //iEvent.getByToken(hgcalLayerClusterTimeToken_, lcTimeHandle);
  //const auto& lcTime = *lcTimeHandle;

  // init vars
  //recHitTools->getEventSetup(iSetup);

  if (debug) std::cout << " - Processing event " << iEvent.id().event() << std::endl;

  initialiseLCTreeVariables();

  // get CaloParticles
  std::vector<caloparticle> caloparticles;
  std::vector<simcluster> simclust;
  int idx = -1;
  for (const auto& it_cp : cps) {
    const CaloParticle& cp = ((it_cp));

    if ((cp.eventId().event() != 0) || (cp.eventId().bunchCrossing() != 0)) {
      continue;
    }

    // allow CP only within HGCAL volume
    if ((abs(cp.eta()) < 1.2) || (abs(cp.eta()) > 3.5)) {
      continue;
    }
    idx++;

    caloparticle tmpcp_;
    tmpcp_.idx_ = idx;
    tmpcp_.pdgid_ = cp.pdgId();
    tmpcp_.energy_ = cp.energy();
    tmpcp_.pt_ = cp.pt();
    tmpcp_.eta_ = cp.eta();
    tmpcp_.phi_ = cp.phi();
    //std::cout << tmpcp_.print() << std::endl;

    // get the simclusters
    const SimClusterRefVector& simclusters = cp.simClusters();
    tmpcp_.nSC_ = simclusters.size();
    for (const auto& it_simc : simclusters) {
      const SimCluster& simc = (*(it_simc));
      const auto& sc_haf = simc.hits_and_fractions();

      simcluster tmpsc_;
      tmpsc_.idx_ = idx;
      tmpsc_.pdgid_ = simc.pdgId();
      tmpsc_.energy_ = simc.energy();
      tmpsc_.pt_ = simc.pt();
      tmpsc_.eta_ = simc.eta();
      tmpsc_.phi_ = simc.phi();
      simclust.push_back(tmpsc_);
      // get the rechits
      for (const auto& it_sc_haf : sc_haf) {
        DetId detid_ = (it_sc_haf.first);

        // need to map RecHits to the SimCluster
        // SimHits are not stored
        auto const itcheck = hitMap.find(detid_);
        // we need this check because some DetIDs assigned to CaloParticle do not always have
        // a RecHit -- due to thresholds or whatever
        if (itcheck != hitMap.end()) {
          const HGCRecHit* hit = itcheck->second;
	  std::pair<std::map<DetId,double>::iterator,bool> lInsert = tmpcp_.rechitsmap_.insert(std::pair<DetId,double>(detid_,(it_sc_haf.second) * hit->energy()));
	  if (!lInsert.second) lInsert.first->second += (it_sc_haf.second) * hit->energy();
        }  //  end of if(itcheck != hitMap.end())
      }    // end of looping over the rechits
    }      // end of looping over the sim clusters
    
    caloparticles.push_back(tmpcp_);
    
  }  // end of looping over the calo particles

  // get the relevant trackster collection
  // LG: For now always use the MergedTrackster
  //   int cp_pdgid_ = 0; if (caloparticles.size()>0) { cp_pdgid_ = caloparticles.at(0).pdgid_; }
  //   std::vector<ticl::Trackster> tracksters = getTracksterCollection(cp_pdgid_,emMCs, mipMCs, hadMCs, mergedMCs);

  // loop over the caloparticles and then find the closest trackster to it

  // keep tracksters [and the corresponding LC] that
  // have at least some ammount of energy from the CP
  //   std::vector<trackster> trksterCollection; trksterCollection.clear();

  //double trackstersEnFromCP = 0.9;
  irun = iEvent.id().run();
  ievent = iEvent.id().event();
  ilumiblock = iEvent.id().luminosityBlock();
  

  

  all_lc_energy.clear();
  all_lc_eta.clear();
  all_lc_phi.clear();
  all_lc_seedEnergy.clear();
  all_lc_seedEta.clear();
  all_lc_seedPhi.clear();
  all_lc_x.clear();
  all_lc_y.clear();
  all_lc_z.clear();
  all_lc_layer.clear();
  all_lc_nrechits.clear();
  all_lc_mult.clear();

  std::vector<layercluster> alllcs;
  getLCs(lcs,alllcs,hitMap);
  std::sort(alllcs.begin(), alllcs.end(), sortLCsByLayerAndEnergy);
  
  std::map<int,int>lMapLC;
  nAllLC = alllcs.size();
  all_lc_energy.reserve(nAllLC);
  all_lc_eta.reserve(nAllLC);
  all_lc_phi.reserve(nAllLC);
  all_lc_seedEnergy.reserve(nAllLC);
  all_lc_seedEta.reserve(nAllLC);
  all_lc_seedPhi.reserve(nAllLC);
  all_lc_x.reserve(nAllLC);
  all_lc_y.reserve(nAllLC);
  all_lc_z.reserve(nAllLC);
  all_lc_layer.reserve(nAllLC);
  all_lc_nrechits.reserve(nAllLC);
  all_lc_mult.reserve(nLayers);
  
  for (auto const& lc : alllcs) {
    all_lc_energy.push_back(lc.energy_);
    all_lc_eta.push_back(lc.eta_);
    all_lc_phi.push_back(lc.phi_);
    all_lc_seedEnergy.push_back(lc.seedEnergy_);
    all_lc_seedEta.push_back(lc.seedEta_);
    all_lc_seedPhi.push_back(lc.seedPhi_);
    all_lc_x.push_back(lc.x_);
    all_lc_y.push_back(lc.y_);
    all_lc_z.push_back(lc.z_);
    all_lc_layer.push_back(lc.layer_);
    std::pair<std::map<int,int>::iterator,bool> isInserted =  lMapLC.insert(std::pair<int,int>(lc.layer_,1));
    if (!isInserted.second) isInserted.first->second += 1;
    
    all_lc_nrechits.push_back(lc.nrechits_);
  }
  
  for (unsigned iL(0); iL<nLayers;++iL){
    std::map<int,int>::iterator lEle = lMapLC.find(iL+1);
    if (lEle !=lMapLC.end()) all_lc_mult.push_back(lEle->second);
    else all_lc_mult.push_back(0);
  }
  
  treeAllLC->Fill();

  for (unsigned iT(0); iT<nIters; ++iT){
    if (debug) std::cout << " -- Processing iter " << iterTypeVec[iT] << std::endl;

    //fill this here - duplicated across tracksters trees...
    lighttreeVec[iT].initialiseTreeVariables((size_t)irun,(int)ievent,(size_t)ilumiblock); 
    lighttreeVec[iT].fillCPinfo(caloparticles);
    lighttreeVec[iT].fillSCinfo(simclust);
    
    auto const& tracksters = *(trkstersVec[iT]);//dummyTrksters;//mergeTrksters;  // if we need a different trackster collection this should go in the loop
    //auto const& tracksters = dummyTrksters;//mergeTrksters;  // if we need a different trackster collection this should go in the loop
    
    // find the tracksters within some DR from the CP

    float maxDrTracksterCP = 0.3;
    int itrksterMin_ = -1;
    for (unsigned int icp = 0; icp < caloparticles.size(); ++icp) {
      
      //AM@TODO
      //Add also CPidx and SCidx of closest CP/SC to trackster...

      std::vector<unsigned> closestTracksters =
	getClosestTrackstersToCPByIndex(caloparticles[icp], lcs, tracksters, maxDrTracksterCP,hitMap);

      if (debug) std::cout << " --- CaloParticle " << icp << " E=" << caloparticles[icp].energy_<< " nb closest tracksters = " << closestTracksters.size() << std::endl;

      // for those tracksters closest to the CP calculate the energy that is associated to the CP
      // and select the one with the closest Energy to the CP
      double trksterCPEnDiffMin_ = std::numeric_limits<double>::max();
      for (unsigned itrkster : closestTracksters) {
	double trksterEnFromCP_ = getTrksterEnFromCP(tracksters[itrkster], lcs, caloparticles[icp]);
	double trksterCPEnDiff = abs(caloparticles[icp].energy_ - trksterEnFromCP_) / (caloparticles[icp].energy_);
	if (debug>1) std::cout << " ---- trackster idx " << itrkster << " with E " << tracksters[itrkster].raw_energy() << " EfromCP = " << trksterEnFromCP_ << std::endl;
	if (trksterCPEnDiff < trksterCPEnDiffMin_) {
	  trksterCPEnDiffMin_ = trksterCPEnDiff;
	  //if (trksterCPEnDiff < trackstersEnFromCP) {
          itrksterMin_ = itrkster;
          lighttreeVec[iT].fillCPEfraction(trksterCPEnDiff,icp);
	  //}
	  //else std::cout << " Don't record missingE..." << std::endl;
	}
      }
      if (debug>0) std::cout << " --- Chose index " << itrksterMin_ << " !" << std::endl;

    }//loop over CP
      
      
    std::vector<std::vector<layercluster> > lcsFromTrksters;
    lcsFromTrksters.reserve(tracksters.size());
    for (unsigned its = 0; its < tracksters.size(); ++its) {
      std::vector<layercluster> lcsFromTrkster;
      getLCsFromTrkster(tracksters[its], lcs, lcsFromTrkster,hitMap);
      std::sort(lcsFromTrkster.begin(), lcsFromTrkster.end(), sortLCsByEnergyAndLayer);
      lcsFromTrksters.push_back(lcsFromTrkster);
    }
    
    lighttreeVec[iT].fillTSinfo(tracksters,nLayers,lcsFromTrksters,layerClusterHandle,recSimColl);
      
    if (itrksterMin_>=0 && fillTripletsInfo) fillDoubletsInfo(tracksters[itrksterMin_],lcs,lighttreeVec[iT]);
    
    lighttreeVec[iT].fillOutputTree();
  }//loop on iterations
    
}//analyze

// ------------ method called once each job just before starting event loop  ------------
void TiCLTreeProducer::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void TiCLTreeProducer::endJob() {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TiCLTreeProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TiCLTreeProducer);
