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

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "RecoEgamma/EgammaMCTools/interface/ElectronMCTruthFinder.h"
#include "RecoEgamma/EgammaMCTools/interface/ElectronMCTruth.h"
#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruthFinder.h"
#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruth.h"

#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticleFwd.h"

#include "SimDataFormats/CaloAnalysis/interface/SimClusterFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"

#include "SimDataFormats/Associations/interface/LayerClusterToSimClusterAssociator.h"
#include "SimDataFormats/Associations/interface/LayerClusterToSimClusterAssociatorBaseImpl.h"

#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "RecoLocalCalo/HGCalRecProducers/interface/HGCalClusteringAlgoBase.h"

#include "DataFormats/ParticleFlowReco/interface/HGCalMultiCluster.h"
#include "DataFormats/Common/interface/ValueMap.h"

// from HGC Validator code
#include "Validation/HGCalValidation/interface/HGVHistoProducerAlgo.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "DataFormats/HGCalReco/interface/Trackster.h"

//ROOT includes
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"

#include "Math/Vector4D.h"

//using namespace ROOT::VecOps;
//using simVertexRVec = ROOT::VecOps::RVec<SimVertex>;

#include <string>
#include <vector>

#include "ICTICLAnalysis/TiCLTreeProducer/interface/CommonDataFormats.h"


#include "ICTICLAnalysis/TiCLTreeProducer/interface/LightTree.h"

#include "TSystem.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "TH2F.h"


bool sortLCsByEnergyAndLayer(const layercluster& a, const layercluster& b) {
  return (a.energy_ > b.energy_) || ((a.energy_ == b.energy_) && (a.layer_ < b.layer_));
}
bool sortLCsByLayerAndEnergy(const layercluster& a, const layercluster& b) {
  return (a.type_ == b.type_) &&
    (
     (a.layer_ < b.layer_) ||
     ((a.layer_ == b.layer_) && (a.energy_ > b.energy_))
     );
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
  std::vector<unsigned> getClosestTrackstersToCPByIndex(const std::string & iterName,
							const caloparticle& cp,
							const reco::CaloClusterCollection& lcs,
							const std::vector<ticl::Trackster>& trkster,
							float maxDrTrksterCP,
							const std::map<DetId, const HGCRecHit*>& hitMap);

  double getTrksterEnFromCP(const ticl::Trackster& trkster,
                           const reco::CaloClusterCollection& lcs,
                           const caloparticle& cp);
  void fillLCvector(const std::string & iterName,
		    std::vector<layercluster> & aLCvec,
		    const reco::BasicCluster& aLC,
		    const unsigned & aIdx,
		    const double & tsMult,
		    const std::map<DetId, const HGCRecHit*>& hitMap,
		    const int aType);
  void getLCs(const reco::CaloClusterCollection& lcs, 
	      std::vector<layercluster> & out,
	      const std::map<DetId, const HGCRecHit*>& hitMap,
	      const int aType);
  void getLCsFromTrkster(const std::string & iterName,
			 const ticl::Trackster& trkster, 
			 const reco::CaloClusterCollection& lcs, 
			 std::vector<layercluster> & out,
			 const std::map<DetId, const HGCRecHit*>& hitMap);

  void initialiseLCTreeVariables();
  void initialiseRHTreeVariables(const unsigned nL);

  void fillDoubletsInfo(const ticl::Trackster & thisTrackster,
			const reco::CaloClusterCollection &  layerClusters,
			LightTree & myTree
			);

  std::shared_ptr<hgcal::RecHitTools> recHitTools;
  
  // ----------member data ---------------------------
  const edm::ESGetToken<CaloGeometry, CaloGeometryRecord> tok_geom_;

  edm::EDGetTokenT<std::vector<CaloParticle>> caloParticlesToken_;
  edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticlesToken_;
  edm::EDGetTokenT<std::vector<SimTrack>> simTracksToken_;
  edm::EDGetTokenT<std::vector<SimVertex>> simVerticesToken_;
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
  edm::EDGetTokenT<reco::CaloClusterCollection> hgcalLayerClustersNEWToken_;
  edm::EDGetTokenT<edm::ValueMap<std::pair<float, float> > > hgcalLayerClusterTimeNEWToken_;

  PhotonMCTruthFinder photonMCTruthFinder_;

  unsigned debug;

  const unsigned nLEE = 28;

  unsigned nIters;

  edm::RunNumber_t irun;
  edm::EventNumber_t ievent;
  edm::LuminosityBlockNumber_t ilumiblock;

  std::vector<LightTree> lighttreeVec;
  std::vector<std::string> iterTypeVec;

  int fillTripletsInfo;

  TTree *treeAllLC = 0;

  int nAllLC;
  std::vector<int> all_lc_type;
  std::vector<double> all_lc_energy;
  std::vector<double> all_lc_eta;
  std::vector<double> all_lc_phi;
  std::vector<double> all_lc_seedEnergy;
  std::vector<double> all_lc_seedArea;
  std::vector<double> all_lc_seedx;
  std::vector<double> all_lc_seedy;
  std::vector<double> all_lc_seedu;
  std::vector<double> all_lc_seedv;
  std::vector<double> all_lc_eminRH;
  std::vector<double> all_lc_areaminRH;
  std::vector<double> all_lc_seedEta;
  std::vector<double> all_lc_seedPhi;
  std::vector<double> all_lc_x;
  std::vector<double> all_lc_y;
  std::vector<double> all_lc_z;
  std::vector<int> all_lc_isHD;
  std::vector<int> all_lc_isSi;
  std::vector<int> all_lc_layer;
  std::vector<int> all_lc_nrechits;
  std::vector<int> all_lc_nrechitsHD;
  std::vector<double> all_lc_efracHD;
  std::vector<int> all_lc_mult;

  TTree *treeAllRH = 0;
  int rh_evt;
  std::vector<int> rh_layer;
  unsigned SNthresh[3] = {0,3,5};
  std::vector<int> rh_n[3];
  std::vector<int> rh_nHD[3];
  std::vector<int> rh_nScint[3];
  std::vector<double> rh_E[3];
  std::vector<double> rh_EHD[3];
  std::vector<double> rh_EScint[3];

  TH2F *h_rh_SoN_Scint = 0;
  TH2F *h_rh_SoN_LD = 0;
  TH2F *h_rh_SoN_HD = 0;
  
};


TiCLTreeProducer::TiCLTreeProducer(const edm::ParameterSet& iConfig):
  tok_geom_(esConsumes<CaloGeometry, CaloGeometryRecord>()),
  caloParticlesToken_(consumes<std::vector<CaloParticle>>(iConfig.getParameter<edm::InputTag>("caloParticles"))),
  genParticlesToken_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticles"))),
  simTracksToken_(consumes<std::vector<SimTrack>>(iConfig.getParameter<edm::InputTag>("simTracks"))),
  simVerticesToken_(consumes<std::vector<SimVertex>>(iConfig.getParameter<edm::InputTag>("simVertices"))),
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
  hgcalLayerClusterTimeToken_(consumes<edm::ValueMap<std::pair<float, float>>>(iConfig.getParameter<edm::InputTag>("layerClusterTime"))),
  hgcalLayerClustersNEWToken_(consumes<reco::CaloClusterCollection>(iConfig.getParameter<edm::InputTag>("hgcalLayerClustersNEW"))),
  hgcalLayerClusterTimeNEWToken_(consumes<edm::ValueMap<std::pair<float, float>>>(iConfig.getParameter<edm::InputTag>("layerClusterTimeNEW"))),
  photonMCTruthFinder_()
  {

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
    lighttreeVec[iT].makeTree(file,iterTypeVec[iT],fillTripletsInfo,debug);
  }
  treeAllLC = file->make<TTree>("treeLC", "treeAllLC");

  treeAllLC->Branch("nAllLC", &nAllLC, "nAllLC/I");
  treeAllLC->Branch("all_lc_type", &all_lc_type);
  treeAllLC->Branch("all_lc_energy", &all_lc_energy);
  treeAllLC->Branch("all_lc_eta", &all_lc_eta);
  treeAllLC->Branch("all_lc_phi", &all_lc_phi);
  treeAllLC->Branch("all_lc_seedEnergy", &all_lc_seedEnergy);
  treeAllLC->Branch("all_lc_seedArea", &all_lc_seedArea);
  treeAllLC->Branch("all_lc_seedx", &all_lc_seedx);
  treeAllLC->Branch("all_lc_seedy", &all_lc_seedy);
  treeAllLC->Branch("all_lc_seedu", &all_lc_seedu);
  treeAllLC->Branch("all_lc_seedv", &all_lc_seedv);
  treeAllLC->Branch("all_lc_eminRH", &all_lc_eminRH);
  treeAllLC->Branch("all_lc_areaminRH", &all_lc_areaminRH);
  treeAllLC->Branch("all_lc_seedEta", &all_lc_seedEta);
  treeAllLC->Branch("all_lc_seedPhi", &all_lc_seedPhi);
  treeAllLC->Branch("all_lc_isSi", &all_lc_isSi);
  treeAllLC->Branch("all_lc_isHD", &all_lc_isHD);
  treeAllLC->Branch("all_lc_layer", &all_lc_layer);
  treeAllLC->Branch("all_lc_nrechits", &all_lc_nrechits);
  treeAllLC->Branch("all_lc_nrechitsHD", &all_lc_nrechitsHD);
  treeAllLC->Branch("all_lc_efracHD", &all_lc_efracHD);
  treeAllLC->Branch("all_lc_mult", &all_lc_mult);
  treeAllLC->Branch("all_lc_x", &all_lc_x);
  treeAllLC->Branch("all_lc_y", &all_lc_y);
  treeAllLC->Branch("all_lc_z", &all_lc_z);

  treeAllRH = file->make<TTree>("treeRH", "treeAllRH");

  treeAllRH->Branch("rh_evt", &rh_evt, "rh_evt/I");
  treeAllRH->Branch("rh_layer", &rh_layer);

  for (unsigned i(0); i<3;++i){
    std::ostringstream lname;
    lname << "_SoN" << SNthresh[i];
    treeAllRH->Branch(("rh_n"+lname.str()).c_str(), &rh_n[i]);
    treeAllRH->Branch(("rh_nHD"+lname.str()).c_str(), &rh_nHD[i]);
    treeAllRH->Branch(("rh_nScint"+lname.str()).c_str(), &rh_nScint[i]);
    treeAllRH->Branch(("rh_E"+lname.str()).c_str(), &rh_E[i]);
    treeAllRH->Branch(("rh_EHD"+lname.str()).c_str(), &rh_EHD[i]);
    treeAllRH->Branch(("rh_EScint"+lname.str()).c_str(), &rh_EScint[i]);
  }

  h_rh_SoN_Scint = file->make<TH2F>("h_rh_SoN_Scint", ";layer;SoN;Scint rechits;", 20,30,50,256,0,32);
  h_rh_SoN_LD = file->make<TH2F>("h_rh_SoN_LD", ";layer;SoN;LD rechits;", 50,0,50,256,0,32);
  h_rh_SoN_HD = file->make<TH2F>("h_rh_SoN_HD", ";layer;SoN;HD rechits;", 35,0,35,256,0,32);
  
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


std::vector<unsigned> TiCLTreeProducer::getClosestTrackstersToCPByIndex(const std::string & iterName,
									const caloparticle& cp,
								   const reco::CaloClusterCollection& lcs,
								   const std::vector<ticl::Trackster>& trksters,
								   float maxDrTrksterCP,
								   const std::map<DetId, const HGCRecHit*>& hitMap) {
  std::vector<unsigned> closestTrksters_;
  unsigned idx = 0;
  for (auto const& t : trksters) {
    std::vector<layercluster> lcsFromTrkster;
    getLCsFromTrkster(iterName,t, lcs, lcsFromTrkster,hitMap);
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
			      const std::map<DetId, const HGCRecHit*>& hitMap,
			      const int aType) {
  
  for (unsigned int ilc = 0; ilc < lcs.size(); ++ilc) {
    const reco::BasicCluster& lc = lcs.at(ilc);
    fillLCvector("AllLC",layerclusters,lc,ilc,-2,hitMap,aType);
  }

}  // end of getLCs


void TiCLTreeProducer::getLCsFromTrkster(const std::string & iterName,
					 const ticl::Trackster& trkster,
					 const reco::CaloClusterCollection& lcs,
					 std::vector<layercluster> & layerclusters,
					 const std::map<DetId, const HGCRecHit*>& hitMap) {
  const std::vector<unsigned int>& lcIdxs = trkster.vertices();

  for (unsigned int ilc = 0; ilc < lcIdxs.size(); ++ilc) {
    const reco::BasicCluster& lc = lcs.at(lcIdxs[ilc]);
    int lType = (iterName.find("NEWLC")==iterName.npos) ? 0 : 1;
    fillLCvector(iterName,layerclusters,lc,lcIdxs[ilc],trkster.vertex_multiplicity(ilc),hitMap,lType);
  }

}  // end of getLCsFromTrkster


void TiCLTreeProducer::fillLCvector(const std::string & iterName,
				    std::vector<layercluster> & aLCvec,
				    const reco::BasicCluster& aLC,
				    const unsigned & aIdx,
				    const double & tsMult,
				    const std::map<DetId, const HGCRecHit*>& hitMap,
				    const int aType){
  auto const& hf = aLC.hitsAndFractions();
  
  int layer_ = 0;
  double EminRH = 1.1*aLC.energy();
  double areaMinRH = 0;
  double echeck = 0;
  int nRHHD = 0;
  double eHD = 0;

  for (unsigned int j = 0; j < hf.size(); j++) {
    const DetId detid_ = hf[j].first;
    layer_ = recHitTools->getLayerWithOffset(detid_);
    auto const hitcheck = hitMap.find(detid_);
    if (hitcheck != hitMap.end()) {
      const HGCRecHit* hit = hitcheck->second;
      double hitE = hit->energy()*hf[j].second;
      echeck += hitE;
      if (hitE < EminRH){
	EminRH = hitE;
	areaMinRH = recHitTools->getCellArea(detid_);
      }

      bool ishd = recHitTools->isSilicon(detid_)?HGCSiliconDetId(detid_).highDensity():false;
      if (ishd) {
	nRHHD++;
	eHD += hitE;
      }
    }
    //break;
  }

  if (abs(echeck-aLC.energy())>0.001) std::cout << " -- Sum of fractions not 1?? SumEfrac = " << echeck << " LC Etot " << aLC.energy() << std::endl;
  
  //for (unsigned int j = 0; j < hf.size(); j++) {
  //std::cout << aLC.printHitAndFraction(j) << std::endl;
  //}
  
  aLCvec.push_back(layercluster());
  auto& layercluster_ = aLCvec.back();
  layercluster_.energy_ = aLC.energy();///(double)trkster.vertex_multiplicity(ilc);
  layercluster_.eta_ = aLC.eta();
  layercluster_.phi_ = aLC.phi();
  layercluster_.type_ = aType;
  layercluster_.eminRH_ = EminRH;
  layercluster_.areaminRH_ = areaMinRH;
  layercluster_.nrechitsHD_ = nRHHD;
  layercluster_.efracHD_ = eHD/layercluster_.energy_;

  //get seed information
  bool isSi = recHitTools->isSilicon(aLC.seed());
  auto const itcheck = hitMap.find(aLC.seed());
  if (itcheck != hitMap.end()) {
    const HGCRecHit* hit = itcheck->second;
    //recHitTools->getLayerWithOffset(itcheck->first) < recHitTools->firstLayerBH();

    layercluster_.seedEnergy_ = hit->energy();///(double)trkster.vertex_multiplicity(ilc);
    layercluster_.seedArea_ = recHitTools->getCellArea(itcheck->first);
    layercluster_.seedx_ = recHitTools->getPosition(itcheck->first).x();
    layercluster_.seedy_ = recHitTools->getPosition(itcheck->first).y();
    layercluster_.seedu_ = isSi ? recHitTools->getCell(itcheck->first).first : -1;
    layercluster_.seedv_ = isSi ? recHitTools->getCell(itcheck->first).second : -1;
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

  //if (iterName == "Trk" && layercluster_.nrechits_<3 && abs(layer_) < recHitTools->firstLayerBH()){
  //std::cout << " -- Arf, " << iterName << " nRH = " << layercluster_.nrechits_ << std::endl;
  //}

  layercluster_.isSi_ = isSi;
  if (layercluster_.isSi_) layercluster_.isHD_ = HGCSiliconDetId(aLC.seed()).highDensity();
  else layercluster_.isHD_ = -1;
  layercluster_.layer_ = abs(layer_);
  layercluster_.idxTracksterLC_ = aIdx;
}

void TiCLTreeProducer::initialiseLCTreeVariables(){

  nAllLC = 0;
  all_lc_type.clear();
  all_lc_energy.clear();
  all_lc_eta.clear();
  all_lc_phi.clear();
  all_lc_seedEnergy.clear();
  all_lc_seedArea.clear();
  all_lc_seedx.clear();
  all_lc_seedy.clear();
  all_lc_seedu.clear();
  all_lc_seedv.clear();
  all_lc_eminRH.clear();
  all_lc_areaminRH.clear();
  all_lc_seedEta.clear();
  all_lc_seedPhi.clear();
  all_lc_layer.clear();
  all_lc_isSi.clear();
  all_lc_isHD.clear();
  all_lc_nrechits.clear();
  all_lc_nrechitsHD.clear();
  all_lc_efracHD.clear();
  all_lc_mult.clear();
  all_lc_x.clear();
  all_lc_y.clear();
  all_lc_z.clear();

};

void TiCLTreeProducer::initialiseRHTreeVariables(const unsigned nL){
  for (unsigned i(0); i<3;++i){
    rh_n[i].clear();
    rh_nHD[i].clear();
    rh_nScint[i].clear();
    rh_E[i].clear();
    rh_EHD[i].clear();
    rh_EScint[i].clear();
    rh_n[i].resize(nL,0);
    rh_nHD[i].resize(nL,0);
    rh_nScint[i].resize(nL,0);
    rh_E[i].resize(nL,0);
    rh_EHD[i].resize(nL,0);
    rh_EScint[i].resize(nL,0);
  }
};

void TiCLTreeProducer::fillDoubletsInfo(const ticl::Trackster & thisTrackster,
			   const reco::CaloClusterCollection &  layerClusters,
			   LightTree & myTree
			   ){

  std::vector< std::vector<Triplet> > tripletVec;
  std::vector<Triplet> dummy;
  dummy.reserve(10);
  tripletVec.resize(nLEE,dummy);
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

  const CaloGeometry* geom = &iSetup.getData(tok_geom_);
  recHitTools->setGeometry(*geom);

  const unsigned nLayers = recHitTools->lastLayerBH();

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

  edm::Handle<std::vector<reco::GenParticle>> genParticles;
  iEvent.getByToken(genParticlesToken_, genParticles);

  edm::Handle<std::vector<SimTrack>> simTracks;
  iEvent.getByToken(simTracksToken_, simTracks);

  edm::Handle<std::vector<SimVertex>> simVertices;
  iEvent.getByToken(simVerticesToken_, simVertices);

  const auto& simClustersToRecoColl = iEvent.get(associatorMapSimToRecoToken_);

  //edm::Handle<hgcal::SimToRecoCollectionWithSimClusters> simToRecoCollHandle;
  //iEvent.getByToken(associatorMapSimToRecoToken_, simToRecoCollHandle);
  //const hgcal::SimToRecoCollectionWithSimClusters& simToRecoColl = *simToRecoCollHandle;
  //auto simRecColl = *simToRecoCollHandle;

  edm::Handle<hgcal::RecoToSimCollectionWithSimClusters> recoToSimCollHandle;
  iEvent.getByToken(associatorMapRecoToSimToken_, recoToSimCollHandle);
  //const hgcal::RecoToSimCollectionWithSimClusters& recoToSimColl = *recoToSimCollHandle;
  auto recSimColl = *recoToSimCollHandle; 

  edm::Handle<reco::CaloClusterCollection> layerClusterOLDHandle;
  iEvent.getByToken(hgcalLayerClustersToken_, layerClusterOLDHandle);
  const reco::CaloClusterCollection& lcsold = *layerClusterOLDHandle;

  edm::Handle<reco::CaloClusterCollection> layerClusterNEWHandle;
  iEvent.getByToken(hgcalLayerClustersNEWToken_, layerClusterNEWHandle);
  const reco::CaloClusterCollection& lcsnew = *layerClusterNEWHandle;

  std::map<DetId, const HGCRecHit*> hitMap;
  fillHitMap(hitMap, *recHitHandleEE, *recHitHandleFH, *recHitHandleBH);

  //edm::Handle<edm::ValueMap<std::pair<float, float>>> lcTimeHandle;
  //iEvent.getByToken(hgcalLayerClusterTimeToken_, lcTimeHandle);
  //const auto& lcTime = *lcTimeHandle;

  // init vars
  //recHitTools->getEventSetup(iSetup);

  if (debug) std::cout << " - Processing event " << iEvent.id().event() << std::endl;

  initialiseLCTreeVariables();
  initialiseRHTreeVariables(nLayers);

  rh_evt = iEvent.id().event();
  rh_layer.resize(nLayers,0);
  for (unsigned iL(0);iL<nLayers;++iL){
    rh_layer[iL] = iL;
  }

  //@@  std::vector<std::vector<double> > newSimLCmult;

  // get photon truth info for conversions
  auto mcTruthPhotons = photonMCTruthFinder_.find(*simTracks, *simVertices);
  PhotonMCTruth mcTruthPho;

  irun = iEvent.id().run();
  ievent = iEvent.id().event();
  ilumiblock = iEvent.id().luminosityBlock();

  // Get x,y,z positions of SimVertices (indices 0 and 1 used to determine distance at CaloFace)
  std::vector<float> sv_posX = {};
  std::vector<float> sv_posY = {};
  std::vector<float> sv_posZ = {};
  for ( const auto& sv : *simVertices ) {
    sv_posX.push_back( sv.position().x() );
    sv_posY.push_back( sv.position().y() );
    sv_posZ.push_back( sv.position().z() );
  }
  
  /*const simVertexRVec & v = *simVertices;

  auto distanceV = v[0].position()-v[1].position();
  float distance = ceilf(distanceV.P() * 100)/100;
  float sv12_posX = (sv_posX)[0] - (sv_posX)[1];
  float sv12_posY = (sv_posY)[0] - (sv_posY)[1];
  float sv12_posZ = (sv_posZ)[0] - (sv_posZ)[1];
  float dist_2d = std::sqrt( sv12_posX*sv12_posX + sv12_posY*sv12_posY );
  float dist_3d = std::sqrt( dist_2d* dist_2d + sv12_posZ*sv12_posZ );

  if (fabs(distance-dist_3d) > 0.01) std::cout << " Marco " << distance 
					       << " Rob " << dist_2d 
					       << " " << dist_3d 
					       << std::endl;
  */

  // get CaloParticles
  std::vector<caloparticle> caloparticles;
  std::vector<float> convAbsDzs;
  std::vector<float> convVtxX;
  std::vector<float> convVtxY;
  std::vector<float> convVtxZ;
  std::vector<simcluster> simclust;
//@@  newSimLCmult.clear();

  if (debug) std::cout << " n Caloparticles" << cps.size() << std::endl;
  
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

    bool mcTruthExists = false;
    float minDR = 10.;
    for( auto mcpho : mcTruthPhotons) {
      if( abs((mcpho.fourMomentum().et() - cp.pt()) / cp.pt()) < 0.5 && deltaR(mcpho.fourMomentum().eta(), mcpho.fourMomentum().phi(), cp.eta(), cp.phi()) < minDR ) {
        minDR = deltaR(mcpho.fourMomentum().eta(), mcpho.fourMomentum().phi(), cp.eta(), cp.phi());
        mcTruthPho = mcpho;
        mcTruthExists = true;
      }
    }

    caloparticle tmpcp_;
    tmpcp_.idx_ = idx;
    tmpcp_.pdgid_ = cp.pdgId();
    tmpcp_.energy_ = cp.energy();
    tmpcp_.pt_ = cp.pt();
    tmpcp_.eta_ = cp.eta();
    tmpcp_.phi_ = cp.phi();

    // get the simclusters
    const SimClusterRefVector& simclusters = cp.simClusters();
    tmpcp_.nSC_ = simclusters.size();

    if (debug) std::cout << "CP particle: " << tmpcp_.print() << std::endl;

    for (const auto& it_simc : simclusters) {
      const SimCluster& simc = (*(it_simc));
      const auto& sc_haf = simc.hits_and_fractions();

      simcluster tmpsc_;
      tmpsc_.idx_ = idx;
      tmpsc_.pdgid_ = simc.pdgId();

      auto st = simc.g4Tracks().at(0);
      auto mom = st.getMomentumAtBoundary();
      auto pos = st.getPositionAtBoundary();

      tmpsc_.energy_ = simc.energy();
      tmpsc_.pt_ = simc.pt();
      tmpsc_.eta_ = simc.eta();
      tmpsc_.phi_ = simc.phi();
      tmpsc_.energyAtB_ = mom.E();
      tmpsc_.ptAtB_ = mom.pt();
      tmpsc_.etaAtB_ = mom.eta();
      tmpsc_.phiAtB_ = mom.phi();
      tmpsc_.xAtB_ = pos.x();
      tmpsc_.yAtB_ = pos.y();
      tmpsc_.zAtB_ = pos.z();

      if (debug) std::cout << tmpsc_.print() << std::endl;

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
    convAbsDzs.push_back( mcTruthExists ? mcTruthPho.vertex().vect().z() : -1. );
    convVtxX.push_back( mcTruthExists ? mcTruthPho.vertex().vect().x() : -1. );
    convVtxY.push_back( mcTruthExists ? mcTruthPho.vertex().vect().y() : -1. );
    convVtxZ.push_back( mcTruthExists ? mcTruthPho.vertex().vect().z() : -1. );
    
//@@    //CAMM
//@@    //@FIXME Probably to be removed when moving away from CMSSW_11_3_0_pre3!
//@@    //hack the multiplicity to have proper truth E fraction (not uint8...) -- fixed now in CMSSW master... 
//@@    
//@@    for (unsigned iT(0); iT<nIters; ++iT){
//@@      
//@@      if (iterTypeVec[iT]!="Sim") continue;
//@@
//@@      auto const& tracksters = *(trkstersVec[iT]);
//@@      //std::map<unsigned,double>checkMult;
//@@      //std::pair<std::map<unsigned,double>::iterator,bool> checkMultInsert;
//@@
//@@      for (unsigned its = 0; its < tracksters.size(); ++its) {
//@@	const auto& scRef = simclusters[tracksters[its].seedIndex()];
//@@	const auto& it = simClustersToRecoColl.find(scRef);
//@@
//@@	if (it == simClustersToRecoColl.end()){
//@@	  std::cout << "SimClustToRecoRef Not found ! " << std::endl;
//@@	  continue;
//@@	}
//@@	const auto& lcVec = it->val;
//@@	if (lcVec.empty()){
//@@	  std::cout << " LCVec empty ! " << std::endl;
//@@	  continue;
//@@	}
//@@	std::vector<double> newLCmult;
//@@	std::vector<layercluster> lcsFromTrkster;
//@@	getLCsFromTrkster(iterTypeVec[iT],tracksters[its], lcs, lcsFromTrkster,hitMap);
//@@	std::sort(lcsFromTrkster.begin(), lcsFromTrkster.end(), sortLCsByEnergyAndLayer);
//@@	for (unsigned ilc(0); ilc<lcsFromTrkster.size(); ++ilc){
//@@	  //std::cout << "--- ilc " << ilc << " original idx " << tracksters[its].vertices()[ilc] << std::endl;
//@@	  double fraction = 1.;
//@@	  for (auto const& [lc, energyScorePair] : lcVec) {
//@@	    //auto const & [lc, energyScorePair] = lcVec[tracksters[its].vertices()[ilc]];
//@@	    if (lc.index() != static_cast<unsigned>(lcsFromTrkster[ilc].idxTracksterLC_)) continue;
//@@	    fraction = energyScorePair.first / lc->energy();
//@@	    //if (fraction <1) std::cout << " ---- lcVec idx " << lc.index() << " original idx " << tracksters[its].vertices()[ilc] << " E " << lc->energy() << " frac " << fraction << std::endl;
//@@	    if (fraction > 0) newLCmult.push_back(1. / fraction);
//@@	    //checkMultInsert = checkMult.insert(std::pair<unsigned,double>(lc.index(),fraction));
//@@	    //if (!checkMultInsert.second) checkMultInsert.first->second += fraction;
//@@	  }
//@@	  //newLCmult.push_back(1. / fraction);
//@@	}
//@@	if (newLCmult.size() != lcsFromTrkster.size()) std::cout << " --- size orig " << lcsFromTrkster.size() << " after : " << newLCmult.size() << std::endl;
//@@	assert(newLCmult.size() == lcsFromTrkster.size());
//@@	newSimLCmult.push_back(newLCmult);
//@@      }
//@@      if (newSimLCmult.size() != tracksters.size()) std::cout << " --- size orig " << tracksters.size() << " after : " << newSimLCmult.size() << std::endl;
//@@
//@@      //std::cout << " Check map when filling vector " << std::endl;
//@@      //std::map<unsigned,double>::iterator checkMultIter = checkMult.begin();
//@@      //for (;checkMultIter!=checkMult.end();++checkMultIter){
//@@      //if ((checkMultIter->second) > 1.000001) std::cout << " event " << iEvent.id().event() << " lc " << checkMultIter->first << " sumFrac " << checkMultIter->second << std::endl;
//@@      //}
//@@    }//loop over iters
    
    
  }  // end of looping over the calo particles


  if (debug) std::cout << "number of selected caloparticles: " <<     caloparticles.size() << std::endl;
  if (caloparticles.size()==0) return;


  auto itcheck = hitMap.begin();
  for (;itcheck != hitMap.end();++itcheck) {
    const HGCRecHit* hit = itcheck->second;
    const DetId & id = itcheck->first;
    const unsigned layer = recHitTools->getLayerWithOffset(id);
    bool isSi = recHitTools->isSilicon(id);
    bool isHD = isSi?HGCSiliconDetId(id).highDensity():false;
    if (!isSi) h_rh_SoN_Scint->Fill(layer,hit->signalOverSigmaNoise());
    else if (isHD) h_rh_SoN_HD->Fill(layer,hit->signalOverSigmaNoise());
    else h_rh_SoN_LD->Fill(layer,hit->signalOverSigmaNoise());
    //reject noise hits
    for (unsigned iSN(0); iSN<3;++iSN){
      if (hit->signalOverSigmaNoise()<=SNthresh[iSN]) continue;
      rh_n[iSN][layer-1] += 1;
      rh_E[iSN][layer-1] += hit->energy();
      if (isHD){
	rh_nHD[iSN][layer-1] += 1;
	rh_EHD[iSN][layer-1] += hit->energy();	
      }
      if (!isSi){
	rh_nScint[iSN][layer-1] += 1;
	rh_EScint[iSN][layer-1] += hit->energy();	
      }
    }
  }

  treeAllRH->Fill();


  
  // get the relevant trackster collection
  // LG: For now always use the MergedTrackster
  //   int cp_pdgid_ = 0; if (caloparticles.size()>0) { cp_pdgid_ = caloparticles.at(0).pdgid_; }
  //   std::vector<ticl::Trackster> tracksters = getTracksterCollection(cp_pdgid_,emMCs, mipMCs, hadMCs, mergedMCs);

  // loop over the caloparticles and then find the closest trackster to it

  // keep tracksters [and the corresponding LC] that
  // have at least some ammount of energy from the CP
  //   std::vector<trackster> trksterCollection; trksterCollection.clear();

  //double trackstersEnFromCP = 0.9;
  

  

  all_lc_type.clear();
  all_lc_energy.clear();
  all_lc_eta.clear();
  all_lc_phi.clear();
  all_lc_seedEnergy.clear();
  all_lc_seedArea.clear();
  all_lc_seedx.clear();
  all_lc_seedy.clear();
  all_lc_seedu.clear();
  all_lc_seedv.clear();
  all_lc_eminRH.clear();
  all_lc_areaminRH.clear();
  all_lc_seedEta.clear();
  all_lc_seedPhi.clear();
  all_lc_x.clear();
  all_lc_y.clear();
  all_lc_z.clear();
  all_lc_isSi.clear();
  all_lc_isHD.clear();
  all_lc_layer.clear();
  all_lc_nrechits.clear();
  all_lc_nrechitsHD.clear();
  all_lc_efracHD.clear();
  all_lc_mult.clear();

  std::vector<layercluster> alllcs;
  //RECO collection
  getLCs(lcsold,alllcs,hitMap,0);
  //TICL collection after modification of LC algo adding density.
  getLCs(lcsnew,alllcs,hitMap,1);
  std::sort(alllcs.begin(), alllcs.end(), sortLCsByLayerAndEnergy);

  //one map for each type of LC
  std::map<int,int>lMapLC[2];
  nAllLC = alllcs.size();
  all_lc_type.reserve(nAllLC);
  all_lc_energy.reserve(nAllLC);
  all_lc_eta.reserve(nAllLC);
  all_lc_phi.reserve(nAllLC);
  all_lc_seedEnergy.reserve(nAllLC);
  all_lc_seedArea.reserve(nAllLC);
  all_lc_seedx.reserve(nAllLC);
  all_lc_seedy.reserve(nAllLC);
  all_lc_seedu.reserve(nAllLC);
  all_lc_seedv.reserve(nAllLC);
  all_lc_eminRH.reserve(nAllLC);
  all_lc_areaminRH.reserve(nAllLC);
  all_lc_seedEta.reserve(nAllLC);
  all_lc_seedPhi.reserve(nAllLC);
  all_lc_x.reserve(nAllLC);
  all_lc_y.reserve(nAllLC);
  all_lc_z.reserve(nAllLC);
  all_lc_isSi.reserve(nAllLC);
  all_lc_isHD.reserve(nAllLC);
  all_lc_layer.reserve(nAllLC);
  all_lc_nrechits.reserve(nAllLC);
  all_lc_nrechitsHD.reserve(nAllLC);
  all_lc_efracHD.reserve(nAllLC);
  all_lc_mult.reserve(nLayers);
  
  for (auto const& lc : alllcs) {
    all_lc_type.push_back(lc.type_);
    all_lc_energy.push_back(lc.energy_);
    all_lc_eta.push_back(lc.eta_);
    all_lc_phi.push_back(lc.phi_);
    all_lc_seedEnergy.push_back(lc.seedEnergy_);
    all_lc_seedArea.push_back(lc.seedArea_);
    all_lc_seedx.push_back(lc.seedx_);
    all_lc_seedy.push_back(lc.seedy_);
    all_lc_seedu.push_back(lc.seedu_);
    all_lc_seedv.push_back(lc.seedv_);
    all_lc_eminRH.push_back(lc.eminRH_);
    all_lc_areaminRH.push_back(lc.areaminRH_);
    all_lc_seedEta.push_back(lc.seedEta_);
    all_lc_seedPhi.push_back(lc.seedPhi_);
    all_lc_x.push_back(lc.x_);
    all_lc_y.push_back(lc.y_);
    all_lc_z.push_back(lc.z_);
    all_lc_isSi.push_back(lc.isSi_);
    all_lc_isHD.push_back(lc.isHD_);
    all_lc_layer.push_back(lc.layer_);
    std::pair<std::map<int,int>::iterator,bool> isInserted =  lMapLC[lc.type_].insert(std::pair<int,int>(lc.layer_,1));
    if (!isInserted.second) isInserted.first->second += 1;
    
    all_lc_nrechits.push_back(lc.nrechits_);
    all_lc_nrechitsHD.push_back(lc.nrechitsHD_);
    all_lc_efracHD.push_back(lc.efracHD_);
  }
  
  for (unsigned iTy(0); iTy<2;++iTy){
    for (unsigned iL(0); iL<nLayers;++iL){
      std::map<int,int>::iterator lEle = lMapLC[iTy].find(iL+1);
      if (lEle !=lMapLC[iTy].end()) all_lc_mult.push_back(lEle->second);
      else all_lc_mult.push_back(0);
    }
  }
  
  treeAllLC->Fill();

  if (debug) std::cout << " - nAllLC = " << nAllLC << std::endl;

  for (unsigned iT(0); iT<nIters; ++iT){
    if (debug) std::cout << " -- Processing iter " << iterTypeVec[iT] << std::endl;

    const reco::CaloClusterCollection& lcs = (iterTypeVec[iT].find("NEWLC")==iterTypeVec[iT].npos) ? lcsold : lcsnew;
    edm::Handle<reco::CaloClusterCollection> layerClusterHandle = (iterTypeVec[iT].find("NEWLC")==iterTypeVec[iT].npos) ? layerClusterOLDHandle : layerClusterNEWHandle;

    //fill this here - duplicated across tracksters trees...
    lighttreeVec[iT].initialiseTreeVariables((size_t)irun,(int)ievent,(size_t)ilumiblock); 
    lighttreeVec[iT].fillSVinfo(sv_posX,sv_posY,sv_posZ);
    lighttreeVec[iT].fillCPinfo(caloparticles, convAbsDzs, convVtxX, convVtxY, convVtxZ);
    lighttreeVec[iT].fillSCinfo(simclust);
    
    auto const& tracksters = *(trkstersVec[iT]);

    //std::cout << " Entry " << ievent 
    //	      << " iteration " << iterTypeVec[iT] 
    //	      << " SimLCmult size " << newSimLCmult.size() 
    //	      << " nTS " << tracksters.size()
    //	      << std::endl;
    
    if (debug) std::cout << " --- Number of tracksters: " << tracksters.size() << std::endl;

    // find the tracksters within some DR from the CP

    float maxDrTracksterCP = 0.4;
    int itrksterMin_ = -1;
    for (unsigned int icp = 0; icp < caloparticles.size(); ++icp) {
      
      //AM@TODO
      //Add also CPidx and SCidx of closest CP/SC to trackster...

      std::vector<unsigned> closestTracksters =
	getClosestTrackstersToCPByIndex(iterTypeVec[iT],caloparticles[icp], lcs, tracksters, maxDrTracksterCP,hitMap);

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
      getLCsFromTrkster(iterTypeVec[iT],tracksters[its], lcs, lcsFromTrkster,hitMap);
      std::sort(lcsFromTrkster.begin(), lcsFromTrkster.end(), sortLCsByEnergyAndLayer);
      lcsFromTrksters.push_back(lcsFromTrkster);
    }
    
//@@    if (iterTypeVec[iT]=="Sim") lighttreeVec[iT].fillTSinfo(ievent,iterTypeVec[iT],tracksters,nLayers,lcsFromTrksters,layerClusterHandle,recSimColl,newSimLCmult);
//@@    else {
    std::vector<std::vector<double> > dummy;

    lighttreeVec[iT].fillTSinfo(ievent,iterTypeVec[iT],tracksters,nLayers,lcsFromTrksters,layerClusterHandle,recSimColl,dummy);
//@@    }
    //if (itrksterMin_>=0 && fillTripletsInfo) fillDoubletsInfo(tracksters[itrksterMin_],lcs,lighttreeVec[iT]);
    
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
