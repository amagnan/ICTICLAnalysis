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

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/LorentzVector.h"

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

#include "TSystem.h"
#include "TVector3.h"
#include "TLorentzVector.h"


class CellAreaChecker : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit CellAreaChecker(const edm::ParameterSet&);
  ~CellAreaChecker() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  virtual void fillHitMap(std::map<DetId, const HGCRecHit*>& hitMap,
                          const HGCRecHitCollection& rechitsEE,
                          const HGCRecHitCollection& rechitsFH,
                          const HGCRecHitCollection& rechitsBH) const;


  void initialiseTreeVariables();

  std::shared_ptr<hgcal::RecHitTools> recHitTools;
  
  // ----------member data ---------------------------
  const edm::ESGetToken<CaloGeometry, CaloGeometryRecord> tok_geom_;

  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsEEToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsFHToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsBHToken_;

  unsigned debug;


  TTree *treeHits;
  size_t run, lumi;
  int event;
  int nHits;
  std::vector<double> hit_energy;
  std::vector<double> hit_eta;
  std::vector<double> hit_phi;
  std::vector<double> hit_area;
  std::vector<double> hit_x;
  std::vector<double> hit_y;
  std::vector<double> hit_u;
  std::vector<double> hit_v;
  std::vector<double> hit_z;
  std::vector<int> hit_isSi;
  std::vector<int> hit_isHD;
  std::vector<int> hit_layer;

};


CellAreaChecker::CellAreaChecker(const edm::ParameterSet& iConfig):
  tok_geom_(esConsumes<CaloGeometry, CaloGeometryRecord>()),
  hgcalRecHitsEEToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsEE"))),
  hgcalRecHitsFHToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsFH"))),
  hgcalRecHitsBHToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsBH")))
{
  
  recHitTools.reset(new hgcal::RecHitTools());
  //now do what ever initialization is needed
  usesResource("TFileService");
  edm::Service<TFileService> file;
  
  debug = iConfig.getParameter<int>("Debug");
  
  treeHits = file->make<TTree>("treeHits", "treeHits");
  treeHits->Branch("run", &run, "run/I");
  treeHits->Branch("event", &event, "event/I");
  treeHits->Branch("lumi", &lumi, "lumi/I");
  treeHits->Branch("nHits", &nHits, "nHits/I");
  treeHits->Branch("hit_energy", &hit_energy);
  treeHits->Branch("hit_eta", &hit_eta);
  treeHits->Branch("hit_phi", &hit_phi);
  treeHits->Branch("hit_area", &hit_area);
  treeHits->Branch("hit_isSi", &hit_isSi);
  treeHits->Branch("hit_isHD", &hit_isHD);
  treeHits->Branch("hit_layer", &hit_layer);
  treeHits->Branch("hit_x", &hit_x);
  treeHits->Branch("hit_y", &hit_y);
  treeHits->Branch("hit_u", &hit_u);
  treeHits->Branch("hit_v", &hit_v);
  treeHits->Branch("hit_z", &hit_z);
  
}

CellAreaChecker::~CellAreaChecker() {}

void CellAreaChecker::fillHitMap(std::map<DetId, const HGCRecHit*>& hitMap,
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
}  // end of CellAreaChecker::fillHitMap


void CellAreaChecker::initialiseTreeVariables(){

  nHits = 0;
  hit_energy.clear();
  hit_eta.clear();
  hit_phi.clear();
  hit_area.clear();
  hit_layer.clear();
  hit_isSi.clear();
  hit_isHD.clear();
  hit_x.clear();
  hit_y.clear();
  hit_u.clear();
  hit_v.clear();
  hit_z.clear();

};

void CellAreaChecker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  const CaloGeometry* geom = &iSetup.getData(tok_geom_);
  recHitTools->setGeometry(*geom);

  //unsigned nLayers = recHitTools->lastLayerBH();

  edm::Handle<HGCRecHitCollection> recHitHandleEE;
  iEvent.getByToken(hgcalRecHitsEEToken_, recHitHandleEE);

  edm::Handle<HGCRecHitCollection> recHitHandleFH;
  iEvent.getByToken(hgcalRecHitsFHToken_, recHitHandleFH);

  edm::Handle<HGCRecHitCollection> recHitHandleBH;
  iEvent.getByToken(hgcalRecHitsBHToken_, recHitHandleBH);


  std::map<DetId, const HGCRecHit*> hitMap;
  fillHitMap(hitMap, *recHitHandleEE, *recHitHandleFH, *recHitHandleBH);


  if (debug) std::cout << " - Processing event " << iEvent.id().event() << std::endl;
  
  initialiseTreeVariables();

  run = iEvent.id().run();
  event = iEvent.id().event();
  lumi = iEvent.id().luminosityBlock();
  
  nHits = hitMap.size();
  hit_energy.reserve(nHits);
  hit_eta.reserve(nHits);
  hit_phi.reserve(nHits);
  hit_area.reserve(nHits);
  hit_x.reserve(nHits);
  hit_y.reserve(nHits);
  hit_u.reserve(nHits);
  hit_v.reserve(nHits);
  hit_z.reserve(nHits);
  hit_isSi.reserve(nHits);
  hit_isHD.reserve(nHits);
  hit_layer.reserve(nHits);
  
  //get seed information
  auto itcheck = hitMap.begin();
  for (;itcheck != hitMap.end();++itcheck) {
    const HGCRecHit* hit = itcheck->second;
    const DetId & id = itcheck->first;

    bool isSi = recHitTools->getLayerWithOffset(id) < recHitTools->firstLayerBH();

    hit_energy.push_back(hit->energy());
    hit_eta.push_back(recHitTools->getEta(id));
    hit_phi.push_back(recHitTools->getPhi(id));
    hit_area.push_back(recHitTools->getCellArea(id));
    hit_x.push_back(recHitTools->getPosition(id).x());
    hit_y.push_back(recHitTools->getPosition(id).y());
    hit_u.push_back(isSi?recHitTools->getCell(id).first:-1);
    hit_v.push_back(isSi?recHitTools->getCell(id).second:-1);
    hit_z.push_back(recHitTools->getPosition(id).z());
    hit_isSi.push_back(isSi);
    hit_isHD.push_back(isSi?HGCSiliconDetId(id).highDensity():-1);
    hit_layer.push_back(recHitTools->getLayerWithOffset(id));
  }
  
  treeHits->Fill();

  if (debug) std::cout << " - nHits = " << nHits << std::endl;


}//analyze

// ------------ method called once each job just before starting event loop  ------------
void CellAreaChecker::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void CellAreaChecker::endJob() {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void CellAreaChecker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CellAreaChecker);
