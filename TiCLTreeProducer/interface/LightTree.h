#ifndef __LightTree__hh
#define __LightTree__hh

#include <string>
#include <vector>

#include <TROOT.h>
//#include <TFile.h>
#include "TTree.h"
#include "TBranch.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "SimDataFormats/Associations/interface/LayerClusterToSimClusterAssociator.h"
#include "SimDataFormats/Associations/interface/LayerClusterToSimClusterAssociatorBaseImpl.h"


// Common Data Formats used by the ntuple
#include "ICTICLAnalysis/TiCLTreeProducer/interface/CommonDataFormats.h"

class LightTree {
 public:
  
  LightTree();
  ~LightTree();
  
  void makeTree(edm::Service<TFileService> & aFile,
		const std::string & aIterName,
		const int aFillTripletsInfo,
		const unsigned aDebug);
  
  
  void initialiseTreeVariables(const size_t irun,
			       const int ievent,
			       const size_t ilumiblock
			       );
  
  void fillTriplets(const std::vector< std::vector<Triplet> > & aTripletVec);
  
  
  void fillCPinfo(const std::vector<caloparticle> & caloparticles, const std::vector<float> dzs);

  void fillSCinfo(const std::vector<simcluster> & simclusters);

  void fillTSinfo(const unsigned evtNum,
		  const std::string iterName,
		  const std::vector<ticl::Trackster> & tracksters,
		  const unsigned nLayers,
		  const std::vector<std::vector<layercluster> > & lcsFromTrkster,
		  const edm::Handle<reco::CaloClusterCollection> & layerClusterHandle,
		  const hgcal::RecoToSimCollectionWithSimClusters & recSimColl,
		  const std::vector<std::vector<double> > & newSimLCmult
		  );

  
  inline void fillCPEfraction(const double & trksterCPEnDiff,
			      const unsigned icp){
    if (cp_missingEnergyFraction.size()<icp+1) return;
    cp_missingEnergyFraction[icp] = trksterCPEnDiff;
  };

  inline void fillOutputTree(){
    outputTree->Fill();
  };
  
 private:

  TTree* outputTree;
  unsigned mDebug;

  unsigned nL;
  size_t run, lumi;
  int event;
  double weight;

  int nTS;
  std::vector<int> ts_CPidx;
  std::vector<int> ts_seedIdx;
  std::vector<double> ts_energy;
  std::vector<double> ts_emEnergy;
  std::vector<double> ts_regEnergy;
  std::vector<double> ts_sigma1;
  std::vector<double> ts_sigma2;
  std::vector<double> ts_sigma3;
  std::vector<double> ts_BCx;
  std::vector<double> ts_BCy;
  std::vector<double> ts_BCz;
  std::vector<double> ts_eta_PCA;
  std::vector<double> ts_phi_PCA;
  std::vector<double> ts_eta_fromLC;
  std::vector<double> ts_phi_fromLC;
  std::vector<double> ts_photon_proba;
  std::vector<double> ts_ele_proba;
  std::vector<double> ts_mu_proba;
  std::vector<double> ts_pi0_proba;
  std::vector<double> ts_chHad_proba;
  std::vector<double> ts_neHad_proba;
  std::vector<double> ts_ambg_proba;
  std::vector<double> ts_unkwn_proba;
  std::vector<int> ts_nLC;
  std::vector<int> ts_firstLayer;
  std::vector<int> ts_lastLayer;
  std::vector<int> ts_outInHopsPerformed;

  int nCP;
  std::vector<double> cp_missingEnergyFraction;
  std::vector<double> cp_energy;
  std::vector<int> cp_pdgid;
  std::vector<int> cp_nSC;
  std::vector<double> cp_pt;
  std::vector<double> cp_eta;
  std::vector<double> cp_phi;
  std::vector<double> cp_convAbsDz;

  int nSC;
  std::vector<int> sc_CPidx;
  std::vector<double> sc_energy;
  std::vector<double> sc_energyAtB;
  std::vector<int> sc_pdgid;
  std::vector<double> sc_pt;
  std::vector<double> sc_eta;
  std::vector<double> sc_phi;
  std::vector<double> sc_ptAtB;
  std::vector<double> sc_etaAtB;
  std::vector<double> sc_phiAtB;
  std::vector<double> sc_xAtB;
  std::vector<double> sc_yAtB;
  std::vector<double> sc_zAtB;

  int nLC;
  std::vector<int> lc_idx;
  std::vector<int> lc_TSidx;
  std::vector<double> lc_energy;
  std::vector<double> lc_eta;
  std::vector<double> lc_phi;
  std::vector<double> lc_seedEnergy;
  std::vector<double> lc_seedEta;
  std::vector<double> lc_seedPhi;
  std::vector<double> lc_x;
  std::vector<double> lc_y;
  std::vector<double> lc_z;
  std::vector<int> lc_algo;
  std::vector<int> lc_isSi;
  std::vector<int> lc_layer;
  std::vector<int> lc_nrechits;
  std::vector<double> lc_tsMult;
  std::vector<int> lc_mult;
  std::vector<int> lc_nSC;
  std::vector<int> lc_pdgid;
  std::vector<int> lc_SCidx[5];
  std::vector<double> lc_SCefrac[5];

  int nTriplets[28];
  std::vector<int> triplets_layerA[28];
  std::vector<int> triplets_layerC[28];
  std::vector<double> triplets_energyA[28];
  std::vector<double> triplets_energyB[28];
  std::vector<double> triplets_energyC[28];
  std::vector<double> triplets_etaB[28];
  std::vector<double> triplets_cosBeta[28];
  std::vector<double> triplets_cosAlphaInner[28];
  std::vector<double> triplets_cosAlphaOuter[28];
  std::vector<int> triplets_inner_in_links[28];
  std::vector<int> triplets_inner_out_links[28];
  std::vector<int> triplets_outer_in_links[28];
  std::vector<int> triplets_outer_out_links[28];

};//class

#endif
