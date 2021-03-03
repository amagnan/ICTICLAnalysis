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


// Common Data Formats used by the ntuple
#include "ICTICLAnalysis/TiCLTreeProducer/interface/CommonDataFormats.h"

class LightTree {
 public:
  
  LightTree();
  ~LightTree();
  
  void makeTree(edm::Service<TFileService> & aFile,
		const std::string & aIterName,
		const bool aFillTripletsInfo);
  
  
  void initialiseTreeVariables(const size_t irun,
			       const int ievent,
			       const size_t ilumiblock
			       );
  
  void fillTriplets(const std::vector< std::vector<Triplet> > & aTripletVec);
  
  
  void fillCPinfo(const std::vector<caloparticle> & caloparticles,
		  const int icp);
  
  void fillTSinfo(const std::vector<ticl::Trackster> & tracksters,
		  const int itrksterMin,
		  const std::vector<layercluster> & lcsFromClosestTrksterToCP);
  
  inline void fillCPEfraction(const double & trksterCPEnDiff){
    cp_missingEnergyFraction = trksterCPEnDiff;
  };
  
  inline void fillOutputTree(){
    outputTree->Fill();
  };
  
 private:

  TTree* outputTree;

  unsigned nL;
  size_t run, lumi;
  int event;
  double weight;

  int nTS;
  int trackster;
  double ts_energy;
  double ts_emEnergy;
  double ts_regEnergy;
  double ts_sigma1;
  double ts_sigma2;
  double ts_sigma3;
  double ts_BCx;
  double ts_BCy;
  double ts_BCz;
  double ts_eta_PCA;
  double ts_phi_PCA;
  double ts_eta_fromLC;
  double ts_phi_fromLC;
  double ts_photon_proba;
  double ts_ele_proba;
  double ts_mu_proba;
  double ts_pi0_proba;
  double ts_chHad_proba;
  double ts_neHad_proba;
  double ts_ambg_proba;
  double ts_unkwn_proba;
  int ts_firstLayer;
  int ts_lastLayer;
  int ts_outInHopsPerformed;

  int nCP;
  double cp_missingEnergyFraction;
  double cp_energy;
  int cp_pdgid;
  double cp_pt;
  double cp_eta;
  double cp_phi;

  
  int nLC;
  std::vector<double> lc_energy;
  std::vector<double> lc_eta;
  std::vector<double> lc_phi;
  std::vector<double> lc_x;
  std::vector<double> lc_y;
  std::vector<double> lc_z;
  std::vector<int> lc_algo;
  std::vector<int> lc_layer;
  std::vector<int> lc_nrechits;
  std::vector<int> lc_tsMult;
  std::vector<int> lc_mult;

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
