#include "ICTICLAnalysis/TiCLTreeProducer/interface/LightTree.h"
#include "SimDataFormats/CaloAnalysis/interface/SimClusterFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"

LightTree::LightTree(){
  //limit triplets to EE only
  nL=28;
};
LightTree::~LightTree(){
  lc_idx.clear();
  lc_TSidx.clear();
  lc_energy.clear();
  lc_eta.clear();
  lc_phi.clear();
  lc_seedEnergy.clear();
  lc_seedEta.clear();
  lc_seedPhi.clear();
  lc_x.clear();
  lc_y.clear();
  lc_z.clear();
  lc_algo.clear();
  lc_layer.clear();
  lc_isSi.clear();
  lc_nrechits.clear();
  lc_tsMult.clear();
  lc_mult.clear();
  lc_nSC.clear();
  lc_pdgid.clear();
  for (unsigned isc(0); isc<5; ++isc){//loop on layers   
    lc_SCidx[isc].clear();
    lc_SCefrac[isc].clear();
  }

  for (unsigned iL(0); iL<nL; ++iL){//loop on layers   
    nTriplets[iL] = 0;
    triplets_layerA[iL].clear();
    triplets_layerC[iL].clear();
    triplets_energyA[iL].clear();
    triplets_energyB[iL].clear();
    triplets_energyC[iL].clear();
    triplets_etaB[iL].clear();
    triplets_cosBeta[iL].clear();
    triplets_cosAlphaInner[iL].clear();
    triplets_cosAlphaOuter[iL].clear();
    triplets_inner_in_links[iL].clear();
    triplets_inner_out_links[iL].clear();
    triplets_outer_in_links[iL].clear();
    triplets_outer_out_links[iL].clear();
  }
};

void LightTree::makeTree(edm::Service<TFileService> & aFile,
			 const std::string & aIterName,
			 const int aFillTripletsInfo,
			 const unsigned aDebug){

  mDebug = aDebug;

  outputTree = aFile->make<TTree>(("TSTree_"+aIterName).c_str(), ("tree"+aIterName).c_str());
    
  outputTree->Branch("run", &run, "run/I");
  outputTree->Branch("event", &event, "event/I");
  outputTree->Branch("lumi", &lumi, "lumi/I");
  outputTree->Branch("weight", &weight, "weight/D");
  outputTree->Branch("nTS", &nTS, "nTS/I");
  outputTree->Branch("ts_CPidx", &ts_CPidx);
  outputTree->Branch("ts_seedIdx", &ts_seedIdx);
  outputTree->Branch("ts_energy", &ts_energy);
  outputTree->Branch("ts_emEnergy", &ts_emEnergy);
  outputTree->Branch("ts_regEnergy", &ts_regEnergy);
  outputTree->Branch("ts_sigma1", &ts_sigma1);
  outputTree->Branch("ts_sigma2", &ts_sigma2);
  outputTree->Branch("ts_sigma3", &ts_sigma3);
  outputTree->Branch("ts_BCx", &ts_BCx);
  outputTree->Branch("ts_BCy", &ts_BCy);
  outputTree->Branch("ts_BCz", &ts_BCz);
  outputTree->Branch("ts_eta_PCA", &ts_eta_PCA);
  outputTree->Branch("ts_phi_PCA", &ts_phi_PCA);
  outputTree->Branch("ts_eta_fromLC", &ts_eta_fromLC);
  outputTree->Branch("ts_phi_fromLC", &ts_phi_fromLC);
  outputTree->Branch("ts_photon_proba", &ts_photon_proba);
  outputTree->Branch("ts_ele_proba", &ts_ele_proba);
  outputTree->Branch("ts_mu_proba", &ts_mu_proba);
  outputTree->Branch("ts_pi0_proba", &ts_pi0_proba);
  outputTree->Branch("ts_chHad_proba", &ts_chHad_proba);
  outputTree->Branch("ts_neHad_proba", &ts_neHad_proba);
  outputTree->Branch("ts_ambg_proba", &ts_ambg_proba);
  outputTree->Branch("ts_unkwn_proba", &ts_unkwn_proba);
  outputTree->Branch("ts_nLC", &ts_nLC);
  outputTree->Branch("ts_firstLayer", &ts_firstLayer);
  outputTree->Branch("ts_lastLayer", &ts_lastLayer);
  outputTree->Branch("ts_outInHopsPerformed", &ts_outInHopsPerformed);
  
  outputTree->Branch("nCP", &nCP, "nCP/I");
  outputTree->Branch("cp_nSC", &cp_nSC);
  outputTree->Branch("cp_missingEnergyFraction", &cp_missingEnergyFraction);
  outputTree->Branch("cp_energy", &cp_energy);
  outputTree->Branch("cp_pt", &cp_pt);
  outputTree->Branch("cp_eta", &cp_eta);
  outputTree->Branch("cp_phi", &cp_phi);
  outputTree->Branch("cp_convAbsDz", &cp_convAbsDz);
  outputTree->Branch("cp_vtxX", &cp_vtxX);
  outputTree->Branch("cp_vtxY", &cp_vtxY);
  outputTree->Branch("cp_vtxZ", &cp_vtxZ);
  outputTree->Branch("cp_pdgid", &cp_pdgid);

  outputTree->Branch("nSV", &nSV, "nSV/I");
  outputTree->Branch("sv_posX", &sv_posX);
  outputTree->Branch("sv_posY", &sv_posY);
  outputTree->Branch("sv_posZ", &sv_posZ);

  outputTree->Branch("nSC", &nSC, "nSC/I");
  outputTree->Branch("sc_CPidx", &sc_CPidx);
  outputTree->Branch("sc_energy", &sc_energy);
  outputTree->Branch("sc_pt", &sc_pt);
  outputTree->Branch("sc_eta", &sc_eta);
  outputTree->Branch("sc_phi", &sc_phi);
  outputTree->Branch("sc_energyAtB", &sc_energyAtB);
  outputTree->Branch("sc_ptAtB", &sc_ptAtB);
  outputTree->Branch("sc_etaAtB", &sc_etaAtB);
  outputTree->Branch("sc_phiAtB", &sc_phiAtB);
  outputTree->Branch("sc_xAtB", &sc_xAtB);
  outputTree->Branch("sc_yAtB", &sc_yAtB);
  outputTree->Branch("sc_zAtB", &sc_zAtB);
  outputTree->Branch("sc_pdgid", &sc_pdgid);
  
  outputTree->Branch("nLC", &nLC, "nLC/I");
  outputTree->Branch("lc_idx", &lc_idx);
  outputTree->Branch("lc_TSidx", &lc_TSidx);
  outputTree->Branch("lc_energy", &lc_energy);
  outputTree->Branch("lc_eta", &lc_eta);
  outputTree->Branch("lc_phi", &lc_phi);
  outputTree->Branch("lc_seedEnergy", &lc_seedEnergy);
  outputTree->Branch("lc_seedEta", &lc_seedEta);
  outputTree->Branch("lc_seedPhi", &lc_seedPhi);
  outputTree->Branch("lc_x", &lc_x);
  outputTree->Branch("lc_y", &lc_y);
  outputTree->Branch("lc_z", &lc_z);
  outputTree->Branch("lc_algo", &lc_algo);
  outputTree->Branch("lc_layer", &lc_layer);
  outputTree->Branch("lc_isSi", &lc_isSi);
  outputTree->Branch("lc_nrechits", &lc_nrechits);
  outputTree->Branch("lc_tsMult", &lc_tsMult);
  outputTree->Branch("lc_mult", &lc_mult);
  outputTree->Branch("lc_nSC", &lc_nSC);
  outputTree->Branch("lc_pdgid", &lc_pdgid);
  for (unsigned isc(0); isc<5; ++isc){//loop on layers   
    std::ostringstream lLabel;
    lLabel << "lc_SCidx_" << isc;
    outputTree->Branch(lLabel.str().c_str(), &lc_SCidx[isc]);
    lLabel.str("");
    lLabel << "lc_SCefrac_" << isc;
    outputTree->Branch(lLabel.str().c_str(), &lc_SCefrac[isc]);
  }

  if (aFillTripletsInfo>0){
    for (unsigned iL(0); iL<nL; ++iL){//loop on layers   
      std::ostringstream lName;
      lName << "nTriplets_" << iL+1;
      outputTree->Branch(lName.str().c_str(), &nTriplets[iL]);
      lName.str("");
      lName << "triplets_layerA_" << iL+1;
      outputTree->Branch(lName.str().c_str(), &triplets_layerA[iL]);
      lName.str("");
      lName << "triplets_layerC_" << iL+1;
      outputTree->Branch(lName.str().c_str(), &triplets_layerC[iL]);
      lName.str("");
      lName << "triplets_energyA_" << iL+1;
      outputTree->Branch(lName.str().c_str(), &triplets_energyA[iL]);
      lName.str("");
      lName << "triplets_energyB_" << iL+1;
      outputTree->Branch(lName.str().c_str(), &triplets_energyB[iL]);
      lName.str("");
      lName << "triplets_energyC_" << iL+1;
      outputTree->Branch(lName.str().c_str(), &triplets_energyC[iL]);
      lName.str("");
      lName << "triplets_etaB_" << iL+1;
      outputTree->Branch(lName.str().c_str(), &triplets_etaB[iL]);
      lName.str("");
      lName << "triplets_cosBeta_" << iL+1;
      outputTree->Branch(lName.str().c_str(), &triplets_cosBeta[iL]);
      lName.str("");
      lName << "triplets_cosAlphaInner_" << iL+1;
      outputTree->Branch(lName.str().c_str(), &triplets_cosAlphaInner[iL]);
      lName.str("");
      lName << "triplets_cosAlphaOuter_" << iL+1;
      outputTree->Branch(lName.str().c_str(), &triplets_cosAlphaOuter[iL]);
      lName.str("");
      lName << "triplets_inner_in_links_" << iL+1;
      outputTree->Branch(lName.str().c_str(), &triplets_inner_in_links[iL]);
      lName.str("");
      lName << "triplets_inner_out_links_" << iL+1;
      outputTree->Branch(lName.str().c_str(), &triplets_inner_out_links[iL]);
      lName.str("");
      lName << "triplets_outer_in_links_" << iL+1;
      outputTree->Branch(lName.str().c_str(), &triplets_outer_in_links[iL]);
      lName.str("");
      lName << "triplets_outer_out_links_" << iL+1;
      outputTree->Branch(lName.str().c_str(), &triplets_outer_out_links[iL]);
    }
  } 
}



void LightTree::initialiseTreeVariables(const size_t irun,
					const int ievent,
					const size_t ilumiblock
					){
  run = irun;
  event = ievent;
  lumi = ilumiblock;
  
  weight = 0;
  
  nTS = 0;
  ts_CPidx.clear();
  ts_seedIdx.clear();
  ts_energy.clear();
  ts_emEnergy.clear();
  ts_regEnergy.clear();
  ts_sigma1.clear();
  ts_sigma2.clear();
  ts_sigma3.clear();
  ts_BCx.clear();
  ts_BCy.clear();
  ts_BCz.clear();
  ts_eta_PCA.clear();
  ts_phi_PCA.clear();
  ts_eta_fromLC.clear();
  ts_phi_fromLC.clear();
  ts_photon_proba.clear();
  ts_ele_proba.clear();
  ts_mu_proba.clear();
  ts_pi0_proba.clear();
  ts_chHad_proba.clear();
  ts_neHad_proba.clear();
  ts_ambg_proba.clear();
  ts_unkwn_proba.clear();
  ts_nLC.clear();
  ts_firstLayer.clear();
  ts_lastLayer.clear();
  ts_outInHopsPerformed.clear();
  
  nCP = 0;
  cp_nSC.clear();
  cp_missingEnergyFraction.clear();
  cp_energy.clear();
  cp_pdgid.clear();
  cp_pt.clear();
  cp_eta.clear();
  cp_phi.clear();
  cp_convAbsDz.clear();
  cp_vtxX.clear();
  cp_vtxY.clear();
  cp_vtxZ.clear();

  nSV = 0;
  sv_posX.clear();
  sv_posY.clear();
  sv_posZ.clear();

  nSC = 0;
  sc_CPidx.clear();
  sc_energy.clear();
  sc_energyAtB.clear();
  sc_pdgid.clear();
  sc_pt.clear();
  sc_eta.clear();
  sc_phi.clear();
  sc_ptAtB.clear();
  sc_etaAtB.clear();
  sc_phiAtB.clear();
  sc_xAtB.clear();
  sc_yAtB.clear();
  sc_zAtB.clear();
    
  nLC = 0;
  lc_idx.clear();
  lc_TSidx.clear();
  lc_energy.clear();
  lc_eta.clear();
  lc_phi.clear();
  lc_seedEnergy.clear();
  lc_seedEta.clear();
  lc_seedPhi.clear();
  lc_x.clear();
  lc_y.clear();
  lc_z.clear();
  lc_algo.clear();
  lc_layer.clear();
  lc_isSi.clear();
  lc_nrechits.clear();
  lc_tsMult.clear();
  lc_mult.clear();
  lc_nSC.clear();
  lc_pdgid.clear();
  for (unsigned isc(0); isc<5; ++isc){//loop on layers   
    lc_SCidx[isc].clear();
    lc_SCefrac[isc].clear();
  }

  for (unsigned iL(0); iL<nL; ++iL){//loop on layers   
    nTriplets[iL] = 0;
    triplets_layerA[iL].clear();
    triplets_layerC[iL].clear();
    triplets_energyA[iL].clear();
    triplets_energyB[iL].clear();
    triplets_energyC[iL].clear();
    triplets_etaB[iL].clear();
    triplets_cosBeta[iL].clear();
    triplets_cosAlphaInner[iL].clear();
    triplets_cosAlphaOuter[iL].clear();
    triplets_inner_in_links[iL].clear();
    triplets_inner_out_links[iL].clear();
    triplets_outer_in_links[iL].clear();
    triplets_outer_out_links[iL].clear();
  }
  
}

void LightTree::fillTriplets(const std::vector< std::vector<Triplet> > & aTripletVec){

    for (unsigned iL(0); iL<nL; ++iL){//loop on layers   
      //std::cout << " -- Filling trackster triplets, layer " << iL+1 << " nTriplets=" << aTripletVec[iL].size() << std::endl;
      nTriplets[iL] = aTripletVec[iL].size();

      triplets_layerA[iL].clear();
      triplets_layerC[iL].clear();
      triplets_energyA[iL].clear();
      triplets_energyB[iL].clear();
      triplets_energyC[iL].clear();
      triplets_etaB[iL].clear();
      triplets_cosBeta[iL].clear();
      triplets_cosAlphaInner[iL].clear();
      triplets_cosAlphaOuter[iL].clear();
      triplets_inner_in_links[iL].clear();
      triplets_inner_out_links[iL].clear();
      triplets_outer_in_links[iL].clear();
      triplets_outer_out_links[iL].clear();
      triplets_layerA[iL].reserve(aTripletVec[iL].size());
      triplets_layerC[iL].reserve(aTripletVec[iL].size());
      triplets_energyA[iL].reserve(aTripletVec[iL].size());
      triplets_energyB[iL].reserve(aTripletVec[iL].size());
      triplets_energyC[iL].reserve(aTripletVec[iL].size());
      triplets_etaB[iL].reserve(aTripletVec[iL].size());
      triplets_cosBeta[iL].reserve(aTripletVec[iL].size());
      triplets_cosAlphaInner[iL].reserve(aTripletVec[iL].size());
      triplets_cosAlphaOuter[iL].reserve(aTripletVec[iL].size());
      triplets_inner_in_links[iL].reserve(aTripletVec[iL].size());
      triplets_inner_out_links[iL].reserve(aTripletVec[iL].size());
      triplets_outer_in_links[iL].reserve(aTripletVec[iL].size());
      triplets_outer_out_links[iL].reserve(aTripletVec[iL].size());
      for ( const auto & iTriplet : aTripletVec[iL]){
	triplets_layerA[iL].push_back(iTriplet.layerA_);
	triplets_layerC[iL].push_back(iTriplet.layerC_);
	triplets_energyA[iL].push_back(iTriplet.eA_);
	triplets_energyB[iL].push_back(iTriplet.eB_);
	triplets_energyC[iL].push_back(iTriplet.eC_);
	triplets_etaB[iL].push_back(iTriplet.etaB_);
	triplets_cosBeta[iL].push_back(iTriplet.cosBeta_);
	triplets_cosAlphaInner[iL].push_back(iTriplet.cosAlphaInner_);
	triplets_cosAlphaOuter[iL].push_back(iTriplet.cosAlphaOuter_);
	triplets_inner_in_links[iL].push_back(iTriplet.inner_in_links_);
	triplets_inner_out_links[iL].push_back(iTriplet.inner_out_links_);
	triplets_outer_in_links[iL].push_back(iTriplet.outer_in_links_);
	triplets_outer_out_links[iL].push_back(iTriplet.outer_out_links_);
      }
    }
    
}


void LightTree::fillCPinfo(const std::vector<caloparticle> & caloparticles,
			   const std::vector<float> dzs,
			   const std::vector<float> cpVtxX,
			   const std::vector<float> cpVtxY,
			   const std::vector<float> cpVtxZ){
  nCP = caloparticles.size();
  for (int icp = 0; icp < nCP; ++icp) {
    cp_nSC.push_back(caloparticles[icp].nSC_);
    cp_energy.push_back(caloparticles[icp].energy_);
    cp_pt.push_back(caloparticles[icp].pt_);
    cp_eta.push_back(caloparticles[icp].eta_);
    cp_phi.push_back(caloparticles[icp].phi_);
    cp_convAbsDz.push_back(dzs[icp]);
    cp_vtxX.push_back(cpVtxX[icp]);
    cp_vtxY.push_back(cpVtxY[icp]);
    cp_vtxZ.push_back(cpVtxZ[icp]);
    cp_pdgid.push_back(caloparticles[icp].pdgid_);
    cp_missingEnergyFraction.push_back(0);
  }
}

void LightTree::fillSVinfo(const std::vector<float>& svPosX,
			   const std::vector<float>& svPosY,
			   const std::vector<float>& svPosZ){
  nSV = svPosX.size();
  for (int idx = 0; idx < nSV; ++idx) {
    sv_posX.push_back(svPosX[idx]);
    sv_posY.push_back(svPosY[idx]);
    sv_posZ.push_back(svPosZ[idx]);
  }
}

void LightTree::fillSCinfo(const std::vector<simcluster> & simclusters){
  nSC = simclusters.size();
  for (int isc = 0; isc < nSC; ++isc) {
    sc_CPidx.push_back(simclusters[isc].idx_);
    sc_energy.push_back(simclusters[isc].energy_);
    sc_pt.push_back(simclusters[isc].pt_);
    sc_eta.push_back(simclusters[isc].eta_);
    sc_phi.push_back(simclusters[isc].phi_);
    sc_energyAtB.push_back(simclusters[isc].energyAtB_);
    sc_ptAtB.push_back(simclusters[isc].ptAtB_);
    sc_etaAtB.push_back(simclusters[isc].etaAtB_);
    sc_phiAtB.push_back(simclusters[isc].phiAtB_);
    sc_xAtB.push_back(simclusters[isc].xAtB_);
    sc_yAtB.push_back(simclusters[isc].yAtB_);
    sc_zAtB.push_back(simclusters[isc].zAtB_);
    sc_pdgid.push_back(simclusters[isc].pdgid_);
  }
}

void LightTree::fillTSinfo(const unsigned evtNum,
			   const std::string iterName,
			   const std::vector<ticl::Trackster> & tracksters,
			   const unsigned nLayers,
			   const std::vector<std::vector<layercluster>> & lcsFromTrkster,
			   const edm::Handle<reco::CaloClusterCollection> & layerClusterHandle,
			   const hgcal::RecoToSimCollectionWithSimClusters & recSimColl,
			   const std::vector<std::vector<double> > & newSimLCmult)
{
  nTS = tracksters.size();
  nLC = 0;
  if(lcsFromTrkster.size()!=static_cast<unsigned>(nTS)){
    std::cout << " -- problem, event " << evtNum << " nTS = " << nTS << " lcs vector size = " << lcsFromTrkster.size() << std::endl;
    return;
  }
  for (int its = 0; its < nTS; ++its) {
    nLC += lcsFromTrkster[its].size();
  }
  lc_idx.reserve(nLC);
  lc_TSidx.reserve(nLC);
  lc_energy.reserve(nLC);
  lc_eta.reserve(nLC);
  lc_phi.reserve(nLC);
  lc_seedEnergy.reserve(nLC);
  lc_seedEta.reserve(nLC);
  lc_seedPhi.reserve(nLC);
  lc_x.reserve(nLC);
  lc_y.reserve(nLC);
  lc_z.reserve(nLC);
  lc_algo.reserve(nLC);
  lc_layer.reserve(nLC);
  lc_isSi.reserve(nLC);
  lc_nrechits.reserve(nLC);
  lc_tsMult.reserve(nLC);
  lc_mult.reserve(nLayers*nTS);
  lc_nSC.reserve(nLC);
  lc_pdgid.reserve(nLC);
  for (unsigned isc(0); isc<5; ++isc){//loop on layers   
    lc_SCidx[isc].reserve(nLC);
    lc_SCefrac[isc].reserve(nLC);
  }
  
  //std::map<unsigned,double>checkMult;
  //std::map<int,int>lMapCheckMult;
  //lMapCheckMult.clear();

  for (int its = 0; its < nTS; ++its) {
    ts_seedIdx.push_back(tracksters[its].seedIndex());
    ts_emEnergy.push_back(tracksters[its].raw_em_energy());
    ts_energy.push_back(tracksters[its].raw_energy());
    ts_regEnergy.push_back(tracksters[its].regressed_energy());
    ts_sigma1.push_back(tracksters[its].sigmasPCA()[0]);
    ts_sigma2.push_back(tracksters[its].sigmasPCA()[1]);
    ts_sigma3.push_back(tracksters[its].sigmasPCA()[2]);
    ts_BCx.push_back(tracksters[its].barycenter().x());
    ts_BCy.push_back(tracksters[its].barycenter().y());
    ts_BCz.push_back(tracksters[its].barycenter().z());
    ts_eta_PCA.push_back(tracksters[its].eigenvectors()[0].eta());
    ts_phi_PCA.push_back(tracksters[its].eigenvectors()[0].phi());
    //ts_outInHopsPerformed.push_back(tracksters[its].outInHopsPerformed());
    ts_photon_proba.push_back(tracksters[its].id_probability(ticl::Trackster::ParticleType::photon));
    ts_ele_proba.push_back(tracksters[its].id_probability(ticl::Trackster::ParticleType::electron));
    ts_mu_proba.push_back(tracksters[its].id_probability(ticl::Trackster::ParticleType::muon));
    ts_pi0_proba.push_back(tracksters[its].id_probability(ticl::Trackster::ParticleType:: neutral_pion));
    ts_chHad_proba.push_back(tracksters[its].id_probability(ticl::Trackster::ParticleType::charged_hadron));
    ts_neHad_proba.push_back(tracksters[its].id_probability(ticl::Trackster::ParticleType::neutral_hadron));
    ts_ambg_proba.push_back(tracksters[its].id_probability(ticl::Trackster::ParticleType::ambiguous));
    ts_unkwn_proba.push_back(tracksters[its].id_probability(ticl::Trackster::ParticleType::unknown));
    
    std::map<int,int>lMapLC;
    lMapLC.clear();

    if (mDebug>1) {
      std::cout << " -- trackster " << its ;
      std::cout << ", num LC = " << lcsFromTrkster[its].size() << std::endl;
    }
    
    if (lcsFromTrkster[its].size()!=0) {
      ts_eta_fromLC.push_back(lcsFromTrkster[its][0].eta_);
      ts_phi_fromLC.push_back(lcsFromTrkster[its][0].phi_);
    }
    else {
      //std::cout << " -- Problem, event " << evtNum << " trk " << its << " has no LCs " << std::endl;
      ts_eta_fromLC.push_back(10);
      ts_phi_fromLC.push_back(10);
    }
    int firstLay=100;
    int lastLay=0;
    ts_nLC.push_back(lcsFromTrkster[its].size());
    

    unsigned lcNum = 0;
    double checkTSenergy = 0;

    for (auto const& lc : lcsFromTrkster[its]) {
      lc_idx.push_back(lc.idxTracksterLC_);
      lc_TSidx.push_back(its);
      lc_energy.push_back(lc.energy_);
      lc_eta.push_back(lc.eta_);
      lc_phi.push_back(lc.phi_);
      lc_seedEnergy.push_back(lc.seedEnergy_);
      lc_seedEta.push_back(lc.seedEta_);
      lc_seedPhi.push_back(lc.seedPhi_);
      lc_x.push_back(lc.x_);
      lc_y.push_back(lc.y_);
      lc_z.push_back(lc.z_);
      lc_algo.push_back(lc.algo_);
      lc_layer.push_back(lc.layer_);
      lc_isSi.push_back(lc.isSi_);
      std::pair<std::map<int,int>::iterator,bool> isInserted =  lMapLC.insert(std::pair<int,int>(lc.layer_,1));
      if (!isInserted.second) isInserted.first->second += 1;
      //isInserted =  lMapCheckMult.insert(std::pair<int,int>(lc.idxTracksterLC_,1));
      //if (!isInserted.second) isInserted.first->second += 1;
      
      if (lc.layer_<firstLay) firstLay=lc.layer_;
      if (lc.layer_>lastLay) lastLay=lc.layer_;
      lc_nrechits.push_back(lc.nrechits_);


      if (newSimLCmult.size()==tracksters.size()) {
	if (newSimLCmult[its].size() == lcsFromTrkster[its].size()){
	  lc_tsMult.push_back(newSimLCmult[its][lcNum]);
	  checkTSenergy+=lc.energy_/newSimLCmult[its][lcNum];
	}
	else {
	  if (lcNum==0) std::cout << " Event " << evtNum << " Wrong size of newSim vector for TS " << its << " new " << newSimLCmult[its].size() << " orig " << lcsFromTrkster[its].size() << std::endl;
	  lc_tsMult.push_back(1);
	}
	//checkMultInsert = checkMult.insert(std::pair<unsigned,double>(lc.idxTracksterLC_,1./lc_tsMult[lc_tsMult.size()-1]));
	//if (!checkMultInsert.second) checkMultInsert.first->second += 1./lc_tsMult[lc_tsMult.size()-1];
      }
      else {
	if (lcNum==0 && its==0 && newSimLCmult.size()>0) std::cout << " Event " << evtNum << " Wrong size of newSim vector ! new " <<  newSimLCmult.size() << " orig " << tracksters.size() << std::endl;
	lc_tsMult.push_back(lc.tsMult_);
	//CAMM something not right with Sim tracksters due to rounding....
	//if (iterName.find("Sim")!=iterName.npos){
	//if (lc.tsMult_ > 0.01) checkTSenergy+=lc.energy_/lc.tsMult_; 
	//}
	//else 
	checkTSenergy+=lc.energy_/lc.tsMult_;
      }


      const edm::Ref<reco::CaloClusterCollection> lcRef(layerClusterHandle,lc.idxTracksterLC_);
      const auto& scsIt = recSimColl.find(lcRef);
      if (scsIt == recSimColl.end()) { 
	lcNum++;
	lc_nSC.push_back(0);
	continue;
      }
      // loop over the SimClusters contributing to this LC 
      const auto& scs = scsIt->val;
      std::vector<std::pair<unsigned,double> > scVec;
      for (const auto& scPair : scs) {
	scVec.push_back(std::pair<unsigned,double>(scPair.first.index(),scPair.second));
      }
      std::sort(scVec.begin(),scVec.end(),sortSCbyEFrac);

      lc_nSC.push_back(scVec.size());
      if (scVec.size()>0) {
	lc_pdgid.push_back(sc_pdgid[scVec[0].first]);
      } else {	
	lc_pdgid.push_back(0);
      }
      if (mDebug>2) {
	std::cout << " LC Id in Trackster = " << lcNum 		  
		  << " , LC Id in global LC collection = " << lc.idxTracksterLC_
		  << " , LC layer = " << lc.layer_
		  << " , E(LC) = " << lc.energy_
		  << " nSimClus " << scVec.size() << " (idx,E) = ";
      }
      unsigned scIter = 0;
      for (const auto& scPair : scVec) {
	//save only up to 5...
	if (scIter<5) { 
	  if (mDebug>2) std::cout << "(" << scPair.first
				  << "," << scPair.second << ") "; 
	  lc_SCidx[scIter].push_back(scPair.first);
	  lc_SCefrac[scIter].push_back(scPair.second);
	  scIter++;
	}
	//}
      } // end of looping over the SimClusters contributing to this LC
      if (mDebug>2) std::cout << std::endl;
      
      lcNum++;
      
    }//loop over LCs

    for (unsigned iL(0); iL<nLayers;++iL){
      std::map<int,int>::iterator lEle = lMapLC.find(iL+1);
      if (lEle !=lMapLC.end()) lc_mult.push_back(lEle->second);
      else lc_mult.push_back(0);
    }
    
    lMapLC.clear();
    
    ts_firstLayer.push_back(firstLay);
    ts_lastLayer.push_back(lastLay);

    if (fabs(checkTSenergy-tracksters[its].raw_energy())>0.001){
      if (iterName.find("Sim")==iterName.npos || evtNum==1) {
	std::cout << " -- Inconsistent total energy ! Event " << evtNum
		  << " iter " << iterName << " Sum LC = " 
		  << checkTSenergy << " TS energy = "
		  << tracksters[its].raw_energy()
		  << std::endl;
      }
    }
    
  }//loop on TS

  /*
  if (nTS==1 && lMapCheckMult.size()>0){
    //std::cout << " Check map mult: " << std::endl;
    int ntot = 0;
    std::map<int,int>::iterator checkMultIter = lMapCheckMult.begin();
    for (;checkMultIter!=lMapCheckMult.end();++checkMultIter){
      if (checkMultIter->second >1. )
	ntot+=checkMultIter->second;
    }
    if (ntot != 0) std::cout << " -- event " << evtNum-1
			     << " iter " << iterName
			     << " ntot " << ntot
			     << std::endl;
    
   }
  */   
      


  
}




