#include "ICTICLAnalysis/TiCLTreeProducer/interface/LightTree.h"

LightTree::LightTree(){
  nL=28;
};
LightTree::~LightTree(){
  lc_energy.clear();
  lc_eta.clear();
  lc_phi.clear();
  lc_x.clear();
  lc_y.clear();
  lc_z.clear();
  lc_algo.clear();
  lc_layer.clear();
  lc_nrechits.clear();
  lc_tsMult.clear();
  lc_mult.clear();
  
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
			 const int aFillTripletsInfo){


  outputTree = aFile->make<TTree>("TSTree", ("tree"+aIterName).c_str());
    
  outputTree->Branch("run", &run, "run/I");
  outputTree->Branch("event", &event, "event/I");
  outputTree->Branch("lumi", &lumi, "lumi/I");
  outputTree->Branch("weight", &weight, "weight/D");
  outputTree->Branch("nTS", &nTS, "nTS/I");
  outputTree->Branch("trackster", &trackster, "trackster/I");
  outputTree->Branch("ts_energy", &ts_energy, "ts_energy/D");
  outputTree->Branch("ts_emEnergy", &ts_emEnergy, "ts_emEnergy/D");
  outputTree->Branch("ts_regEnergy", &ts_regEnergy, "ts_regEnergy/D");
  outputTree->Branch("ts_sigma1", &ts_sigma1, "ts_sigma1/D");
  outputTree->Branch("ts_sigma2", &ts_sigma2, "ts_sigma2/D");
  outputTree->Branch("ts_sigma3", &ts_sigma3, "ts_sigma3/D");
  outputTree->Branch("ts_BCx", &ts_BCx, "ts_BCx/D");
  outputTree->Branch("ts_BCy", &ts_BCy, "ts_BCy/D");
  outputTree->Branch("ts_BCz", &ts_BCz, "ts_BCz/D");
  outputTree->Branch("ts_eta_PCA", &ts_eta_PCA, "ts_eta_PCA/D");
  outputTree->Branch("ts_phi_PCA", &ts_phi_PCA, "ts_phi_PCA/D");
  outputTree->Branch("ts_eta_fromLC", &ts_eta_fromLC, "ts_eta_fromLC/D");
  outputTree->Branch("ts_phi_fromLC", &ts_phi_fromLC, "ts_phi_fromLC/D");
  outputTree->Branch("ts_photon_proba", &ts_photon_proba, "ts_photon_proba/D");
  outputTree->Branch("ts_ele_proba", &ts_ele_proba, "ts_ele_proba/D");
  outputTree->Branch("ts_mu_proba", &ts_mu_proba, "ts_mu_proba/D");
  outputTree->Branch("ts_pi0_proba", &ts_pi0_proba, "ts_pi0_proba/D");
  outputTree->Branch("ts_chHad_proba", &ts_chHad_proba, "ts_chHad_proba/D");
  outputTree->Branch("ts_neHad_proba", &ts_neHad_proba, "ts_neHad_proba/D");
  outputTree->Branch("ts_ambg_proba", &ts_ambg_proba, "ts_ambg_proba/D");
  outputTree->Branch("ts_unkwn_proba", &ts_unkwn_proba, "ts_unkwn_proba/D");
  outputTree->Branch("ts_firstLayer", &ts_firstLayer, "ts_firstLayer/I");
  outputTree->Branch("ts_lastLayer", &ts_lastLayer, "ts_lastLayer/I");
  outputTree->Branch("ts_outInHopsPerformed", &ts_outInHopsPerformed, "ts_outInHopsPerformed/I");
  
  outputTree->Branch("nCP", &nCP, "nCP/I");
  outputTree->Branch("cp_missingEnergyFraction", &cp_missingEnergyFraction, "cp_missingEnergyFraction/D");
  outputTree->Branch("cp_energy", &cp_energy, "cp_energy/D");
  outputTree->Branch("cp_pt", &cp_pt, "cp_pt/D");
  outputTree->Branch("cp_eta", &cp_eta, "cp_eta/D");
  outputTree->Branch("cp_phi", &cp_phi, "cp_phi/D");
  outputTree->Branch("cp_pdgid", &cp_pdgid, "cp_pdgid/I");
  
  outputTree->Branch("nLC", &nLC, "nLC/I");
  outputTree->Branch("lc_energy", &lc_energy);
  outputTree->Branch("lc_eta", &lc_eta);
  outputTree->Branch("lc_phi", &lc_phi);
  outputTree->Branch("lc_x", &lc_x);
  outputTree->Branch("lc_y", &lc_y);
  outputTree->Branch("lc_z", &lc_z);
  outputTree->Branch("lc_algo", &lc_algo);
  outputTree->Branch("lc_layer", &lc_layer);
  outputTree->Branch("lc_nrechits", &lc_nrechits);
  outputTree->Branch("lc_tsMult", &lc_tsMult);
  outputTree->Branch("lc_mult", &lc_mult);
  
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
  trackster = -1;
  ts_energy = 0;
  ts_emEnergy = 0;
  ts_regEnergy = 0;
  ts_sigma1 = 0;
  ts_sigma2 = 0;
  ts_sigma3 = 0;
  ts_BCx = 0;
  ts_BCy = 0;
  ts_BCz = 0;
    ts_eta_PCA = 0;
    ts_phi_PCA = 0;
    ts_eta_fromLC = 0;
    ts_phi_fromLC = 0;
    ts_photon_proba = 0;
    ts_ele_proba = 0;
    ts_mu_proba = 0;
    ts_pi0_proba = 0;
    ts_chHad_proba = 0;
    ts_neHad_proba = 0;
    ts_ambg_proba = 0;
    ts_unkwn_proba = 0;
    ts_firstLayer = 0;
    ts_lastLayer = 0;
    ts_outInHopsPerformed = 0;
    
    nCP = 0;
    cp_missingEnergyFraction = 0;
    cp_energy = 0;
    cp_pdgid = 0;
    cp_pt = 0;
    cp_eta = 0;
    cp_phi = 0;
    
    nLC = 0;
    lc_energy.clear();
    lc_eta.clear();
    lc_phi.clear();
    lc_x.clear();
    lc_y.clear();
    lc_z.clear();
    lc_algo.clear();
    lc_layer.clear();
    lc_nrechits.clear();
    lc_tsMult.clear();
    lc_mult.clear();

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
		  const int icp){
  nCP = caloparticles.size();
  if (icp>=0){
    cp_energy = caloparticles[icp].energy_;
    cp_pt = caloparticles[icp].pt_;
    cp_eta = caloparticles[icp].eta_;
    cp_phi = caloparticles[icp].phi_;
    cp_pdgid = caloparticles[icp].pdgid_;
  }
}



void LightTree::fillTSinfo(const std::vector<ticl::Trackster> & tracksters,
			   const int itrksterMin,
			   const std::vector<layercluster> & lcsFromClosestTrksterToCP){
  
    nTS = tracksters.size();
    std::map<int,int>lMapLC;
    
    trackster = itrksterMin;
    if (itrksterMin >= 0) { 
	ts_emEnergy = tracksters[itrksterMin].raw_em_energy();
	ts_energy = tracksters[itrksterMin].raw_energy();
	ts_regEnergy = tracksters[itrksterMin].regressed_energy();
	ts_sigma1 = tracksters[itrksterMin].sigmasPCA()[0];
	ts_sigma2 = tracksters[itrksterMin].sigmasPCA()[1];
	ts_sigma3 = tracksters[itrksterMin].sigmasPCA()[2];
	ts_BCx = tracksters[itrksterMin].barycenter().x();
	ts_BCy = tracksters[itrksterMin].barycenter().y();
	ts_BCz = tracksters[itrksterMin].barycenter().z();
	ts_eta_PCA = tracksters[itrksterMin].eigenvectors()[0].eta();
	ts_phi_PCA = tracksters[itrksterMin].eigenvectors()[0].phi();
	ts_outInHopsPerformed = tracksters[itrksterMin].outInHopsPerformed();
	ts_photon_proba = tracksters[itrksterMin].id_probability(ticl::Trackster::ParticleType::photon);
	ts_ele_proba = tracksters[itrksterMin].id_probability(ticl::Trackster::ParticleType::electron);
	ts_mu_proba = tracksters[itrksterMin].id_probability(ticl::Trackster::ParticleType::muon);
	ts_pi0_proba = tracksters[itrksterMin].id_probability(ticl::Trackster::ParticleType:: neutral_pion);
	ts_chHad_proba = tracksters[itrksterMin].id_probability(ticl::Trackster::ParticleType::charged_hadron);
	ts_neHad_proba = tracksters[itrksterMin].id_probability(ticl::Trackster::ParticleType::neutral_hadron);
	ts_ambg_proba = tracksters[itrksterMin].id_probability(ticl::Trackster::ParticleType::ambiguous);
	ts_unkwn_proba = tracksters[itrksterMin].id_probability(ticl::Trackster::ParticleType::unknown);
	
	lMapLC.clear();
      
	ts_eta_fromLC = lcsFromClosestTrksterToCP[0].eta_;
	ts_phi_fromLC = lcsFromClosestTrksterToCP[0].phi_;
	
	int firstLay=30;
	int lastLay=0;
	nLC = lcsFromClosestTrksterToCP.size();
	for (auto const& lc : lcsFromClosestTrksterToCP) {
	  lc_energy.push_back(lc.energy_);
	  lc_eta.push_back(lc.eta_);
	  lc_phi.push_back(lc.phi_);
	  lc_x.push_back(lc.x_);
	  lc_y.push_back(lc.y_);
	  lc_z.push_back(lc.z_);
	  lc_algo.push_back(lc.algo_);
	  lc_layer.push_back(lc.layer_);
	  std::pair<std::map<int,int>::iterator,bool> isInserted =  lMapLC.insert(std::pair<int,int>(lc.layer_,1));
	  if (!isInserted.second) isInserted.first->second += 1;
	  
	  if (lc.layer_<firstLay) firstLay=lc.layer_;
	  if (lc.layer_>lastLay) lastLay=lc.layer_;
	  lc_nrechits.push_back(lc.nrechits_);
	  lc_tsMult.push_back(lc.tsMult_);
	}
	
	for (unsigned iL(0); iL<28;++iL){
	  std::map<int,int>::iterator lEle = lMapLC.find(iL+1);
	  if (lEle !=lMapLC.end()) lc_mult.push_back(lEle->second);
	  else lc_mult.push_back(0);
	}
	
	lMapLC.clear();
	
	ts_firstLayer = firstLay;
	ts_lastLayer = lastLay;
	
    }//if TS found

}




