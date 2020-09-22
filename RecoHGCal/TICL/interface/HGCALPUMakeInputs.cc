#include "HGCALPUMakeInputs.h"

HGCALPUMakeInputs::HGCALPUMakeInputs()

HGCALPUMakeInputs::~HGCALPUMakeInputs(){};

void HGCALPUMakeInputs::makeFeatureMap(const std::vector<reco::CaloCluster>& layerClusters) {

  std::vector<float> lc_x;
  std::vector<float> lc_y;
  std::vector<float> lc_z;
  std::vector<float> lc_eta;
  std::vector<float> lc_phi;
  std::vector<float> lc_energy;
  std::vector<float> lc_nrechits;
  std::vector<float> lc_t;
  std::vector<float> lc_tErr;
  std::vector<int>   lc_layer;

  for (const auto& lc_ptr : layerClusters.ptrs()) {

    const auto& lc = *lc_ptr;
    const auto& lc_time = layerClusterTimeMap[lc_ptr];

    lc_x.push_back(lc.x());
    lc_y.push_back(lc.y());
    lc_z.push_back(lc.z());
    lc_eta.push_back(lc.eta());
    lc_phi.push_back(lc.phi());
    lc_energy.push_back(lc.energy());
    lc_nrechits.push_back(lc.size());
    lc_t.push_back(lc_time.first);
    lc_tErr.push_back(c_time.second);


    int layer_id = -1;
    float hits_energy_from_cp = 0;
    for (const auto& rechit : lc.hitsAndFractions()) {
      auto detid = rechit.first;
      if (layer_id < 0) {
	layer_id = recHitTools->getLayerWithOffset(detid);
      }
      auto rh_energy = hitMap.at(detid)->energy();
      // find energy fraction from CaloParticle
      auto it = cp_hits_and_fractions.find(detid);
      if (it != cp_hits_and_fractions.end()) {
	hits_energy_from_cp += it->second * rh_energy;
      }
    }
    lc_layer.push_back(layer_id);

  }

  // naming from: https://gitlab.cern.ch/hqu/hgcal_reco_analysis/-/blob/dev/pileup/HGCMLAnalyzer/plugins/PileupMitigationGlobalAnalyzer.cc
  feature_map_["lc_x"]        = lc_x;
  feature_map_["lc_y"]        = lc_y;
  feature_map_["lc_z"]        = lc_z;
  feature_map_["lc_eta"]      = lc_eta;
  feature_map_["lc_phi"]      = lc_phi;
  feature_map_["lc_energy"]   = lc_energy;
  feature_map_["lc_nrechits"] = lc_nrechits;
  feature_map_["lc_t"]        = lc_t;
  feature_map_["lc_tErr"]     = lc_tErr;
  feature_map_["lc_layer"]    = lc_layer;

}

