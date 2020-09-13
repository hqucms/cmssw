#include "HGCALPUMakeInputs.h"

HGCALPUMakeInputs::HGCALPUMakeInputs()

HGCALPUMakeInputs::~HGCALPUMakeInputs(){};

void HGCALPUMakeInputs::makeFeatureMap(const std::vector<reco::CaloCluster>& layerClusters) {

  std::vector<float> lc_energy;

  for (const auto& lc_ptr : layerClusters.ptrs()) {
    const auto& lc = *lc_ptr;
    lc_energy.push_back(lc.energy());
  }

  // naming from: https://gitlab.cern.ch/hqu/hgcal_reco_analysis/-/blob/dev/pileup/HGCMLAnalyzer/plugins/PileupMitigationGlobalAnalyzer.cc
  feature_map_["lc_energy"] = lc_energy;


}

