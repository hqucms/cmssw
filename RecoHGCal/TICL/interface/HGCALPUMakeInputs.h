#ifndef RecoHGCal_TICL_HGCALPUMakeFeatures_h
#define RecoHGCal_TICL_HGCALPUMakeFeatures_h

#include <string>
#include <vector>

#include <unordered_map>
#include <memory>
#include <algorithm>

#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"


class HGCALPUMakeInputs {

 public:						       
  HGCALPUMakeInputs();
  ~HGCALPUMakeInputs() override;
 
  void makeFeatureMap(const std::vector<reco::CaloCluster>& layerClusters);

 private:
  std::unordered_map<std::string, std::vector<float>> feature_map_;

}
#endif
