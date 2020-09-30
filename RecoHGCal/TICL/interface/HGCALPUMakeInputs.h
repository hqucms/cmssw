#ifndef RecoHGCal_TICL_HGCALPUMakeFeatures_h
#define RecoHGCal_TICL_HGCALPUMakeFeatures_h

#include <string>
#include <vector>

#include <unordered_map>
#include <memory>
#include <algorithm>

#include <math.h>

#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"


class HGCALPUMakeInputs {

 public:						       
  HGCALPUMakeInputs();
  ~HGCALPUMakeInputs() override;
 
  std::unordered_map<std::string, std::vector<float>> makeFeatureMap(const std::vector<edm::Ptr<reco::CaloCluster>>& layerClusters);

}
#endif
