#ifndef RecoHGCal_TICL_HGCALPUMakeFeatures_h
#define RecoHGCal_TICL_HGCALPUMakeFeatures_h

#include <string>
#include <vector>
#include <unordered_map>

#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

namespace ticl {

  class HGCALPUMakeInputs {
  public:
    HGCALPUMakeInputs() {}
    ~HGCALPUMakeInputs() {}

    std::unordered_map<std::string, std::vector<float>> makeFeatureMap(
        const std::vector<edm::Ptr<reco::CaloCluster>>& layerClusters,
        const edm::ValueMap<std::pair<float, float>>& layerClusterTime,
        const std::shared_ptr<hgcal::RecHitTools> recHitTools);

    constexpr static unsigned max_num_layer_clusters = 80000;
  };
}  // namespace ticl

#endif  // RecoHGCal_TICL_HGCALPUMakeFeatures_h
