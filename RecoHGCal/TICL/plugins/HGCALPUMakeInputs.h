#ifndef RecoHGCal_TICL_HGCALPUMakeFeatures_h
#define RecoHGCal_TICL_HGCALPUMakeFeatures_h

#include <string>
#include <vector>
#include <unordered_map>

#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

namespace ticl {

  class HGCALPUMakeInputs {
  public:
    HGCALPUMakeInputs() {}
    ~HGCALPUMakeInputs() {}

    static std::unordered_map<std::string, std::vector<float>> makeFeatureMap(
        const std::vector<reco::CaloClusterPtr>& layerClusters,
        const edm::ValueMap<std::pair<float, float>>& layerClusterTime,
        const hgcal::RecHitTools& recHitTools,
        bool debug = false);

    constexpr static unsigned max_num_layer_clusters = 80000;
  };
}  // namespace ticl

#endif  // RecoHGCal_TICL_HGCALPUMakeFeatures_h
