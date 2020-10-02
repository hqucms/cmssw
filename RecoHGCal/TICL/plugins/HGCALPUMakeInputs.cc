#include <iostream>
#include <memory>
#include <cmath>
#include <numeric>

#include "HGCALPUMakeInputs.h"

#include "PhysicsTools/NearestNeighbors/interface/NearestNeighbors.h"

using namespace ticl;
using PointCloud = cms::nanoflann::PointCloud<float>;

std::unordered_map<std::string, std::vector<float>> HGCALPUMakeInputs::makeFeatureMap(
    const std::vector<reco::CaloClusterPtr>& layerClusters,
    const edm::ValueMap<std::pair<float, float>>& layerClusterTime,
    const hgcal::RecHitTools& recHitTools) {
  PointCloud::Points points;
  std::unordered_map<std::string, std::vector<float>> fts;

  // features
  for (unsigned idx = 0; idx < layerClusters.size() && idx < max_num_layer_clusters; ++idx) {
    const auto& lc_ptr = layerClusters[idx];
    const auto& lc = *lc_ptr;
    const auto& lc_time = layerClusterTime[lc_ptr];

    points.push_back({{float(std::abs(lc.eta()) - 2.25), float(lc.phi()), float((std::abs(lc.z()) - 344) * 0.02)}});

    fts["lc_abseta"].push_back(std::abs(lc.eta()));
    fts["lc_phi"].push_back(lc.phi());
    fts["lc_absZ"].push_back(std::abs(lc.z()));
    fts["lc_x"].push_back(lc.x());
    fts["lc_y"].push_back(lc.y());
    fts["lc_time"].push_back(std::max(lc_time.first, 0.f));
    fts["lc_timeErr"].push_back(std::max(lc_time.second, 0.f));
    fts["lc_nrechits"].push_back(lc.size());
    fts["lc_logE"].push_back(std::log(lc.energy()));
    fts["lc_layer"].push_back(recHitTools.getLayerWithOffset(lc.hitsAndFractions().at(0).first));
  }

  // indices to sort the lc in decreasing energy
  // since the given `layerClusters` is already sorted by energy
  // this is just a simple sequence of 0 to max_num_layer_clusters
  fts.emplace("lc_perm_idx", std::vector<float>(max_num_layer_clusters));
  std::iota(fts["lc_perm_idx"].begin(), fts["lc_perm_idx"].end(), 0);

  // kNN indices
  auto knn = [&](size_t num_neighbors, size_t max_support_size, size_t max_query_size = 0) {
    if (max_query_size == 0)
      max_query_size = max_support_size;

    PointCloud::Points supports(points.begin(), points.begin() + std::min(max_support_size, points.size()));
    PointCloud::Points queries(points.begin(), points.begin() + std::min(max_query_size, points.size()));
    auto result = PointCloud::knn<float>(supports, queries, num_neighbors);  // queries.size() * num_neighbors
    std::vector<float> output(max_query_size * num_neighbors, max_support_size - 1);
    std::copy(result.begin(), result.begin() + std::min(result.size(), output.size()), output.begin());

    std::cout << "[knn] k=" << num_neighbors << ", num_points=" << points.size()
              << ", max_support_size=" << max_support_size << ", max_query_size=" << max_query_size
              << ", knn_result_size=" << result.size() << ", final_output_size=" << output.size() << std::endl;
    return output;
  };

  fts["encode_knnIdx1"] = knn(16, max_num_layer_clusters);
  fts["encode_knnIdx2"] = knn(16, max_num_layer_clusters / 4);
  fts["encode_knnIdx3"] = knn(16, max_num_layer_clusters / 16);
  fts["encode_knnIdx4"] = knn(16, max_num_layer_clusters / 64);
  fts["decode_knnIdx1"] = knn(1, max_num_layer_clusters / 256, max_num_layer_clusters / 64);
  fts["decode_knnIdx2"] = knn(1, max_num_layer_clusters / 64, max_num_layer_clusters / 16);
  fts["decode_knnIdx3"] = knn(1, max_num_layer_clusters / 16, max_num_layer_clusters / 4);
  fts["decode_knnIdx4"] = knn(1, max_num_layer_clusters / 4, max_num_layer_clusters);

  return fts;
}
