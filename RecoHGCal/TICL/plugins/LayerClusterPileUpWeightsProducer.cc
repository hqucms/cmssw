#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "HGCALPUMakeInputs.h"
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <memory>
#include <nlohmann/json.hpp>

using namespace cms::Ort;
using namespace ticl;

// struct to hold preprocessing parameters
struct PreprocessParams {
  struct VarInfo {
    VarInfo() {}
    VarInfo(float median, float norm_factor, float replace_inf_value, float lower_bound, float upper_bound, float pad)
        : center(median),
          norm_factor(norm_factor),
          replace_inf_value(replace_inf_value),
          lower_bound(lower_bound),
          upper_bound(upper_bound),
          pad(pad) {}
    float center = 0;
    float norm_factor = 1;
    float replace_inf_value = 0;
    float lower_bound = -5;
    float upper_bound = 5;
    float pad = 0;
  };

  unsigned min_length = 0;
  unsigned max_length = 0;
  std::vector<std::string> var_names;
  std::unordered_map<std::string, VarInfo> var_info_map;

  VarInfo info(const std::string &name) const { return var_info_map.at(name); }
};

class LayerClusterPileUpWeightsProducer : public edm::stream::EDProducer<edm::GlobalCache<ONNXRuntime>> {
public:
  explicit LayerClusterPileUpWeightsProducer(const edm::ParameterSet &, const ONNXRuntime *);
  ~LayerClusterPileUpWeightsProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &);

  static std::unique_ptr<ONNXRuntime> initializeGlobalCache(const edm::ParameterSet &);
  static void globalEndJob(const ONNXRuntime *);

private:
  void produce(edm::Event &, const edm::EventSetup &) override;

  std::vector<float> center_norm_pad(const std::vector<float> &input,
                                     float center,
                                     float scale,
                                     unsigned min_length,
                                     unsigned max_length,
                                     float pad_value = 0,
                                     float replace_inf_value = 0,
                                     float min = 0,
                                     float max = -1);
  void make_inputs(const std::unordered_map<std::string, std::vector<float>> &taginfo);

  edm::EDGetTokenT<reco::CaloClusterView> layerClustersToken_;
  edm::EDGetTokenT<edm::ValueMap<std::pair<float, float>>> layerClusterTimeToken_;
  hgcal::RecHitTools recHitTools_;

  std::vector<std::string> output_names_;           // names of the output scores
  std::vector<std::string> input_names_;            // names of each input group - the ordering is important!
  std::vector<std::vector<int64_t>> input_shapes_;  // shapes of each input group (-1 for dynamic axis)
  std::vector<unsigned> input_sizes_;               // total length of each input vector
  std::unordered_map<std::string, PreprocessParams> prep_info_map_;  // preprocessing info for each input group

  FloatArrays data_;

  bool debug_ = false;
};

LayerClusterPileUpWeightsProducer::LayerClusterPileUpWeightsProducer(const edm::ParameterSet &iConfig,
                                                                     const ONNXRuntime *cache)
    : layerClustersToken_(consumes<reco::CaloClusterView>(iConfig.getParameter<edm::InputTag>("layerClusters"))),
      layerClusterTimeToken_(
          consumes<edm::ValueMap<std::pair<float, float>>>(iConfig.getParameter<edm::InputTag>("layerClusterTime"))),
      output_names_(iConfig.getParameter<std::vector<std::string>>("output_names")),
      debug_(iConfig.getUntrackedParameter<bool>("debugMode", false)) {
  // load preprocessing info
  std::ifstream ifs(iConfig.getParameter<edm::FileInPath>("preprocess_json").fullPath());
  nlohmann::json js = nlohmann::json::parse(ifs);
  js.at("input_names").get_to(input_names_);
  for (const auto &group_name : input_names_) {
    const auto &group_pset = js.at(group_name);
    auto &prep_params = prep_info_map_[group_name];
    group_pset.at("var_names").get_to(prep_params.var_names);
    if (group_pset.contains("var_length")) {
      prep_params.min_length = group_pset.at("var_length");
      prep_params.max_length = prep_params.min_length;
    } else {
      prep_params.min_length = group_pset.at("min_length");
      prep_params.max_length = group_pset.at("max_length");
      input_shapes_.push_back({1, (int64_t)prep_params.var_names.size(), -1});
    }
    const auto &var_info_pset = group_pset.at("var_infos");
    for (const auto &var_name : prep_params.var_names) {
      const auto &var_pset = var_info_pset.at(var_name);
      double median = var_pset.at("median");
      double norm_factor = var_pset.at("norm_factor");
      double replace_inf_value = var_pset.at("replace_inf_value");
      double lower_bound = var_pset.at("lower_bound");
      double upper_bound = var_pset.at("upper_bound");
      double pad = var_pset.contains("pad") ? double(var_pset.at("pad")) : 0;
      prep_params.var_info_map[var_name] =
          PreprocessParams::VarInfo(median, norm_factor, replace_inf_value, lower_bound, upper_bound, pad);
    }

    // create data storage with a fixed size vector initilized w/ 0
    const auto &len = input_sizes_.emplace_back(prep_params.max_length * prep_params.var_names.size());
    data_.emplace_back(len, 0);
  }

  if (debug_) {
    for (unsigned i = 0; i < input_names_.size(); ++i) {
      const auto &group_name = input_names_.at(i);
      if (!input_shapes_.empty()) {
        std::cout << group_name << "\nshapes: ";
        for (const auto &x : input_shapes_.at(i)) {
          std::cout << x << ", ";
        }
      }
      std::cout << "\nvariables: ";
      for (const auto &x : prep_info_map_.at(group_name).var_names) {
        std::cout << x << ", ";
      }
      std::cout << "\n";
    }
    std::cout << "output_names: ";
    for (const auto &flav_name : output_names_) {
      std::cout << flav_name << ", ";
    }
    std::cout << "\n";
  }

  produces<std::vector<float>>();
}

LayerClusterPileUpWeightsProducer::~LayerClusterPileUpWeightsProducer() {}

void LayerClusterPileUpWeightsProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  // layerClusterPileUpWeightsProducer
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("layerClusters", edm::InputTag("hgcalLayerClusters", "", "RECO"));
  desc.add<edm::InputTag>("layerClusterTime", edm::InputTag("hgcalLayerClusters", "timeLayerCluster", "RECO"));
  desc.add<edm::FileInPath>(
      "preprocess_json",
      edm::FileInPath("RecoHGCal/TICL/data/LayerClusterPUId/preprocess.json"));  //preproccesses json
  desc.add<edm::FileInPath>("model_path",
                            edm::FileInPath("RecoHGCal/TICL/data/LayerClusterPUId/RandLANet.onnx"));  // model
  desc.add<std::vector<std::string>>("output_names", std::vector<std::string>{"softmax"});
  desc.addOptionalUntracked<bool>("debugMode", false);

  descriptions.add("hgcalLayerClusterPileupWeights", desc);
}

std::unique_ptr<ONNXRuntime> LayerClusterPileUpWeightsProducer::initializeGlobalCache(const edm::ParameterSet &iConfig) {
  return std::make_unique<ONNXRuntime>(iConfig.getParameter<edm::FileInPath>("model_path").fullPath());
}

void LayerClusterPileUpWeightsProducer::globalEndJob(const ONNXRuntime *cache) {}

void LayerClusterPileUpWeightsProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  const auto &layerClusters = iEvent.get(layerClustersToken_);
  const auto &layerClusterTimeMap = iEvent.get(layerClusterTimeToken_);

  edm::ESHandle<CaloGeometry> geom;
  iSetup.get<CaloGeometryRecord>().get(geom);
  recHitTools_.setGeometry(*geom);

  auto results = std::make_unique<std::vector<float>>(layerClusters.size(), 0);

  for (float z_sign : {-1, 1}) {
    // do one endcap at a time
    std::vector<reco::CaloClusterPtr> lc_ptrs;
    for (const auto &lc : layerClusters.ptrs()) {
      if (lc->eta() * z_sign > 1) {
        lc_ptrs.push_back(lc);
      }
    }
    std::sort(lc_ptrs.begin(), lc_ptrs.end(), [&](const reco::CaloClusterPtr &a, const reco::CaloClusterPtr &b) {
      return a->energy() > b->energy();
    });  // sort by LC energy

    auto features = HGCALPUMakeInputs::makeFeatureMap(lc_ptrs, layerClusterTimeMap, recHitTools_);

    // pre-processing
    make_inputs(features);

    // run prediction and get outputs
    auto outputs = globalCache()->run(input_names_, data_, input_shapes_)[0];

    // lc[i] <==> outputs[i]
    // .key() gives the index in the original input collection
    for (unsigned i = 0; i < std::min(lc_ptrs.size(), outputs.size()); ++i) {
      (*results)[lc_ptrs[i].key()] = outputs[i];
    }

  }  // end of looping over the two Hgcal endcaps

  if (debug_) {
    std::cout << "=== " << iEvent.id().run() << ":" << iEvent.id().luminosityBlock() << ":" << iEvent.id().event()
              << " ===" << std::endl;
    for (unsigned idx = 0; idx < layerClusters.size(); ++idx) {
      const auto &lc = layerClusters[idx];
      std::cout << " - LC #" << idx << ", x=" << lc.x() << ", y=" << lc.y() << ", z=" << lc.z()
                << ", output=" << (*results)[idx] << std::endl;
    }
  }

  iEvent.put(std::move(results));
}

std::vector<float> LayerClusterPileUpWeightsProducer::center_norm_pad(const std::vector<float> &input,
                                                                      float center,
                                                                      float norm_factor,
                                                                      unsigned min_length,
                                                                      unsigned max_length,
                                                                      float pad_value,
                                                                      float replace_inf_value,
                                                                      float min,
                                                                      float max) {
  // do variable shifting/scaling/padding/clipping in one go

  assert(min <= pad_value && pad_value <= max);
  assert(min_length <= max_length);

  unsigned target_length = std::clamp((unsigned)input.size(), min_length, max_length);
  std::vector<float> out(target_length, pad_value);
  for (unsigned i = 0; i < input.size() && i < target_length; ++i) {
    auto val = std::isfinite(input[i]) ? input[i] : replace_inf_value;
    out[i] = std::clamp((val - center) * norm_factor, min, max);
  }
  return out;
}

void LayerClusterPileUpWeightsProducer::make_inputs(
    const std::unordered_map<std::string, std::vector<float>> &features) {
  auto get = [&features](const std::string &name) {
    auto item = features.find(name);
    if (item != features.end()) {
      return item->second;
    } else {
      throw cms::Exception("InvalidArgument") << "Feature " << name << " does not exist!";
    }
  };

  for (unsigned igroup = 0; igroup < input_names_.size(); ++igroup) {
    const auto &group_name = input_names_[igroup];
    const auto &prep_params = prep_info_map_.at(group_name);
    auto &group_values = data_[igroup];
    group_values.resize(input_sizes_[igroup]);
    // first reset group_values to 0
    std::fill(group_values.begin(), group_values.end(), 0);
    unsigned curr_pos = 0;
    // transform/pad
    for (unsigned i = 0; i < prep_params.var_names.size(); ++i) {
      const auto &varname = prep_params.var_names[i];
      const auto &raw_value = get(varname);
      const auto &info = prep_params.info(varname);
      auto val = center_norm_pad(raw_value,
                                 info.center,
                                 info.norm_factor,
                                 prep_params.min_length,
                                 prep_params.max_length,
                                 info.pad,
                                 info.replace_inf_value,
                                 info.lower_bound,
                                 info.upper_bound);
      std::copy(val.begin(), val.end(), group_values.begin() + curr_pos);
      curr_pos += val.size();
      if (i == 0 && (!input_shapes_.empty())) {
        input_shapes_[igroup][2] = val.size();
      }

      if (debug_) {
        std::cout << " -- var=" << varname << ", center=" << info.center << ", scale=" << info.norm_factor
                  << ", replace=" << info.replace_inf_value << ", pad=" << info.pad << std::endl;
        for (const auto &v : val) {
          std::cout << v << ",";
        }
        std::cout << std::endl;
      }
    }
    group_values.resize(curr_pos);
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(LayerClusterPileUpWeightsProducer);
