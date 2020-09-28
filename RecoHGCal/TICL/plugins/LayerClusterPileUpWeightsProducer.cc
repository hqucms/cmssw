/*
Inputs: 
 - layerClusters
 - timing & error map

Output:
 - std::vector<float> : weight for each LayerCluster

Steps:
 - make_inputs:
   - all LCs -> LCs +z, LCs -z
   - make_inputs_endcap(vector<LC*>)
     - sorted by energy
     - feature arrays 
       - obj -> array
       - center_norm_pad
     - coordinates array 
       - obj -> array
       - kNN
   - predict()
*/

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

// FIXME
//#include "DataFormats/BTauReco/interface/DeepBoostedJetFeatures.h"
#include "RecoHGCal/TICL/interface/HGCALPUMakeInputs.h"
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <nlohmann/json.hpp>

using namespace cms::Ort;

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
  //  typedef std::vector<reco::DeepBoostedJetTagInfo> TagInfoCollection;
  //  typedef reco::JetTagCollection JetTagCollection;

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
  //void make_inputs(const reco::DeepBoostedJetTagInfo &taginfo); // from HQ
  void make_inputs(const std::unordered_map<std::string, std::vector<float>>  &taginfo);

  const edm::EDGetTokenT<TagInfoCollection> src_;
  std::vector<std::string> output_names_;           // names of the output scores
  std::vector<std::string> input_names_;            // names of each input group - the ordering is important!
  std::vector<std::vector<int64_t>> input_shapes_;  // shapes of each input group (-1 for dynamic axis)
  std::vector<unsigned> input_sizes_;               // total length of each input vector
  std::unordered_map<std::string, PreprocessParams> prep_info_map_;  // preprocessing info for each input group

  FloatArrays data_;

  bool debug_ = false;
};

LayerClusterPileUpWeightsProducer::LayerClusterPileUpWeightsProducer(const edm::ParameterSet &iConfig, const ONNXRuntime *cache)
    : src_(consumes<TagInfoCollection>(iConfig.getParameter<edm::InputTag>("src"))),
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

  // get output names from output_names
  for (const auto &flav_name : output_names_) {
    produces<JetTagCollection>(flav_name);
  }
}

LayerClusterPileUpWeightsProducer::~LayerClusterPileUpWeightsProducer() {}

void LayerClusterPileUpWeightsProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  // layerClusterPileUpWeightsProducer
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("pfDeepBoostedJetTagInfos"));
  desc.add<edm::FileInPath>("preprocess_json",
                            edm::FileInPath("RecoBTag/Combined/data/DeepBoostedJet/V02/full/resnet.onnx"));
  desc.add<edm::FileInPath>("model_path",
                            edm::FileInPath("RecoBTag/Combined/data/DeepBoostedJet/V02/full/resnet.onnx"));
  desc.add<std::vector<std::string>>("output_names", std::vector<std::string>{"softmax"});
  desc.addOptionalUntracked<bool>("debugMode", false);

  descriptions.addWithDefaultLabel(desc);
}

std::unique_ptr<ONNXRuntime> LayerClusterPileUpWeightsProducer::initializeGlobalCache(const edm::ParameterSet &iConfig) {
  return std::make_unique<ONNXRuntime>(iConfig.getParameter<edm::FileInPath>("model_path").fullPath());
}

void LayerClusterPileUpWeightsProducer::globalEndJob(const ONNXRuntime *cache) {}

void LayerClusterPileUpWeightsProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  edm::Handle<TagInfoCollection> tag_infos;
  iEvent.getByToken(src_, tag_infos);

  // get the unordered map:
  std::unordered_map<std::string, std::vector<float>> inputs_ = makeFeatureMap(const std::vector<reco::CaloCluster>& layerClusters);

  // initialize output collection
  //  std::vector<std::unique_ptr<JetTagCollection>> output_tags;
  std::vector<std::unique_ptr<std::vector<reco::CaloCluster>> output_tags;
  
  
  // what should I do with this part?
  if (!tag_infos->empty()) {
    auto jet_ref = tag_infos->begin()->jet();
    auto ref2prod = edm::makeRefToBaseProdFrom(jet_ref, iEvent);
    for (std::size_t i = 0; i < output_names_.size(); i++) {
      output_tags.emplace_back(std::make_unique<JetTagCollection>(ref2prod));
    }
  } else {
    for (std::size_t i = 0; i < output_names_.size(); i++) {
      output_tags.emplace_back(std::make_unique<JetTagCollection>());
    }
  }


  for (unsigned lc_n = 0; lc_n < tag_infos->size(); ++lc_n) {
    const auto &taginfo = (*tag_infos)[lc_n];
    std::vector<float> outputs(output_names_.size(), 0);  // init as all zeros

    if (!taginfo.features().empty()) {
      // convert inputs
      make_inputs(taginfo);
      // run prediction and get outputs
      outputs = globalCache()->run(input_names_, data_, input_shapes_)[0];
      assert(outputs.size() == output_names_.size());
    }

    const auto &lc_ref = tag_infos->at(lc_n).jet();
    for (std::size_t flav_n = 0; flav_n < output_names_.size(); flav_n++) {
      (*(output_tags[flav_n]))[lc_ref] = outputs[flav_n];
    }
  }

  /*
  // from HQ
  for (unsigned jet_n = 0; jet_n < tag_infos->size(); ++jet_n) {
    const auto &taginfo = (*tag_infos)[jet_n];
    std::vector<float> outputs(output_names_.size(), 0);  // init as all zeros

    if (!taginfo.features().empty()) {
      // convert inputs
      make_inputs(taginfo);
      // run prediction and get outputs
      outputs = globalCache()->run(input_names_, data_, input_shapes_)[0];
      assert(outputs.size() == output_names_.size());
    }

    const auto &jet_ref = tag_infos->at(jet_n).jet();
    for (std::size_t flav_n = 0; flav_n < output_names_.size(); flav_n++) {
      (*(output_tags[flav_n]))[jet_ref] = outputs[flav_n];
    }
  }
  */
  if (debug_) {
    std::cout << "=== " << iEvent.id().run() << ":" << iEvent.id().luminosityBlock() << ":" << iEvent.id().event()
              << " ===" << std::endl;
    for (unsigned jet_n = 0; jet_n < tag_infos->size(); ++jet_n) {
      const auto &jet_ref = tag_infos->at(jet_n).jet();
      std::cout << " - Jet #" << jet_n << ", pt=" << jet_ref->pt() << ", eta=" << jet_ref->eta()
                << ", phi=" << jet_ref->phi() << std::endl;
      for (std::size_t flav_n = 0; flav_n < output_names_.size(); ++flav_n) {
        std::cout << "    " << output_names_.at(flav_n) << " = " << (*(output_tags.at(flav_n)))[jet_ref] << std::endl;
      }
    }
  }

  // put into the event
  for (std::size_t flav_n = 0; flav_n < output_names_.size(); ++flav_n) {
    iEvent.put(std::move(output_tags[flav_n]), output_names_[flav_n]);
  }
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
    out[i] = std::clamp((catch_infs(input[i], replace_inf_value) - center) * norm_factor, min, max);
  }
  return out;
}

//void LayerClusterPileUpWeightsProducer::make_inputs(const reco::DeepBoostedJetTagInfo &taginfo) {
void LayerClusterPileUpWeightsProducer::make_inputs(std::unordered_map<std::string, std::vector<float>> &taginfo) {
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
      const auto &raw_value = taginfo.features().get(varname);
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
