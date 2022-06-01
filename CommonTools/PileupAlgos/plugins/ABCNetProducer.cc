//========================================================================================
// Class:      ABCNetProducer                     ---------------------------------------
//----------------------------------------------------------------------------------------
/**\class ABCNetProducer ABCNetProducer.cc PileupAlgos/plugins/ABCNetProducer.cc
------------------------------------------------------------------------------------------
 Description: This class produces ABCNet weights to be used to mitigate pileup  ---
 -----------------------------------------------------------------------------------------
 Implementation:                                                                       ---
     This EDProducer is meant to be used with CMSSW >= 10_6_26                         ---
*/
//========================================================================================
// Authors:  Fabio Iemmi (IHEP)                                      ---------------------
//         Created:  MON, 23 May 2022 10:00:28 GMT  ---------------------------------------
//========================================================================================

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h" // for edm::ParameterSet
#include "DataFormats/PatCandidates/interface/PackedCandidate.h" //for pat
#include "DataFormats/Candidate/interface/Candidate.h" // for reco
#include "CommonTools/PileupAlgos/interface/PuppiCandidate.h"// for puppi candidates
#include "DataFormats/Math/interface/LorentzVector.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h" // to use TensorFlow
#include "FWCore/Framework/interface/stream/EDAnalyzer.h"

#include "ABCNetMakeInputs.h"

#include <limits>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <memory>
#include <nlohmann/json.hpp>

using namespace abcnet;

// struct to hold preprocessing parameters
struct PreprocessParams {
  struct VarInfo {
    VarInfo() {}
    VarInfo(float upper_bound, float lower_bound, float norm_factor, float replace_nan_value, float pad_value)
      : upper_bound(upper_bound),
        lower_bound(lower_bound),
	norm_factor(norm_factor),
	replace_nan_value(replace_nan_value),
	pad_value(pad_value) {}
    float upper_bound = 100000;
    float lower_bound = -100000;
    float norm_factor = 1;
    float replace_nan_value = 0;
    float pad_value = 0;
  };

  std::vector<std::string> var_names;
  std::unordered_map<std::string, VarInfo> var_info_map;

  VarInfo info(const std::string &name) const { return var_info_map.at(name); }

};

struct ABCNetTFCache {
  ABCNetTFCache() : graphDef(nullptr) {}
  std::atomic<tensorflow::GraphDef*> graphDef;
};

///////////////////////////////////////////////////////////////////////////////////
// Declaring class ----------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////////

class ABCNetProducer : public edm::stream::EDProducer<edm::GlobalCache<ABCNetTFCache>> {

public:
  explicit ABCNetProducer(const edm::ParameterSet&, const ABCNetTFCache*);   
  ~ABCNetProducer();

  static void fillDescriptions(edm::ConfigurationDescriptions &);

  static std::unique_ptr<ABCNetTFCache> initializeGlobalCache(const edm::ParameterSet&);
  static void globalEndJob(const ABCNetTFCache*);

private:
  void produce(edm::Event &, const edm::EventSetup &) override;
  std::vector<float> minmax_scale(const std::vector<float> &input,
				  float upper_bound = 100000,
				  float lower_bound = 100000,
				  float norm_factor = 1,
				  float pad_value = 0,
				  float replace_nan_value = 0
				  );
  void preprocess(std::unordered_map<std::string, std::vector<float>> &featureMap, bool debug);
  // tokens
  edm::EDGetTokenT<reco::CandidateView> tokenPFCandidates_;
  std::vector<std::string> input_names_; //names of the input features. Ordering matters!
  std::unordered_map<std::string, PreprocessParams> prep_info_map_;  //preprocessing info for each input feature
  //session for TF evaluation
  tensorflow::Session* session_;
  std::string input_tensor_name_;
  std::string output_tensor_name_;
  bool debug_;
  constexpr static unsigned max_num_PFCandidates = 4000;
  
  edm::EDPutTokenT<edm::ValueMap<float>> ABCNetOut_;
};

//////////////////////////////////////////////////////////////////////////////////
// Defining methods --------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////////

// constructors
ABCNetProducer::ABCNetProducer(const edm::ParameterSet& iConfig, const ABCNetTFCache* cache):
  tokenPFCandidates_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("candName"))),
  session_(nullptr),
  input_tensor_name_(iConfig.getParameter<std::string>("input_tensor_name")),
  output_tensor_name_(iConfig.getParameter<std::string>("output_tensor_name")),
  debug_(iConfig.getParameter<bool>("debug"))
{
  //parse data from preprocessing JSON file
  std::ifstream ifs(iConfig.getParameter<edm::FileInPath>("preprocess_json").fullPath());
  nlohmann::json js = nlohmann::json::parse(ifs);
  js.at("input_names").get_to(input_names_); //store names of input features from JSON inside input_names_
  const auto &var_info_pset = js.at("var_info"); //get the var_info block from JSON
  
  for (const auto & input_name : input_names_) {
    const auto &var_pset = var_info_pset.at(input_name);
    double upper_bound = var_pset.at("upper_bound");
    double lower_bound = var_pset.at("lower_bound");
    double norm_factor = var_pset.at("norm_factor");
    double replace_nan_value = var_pset.at("replace_nan_value");
    double pad_value = var_pset.at("pad_value");
    auto &prep_params = prep_info_map_[input_name];
    prep_params.var_info_map[input_name] = PreprocessParams::VarInfo(upper_bound,
								     lower_bound,
								     norm_factor,
								     replace_nan_value,
								     pad_value
								     );
  }
  
  //create the TF session using the meta graph from the cache
  session_ = tensorflow::createSession(cache->graphDef);

  // Produce a ValueMap of floats linking each PF candidate with its ABCNet weight
  ABCNetOut_ = produces<edm::ValueMap<float>> ();

};

// destructor
ABCNetProducer::~ABCNetProducer() {
};

void ABCNetProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {};

std::unique_ptr<ABCNetTFCache> ABCNetProducer::initializeGlobalCache(const edm::ParameterSet& iConfig) {
  tensorflow::setLogging("3");
  // get the pb file
  std::string pbFile = iConfig.getParameter<edm::FileInPath>("graph_path").fullPath();
  // load the graph definition and save it in the cache
  ABCNetTFCache* cache = new ABCNetTFCache();
  cache->graphDef = tensorflow::loadGraphDef(pbFile);
  return std::unique_ptr<ABCNetTFCache>(cache);
};

void ABCNetProducer::globalEndJob(const ABCNetTFCache* cache) {
  if (cache->graphDef != nullptr) {
    delete cache->graphDef;
  }
};

std::vector<float> ABCNetProducer::minmax_scale(const std::vector<float> &input,
						float upper_bound,
						float lower_bound,
						float norm_factor,
						float pad_value,
						float replace_nan_value
						) {
  unsigned target_length = std::clamp((unsigned)input.size(), max_num_PFCandidates, max_num_PFCandidates);
  std::vector<float> out(target_length, pad_value);
  for (unsigned i = 0; i < input.size() && i < target_length; ++i) {
    //first of all, check if input is inf/nan and replace it if needed
    auto val = std::isfinite(input[i]) ? input[i] : replace_nan_value;
    val = std::isnan(input[i]) ? replace_nan_value : input[i];
    //second of all, perform the actual min-max scaling
    val = (val>upper_bound) ? upper_bound : val;
    val = (val<lower_bound) ? lower_bound : val;
    out[i] = val/norm_factor; //need to check this value
  }
  return out;
};

void ABCNetProducer::preprocess(std::unordered_map<std::string, std::vector<float>> & featureMap, bool debug) {

  for (const auto & input_name : input_names_) {
    //std::cout << input_name << std::endl;
    const auto & preprocessing_params = prep_info_map_[input_name].var_info_map[input_name];
    //std::cout << preprocessing_params.upper_bound << std::endl;
    featureMap[input_name] = minmax_scale(featureMap[input_name],
					  preprocessing_params.upper_bound,
					  preprocessing_params.lower_bound,
					  preprocessing_params.norm_factor,
					  preprocessing_params.pad_value,
					  preprocessing_params.replace_nan_value
					  );
  }
  
  if(debug) {
    std::cout << "*** NOW CHECKING THE SANITY OF THE PREPROCESSING STEP ***" << std::endl;
    for(auto const & feat : featureMap) {
      bool issue = false;
      std::cout << feat.first << std::endl;
      if (feat.first == "PFCandEta" || feat.first == "PFCandPhi" || feat.first == "PFCandLogPt" || feat.first == "PFCandLogE" || feat.first == "PFCandCharge") {
	std::cout << "This feature is expected to take values outside [-1;1] interval. Skipping..." << std::endl;
	continue;
      } 
      for (auto const & element : feat.second) {
	if (element < -1 || element > 1) {
	  std::cout << "*** WARNING *** Found element with value outside [-1;1] interval" << std::endl;
	  issue = true;
	  break;
	}
      }
      if (!issue) std::cout << "No issues with input feature detected" << std::endl;
    }
  }
  
};

void ABCNetProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // Get PFCandidate Collection
  edm::Handle<reco::CandidateView> PFCandidates;
  iEvent.getByToken(tokenPFCandidates_, PFCandidates);
  const reco::CandidateView *pfCol = PFCandidates.product();
  
  //make feature map and preprocess features
  std::vector<size_t> indices; 
  auto features = ABCNetMakeInputs::makeFeatureMap(pfCol, indices, debug_);
  ABCNetProducer::preprocess(features, debug_);
  
  //fill the input tensor
  tensorflow::Tensor inputs (tensorflow::DT_FLOAT, { 1, 4000, 19 });
  inputs.flat<float>().setZero();
  for (int j = 0; j < 19; j++) { //may need to find better solution than hard-coding 19
    for (int i = 0; i < 4000; i++) { //may need to find better solution than hard-coding 4000
      inputs.tensor<float,3>()(0,i,j) = float(features[input_names_.at(j)].at(i)); //looks suboptimal; is there way of filling the tensor avoiding nested loops?
    }
  }
  
  std::vector<tensorflow::Tensor> outputs;
  tensorflow::run(session_, { { input_tensor_name_, inputs } }, { output_tensor_name_ }, &outputs);
  //std::cout << "PRINTING NETWORK OUTPUTS" << std::endl;
  //for (int i = 0; i < 4000; i++) std::cout << outputs.at(0).tensor<float,3>()(0, i, 19) << " ";
  
  //initialize container for ABCNet weights
  std::vector<float> weights;
  //throw random numbers in [0,1] as ABCNet weights for now
  //srand(100);
 
  int PFCounter = 0;
  for(auto const& aPF : *pfCol) {
    const pat::PackedCandidate *lPack = dynamic_cast<const pat::PackedCandidate*>(&aPF);
    float abcweight = -1.;
    if(lPack == nullptr) { 
      // throw error
      throw edm::Exception(edm::errors::LogicError,"ABCNetProducer: cannot get weights since inputs are not PackedCandidates");
    }
    else{
      if (indices.at(PFCounter) >= 4000) { //if the particle wasn't considered in evaluation, assign 0 weight
	abcweight = 0.0;
      }
      else abcweight = outputs.at(0).tensor<float,3>()(0, indices.at(PFCounter), 19);
    }
    weights.push_back(abcweight);
    PFCounter++;
  }

  //std::unique_ptr<edm::ValueMap<float>> ABCNetOut(new edm::ValueMap<float>());
  edm::ValueMap<float> ABCNetOut;
  //edm::ValueMap<float>::Filler  ABCNetFiller(*ABCNetOut);
  edm::ValueMap<float>::Filler ABCNetFiller(ABCNetOut);
  ABCNetFiller.insert(PFCandidates,weights.begin(),weights.end());
  ABCNetFiller.fill();
  iEvent.emplace(ABCNetOut_, ABCNetOut);
  //iEvent.put(std::move(ABCNetOut), "weights");

}
//define this as a plug-in
DEFINE_FWK_MODULE(ABCNetProducer);
