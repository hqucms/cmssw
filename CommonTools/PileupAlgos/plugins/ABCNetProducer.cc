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
    VarInfo(float upper_bound, float lower_bound, float norm_factor, float replace_nan_value)
      : upper_bound(upper_bound),
        lower_bound(lower_bound),
	norm_factor(norm_factor),
	replace_nan_value(replace_nan_value) {}
    float upper_bound = 100000;
    float lower_bound = -100000;
    float norm_factor = 1;
    float replace_nan_value = 0;
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
  typedef math::XYZTLorentzVector LorentzVector;
  typedef std::vector< pat::PackedCandidate > PackedOutputCollection;
  explicit ABCNetProducer(const edm::ParameterSet&, const ABCNetTFCache*);   
  ~ABCNetProducer();

  static void fillDescriptions(edm::ConfigurationDescriptions &);

  static std::unique_ptr<ABCNetTFCache> initializeGlobalCache(const edm::ParameterSet&);
  static void globalEndJob(const ABCNetTFCache*);

private:
  std::unique_ptr< PackedOutputCollection > fPackedPuppiCandidates;
  void produce(edm::Event &, const edm::EventSetup &) override;
  // tokens
  edm::EDGetTokenT<reco::CandidateView> tokenPFCandidates_;
  std::vector<std::string> input_names_; //names of the input features. Ordering matters!
  std::unordered_map<std::string, PreprocessParams> prep_info_map_;  //preprocessing info for each input feature
  edm::EDPutTokenT<edm::ValueMap<float>> ABCNetOut_;
};

//////////////////////////////////////////////////////////////////////////////////
// Defining methods --------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////////

// constructors
ABCNetProducer::ABCNetProducer(const edm::ParameterSet& iConfig, const ABCNetTFCache* cache):
  tokenPFCandidates_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("candName")))
{
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
    auto &prep_params = prep_info_map_[input_name];
    prep_params.var_info_map[input_name] = PreprocessParams::VarInfo(upper_bound, lower_bound, norm_factor, replace_nan_value);
  }
  
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

void ABCNetProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // Get PFCandidate Collection
  edm::Handle<reco::CandidateView> PFCandidates;
  iEvent.getByToken(tokenPFCandidates_, PFCandidates);
  const reco::CandidateView *pfCol = PFCandidates.product();
  auto features = ABCNetMakeInputs::makeFeatureMap(pfCol, false);

  //initialize container for ABCNet weights
  std::vector<float> weights;
  //throw random numbers in [0,1] as ABCNet weights for now
  srand(100);
 
  for(auto const& aPF : *pfCol) {  
    const pat::PackedCandidate *lPack = dynamic_cast<const pat::PackedCandidate*>(&aPF);
    float abcweight = -1.;
    if(lPack == nullptr) { 
      // throw error
      throw edm::Exception(edm::errors::LogicError,"ABCNetProducer: cannot get weights since inputs are not PackedCandidates");
    }
    else{
      abcweight = (float) rand()/RAND_MAX;
    }
    weights.push_back(abcweight);
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
