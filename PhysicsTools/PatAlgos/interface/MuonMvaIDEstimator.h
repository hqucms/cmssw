#ifndef __PhysicsTools_PatAlgos_MuonMvaIDEstimator__
#define __PhysicsTools_PatAlgos_MuonMvaIDEstimator__

#include <memory>
#include <string>
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/stream/EDAnalyzer.h"
//

namespace pat {
  class Muon;
}

namespace edm {
  class FileInPath;
}

namespace pat {
class MuonMvaIDEstimator : public edm::stream::EDProducer<edm::GlobalCache<cms::Ort::ONNXRuntime>> {
//class MuonMvaIDEstimator {
	public:
           explicit MuonMvaIDEstimator(const edm::ParameterSet&, const cms::Ort::ONNXRuntime *);
	   ~MuonMvaIDEstimator();
	 
	   static void fillDescriptions(edm::ConfigurationDescriptions &);
	   static std::unique_ptr<cms::Ort::ONNXRuntime> initializeGlobalCache(const edm::ParameterSet &);
	   static void globalEndJob(const cms::Ort::ONNXRuntime *);
	   std::vector<float> computeMVAID(const pat::Muon& imuon) const;
	   
	private:
	   const pat::Muon& imuon;
	   std::vector<std::string> flav_names_;             // names of the output scores
	   std::vector<std::string> input_names_;            // names of each input group - the ordering is important!
	   std::vector<float> outputs; 
           //std::unique_ptr<const cms::Ort::ONNXRuntime> outputs;
	 };
};	 
#endif
