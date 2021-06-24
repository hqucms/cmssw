#ifndef __PhysicsTools_PatAlgos_MuonMvaIDEstimator__
#define __PhysicsTools_PatAlgos_MuonMvaIDEstimator__

#include <memory>
#include <string>

namespace pat {
  class Muon;
}

namespace edm {
  class FileInPath;
}

namespace pat {
class MuonMvaIDEstimator : public edm::stream::EDProducer< edm::GlobalCache<ONNXRuntime> > {
	public:
	   explicit MuonMvaIDEstimator(const edm::ParameterSet&, const ONNXRuntime *);
	   ~MuonMvaIDEstimator();
	 
	   static void fillDescriptions(edm::ConfigurationDescriptions &);
	   static std::unique_ptr<ONNXRuntime> initializeGlobalCache(const edm::ParameterSet &);
	   static void globalEndJob(const ONNXRuntime *);
	   float produce(const pat::Muon& imuon) const;
	   
	private:
	   const pat::Muon& imuon;
	   std::vector<std::string> flav_names_;             // names of the output scores
	   std::vector<std::string> input_names_;            // names of each input group - the ordering is important!
	   FloatArrays input_values_;
	 };
};	 
#endif
