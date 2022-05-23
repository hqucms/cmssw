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
//         Created:  THU, 5 May 2022 11:00:28 GMT  ---------------------------------------
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

#include <limits>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <memory>
#include <nlohmann/json.hpp>

///////////////////////////////////////////////////////////////////////////////////
// Declaring class ----------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////////

class ABCNetProducer : public edm::stream::EDProducer<> {

public:
  typedef math::XYZTLorentzVector LorentzVector;
  typedef std::vector< pat::PackedCandidate > PackedOutputCollection;
  explicit ABCNetProducer(const edm::ParameterSet&);   
  ~ABCNetProducer();

  static void fillDescriptions(edm::ConfigurationDescriptions &);

private:
  std::unique_ptr< PackedOutputCollection > fPackedPuppiCandidates;
  void produce(edm::Event &, const edm::EventSetup &) override;
  // tokens
  edm::EDGetTokenT<reco::CandidateView> tokenPFCandidates_;
};

//////////////////////////////////////////////////////////////////////////////////
// Defining methods --------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////////

// constructors
ABCNetProducer::ABCNetProducer(const edm::ParameterSet& iConfig):
  tokenPFCandidates_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("candName")))
{
  // Produce a ValueMap of floats linking each PF candidate with its ABCNet weight
  produces<edm::ValueMap<float> > ();

};

// destructor
ABCNetProducer::~ABCNetProducer() {
};

void ABCNetProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {};

void ABCNetProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  std::cout << "I AM ABCNET PRODUCER, I AM PRODUCING" << std::endl;
  // Get PFCandidate Collection
  edm::Handle<reco::CandidateView> PFCandidates;
  iEvent.getByToken(tokenPFCandidates_, PFCandidates);
  const reco::CandidateView *pfCol = PFCandidates.product();
  
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
  
  std::unique_ptr<edm::ValueMap<float>> ABCNetOut(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler  ABCNetFiller(*ABCNetOut);
  ABCNetFiller.insert(PFCandidates,weights.begin(),weights.end());
  ABCNetFiller.fill();
  
  iEvent.put(std::move(ABCNetOut), "weights");
  
}
//define this as a plug-in
DEFINE_FWK_MODULE(ABCNetProducer);
