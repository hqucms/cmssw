#include "PhysicsTools/ONNXRuntime/ONNXRuntime.h"
#include "PhysicsTools/PatAlgos/interface/MuonMvaIDEstimator.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "CommonTools/MVAUtils/interface/GBRForestTools.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"

using namespace pat;
using namespace cms::Ort;

MuonMvaIDEstimator::~MuonMvaIDEstimator() {}
 
 void MuonMvaIDEstimator::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
   // pfDeepBoostedJetTags
   edm::ParameterSetDescription desc;
   desc.add<edm::FileInPath>("mvaIDTrainingFile",
                             edm::FileInPath("RecoMuon/MuonIdentification/data/mvaID.onnx"));
   desc.add<std::vector<std::string>>("flav_names",
                                      std::vector<std::string>{
                                          "probGOOD",
                                          "probBAD",
                                      });
 
   descriptions.addWithDefaultLabel(desc);
 }
 
 std::unique_ptr<cms::Ort::ONNXRuntime> MuonMvaIDEstimator::initializeGlobalCache(const edm::ParameterSet &iConfig) {
   return std::make_unique<cms::Ort::ONNXRuntime>(iConfig.getParameter<edm::FileInPath>("mvaIDTrainingFile").fullPath());
 }
 void MuonMvaIDEstimator::globalEndJob(const cms::Ort::ONNXRuntime *cache) {}
 
 float MuonMvaIDEstimator::computeMVAID(const pat::Muon& muon) const {
   float local_chi2  = muon.combinedQuality().chi2LocalPosition;
   float kink  = muon.combinedQuality().trkKink;   
   float segment_comp =  muon.segmentCompatibility( arbitrationType);  
   float nMatchedStations = muon.numberOfMatchedStations();   
   float pt = muon.pt();
   float eta = muon.eta();
   float global_muon = muon.isGlobalMuon();  
   if (muon.innerTrack().isNonnull()){
       float Valid_pixel  = muon.innerTrack()->hitPattern().numberOfValidPixelHits();
       float tracker_layers  = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement();
       float validFraction   = muon.innerTrack()->validFraction();
       }
   else{
       Valid_pixel = -99.;
       tracker_layers = -99;
       validFraction   = -99;
       }  
   if (muon.globalTrack().isNonnull()){
       float norm_chi2  = muon.globalTrack()->normalizedChi2();
       float n_Valid_hits = muon.globalTrack()->hitPattern().numberOfValidMuonHits(); 
       }
   else{
       float norm_chi2  = muon.innerTrack()->normalizedChi2();
       float n_Valid_hits = muon.innerTrack()->hitPattern().numberOfValidMuonHits();
       }
       
   std::vector<std::string> input_names_ {'global_muon','validFraction','norm_chi2','local_chi2','kink','segment_comp','n_Valid_hits','n_MatchedStations','Valid_pixel','tracker_layers','pt','eta'};
   float input_values_ {global_muon,validFraction,norm_chi2,local_chi2,kink,segment_comp,n_Valid_hits,n_MatchedStations,Valid_pixel,tracker_layers,pt,eta};
   std::vector<std::string> flav_names_{"probGOOD","probBAD"}

   std::vector<float> outputs(flav_names_.size(), 0);  // init as all zeros

   outputs = globalCache()->run(input_names_, input_values_)[0];
   assert(outputs.size() == flav_names_.size());
   return outputs
}
